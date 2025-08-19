
import pandas as pd
import numpy as np 

import plotly.express as px 
import plotly.graph_objects as go

from sklearn.decomposition import PCA
from scipy import stats


def compare_vol_plot(DE_df_list: list, 
                     DE_df_name: list, 
                     fig_title = '', 
                     fig_dimension = None, 
                     fig_fixed_range = False):
    
    assert len(DE_df_list) == len(DE_df_name), print("Please ensure DE_df_list and DE_df_name len matches")
    
    # Initialize ranges for plots to prevent autosizing
    xmin, xmax, ymin, ymax = 0, 0, 0, 0
    
    fig = go.Figure()
    # Add traces of individual DE_df 
    for DE_df, DE_name in zip(DE_df_list, DE_df_name): 
        plot_df = DE_df.copy()
        fig.add_trace(go.Scatter(x = plot_df['logFC'], 
                                y = -np.log10(plot_df['FDR']),
                                text = plot_df['symbol'],
                                mode = 'markers', 
                                name = DE_name,
                                marker = dict(size = 10, 
                                            #    color = 'grey', 
                                            opacity=0.3)
                                )
                    )
        if fig_fixed_range:
            xmin =  min(plot_df['logFC'])*1.10 if min(plot_df['logFC']) < xmin else xmin
            xmax =  max(plot_df['logFC'])*1.10 if max(plot_df['logFC']) > xmax else xmax
            # ymin =  min(-np.log10(plot_df['FDR']))*1.10 if min(-np.log10(plot_df['FDR'])) < ymin else ymin  # ymin will be 0 anyways
            ymax =  max(-np.log10(plot_df['FDR']))*1.10 if max(-np.log10(plot_df['FDR'])) > ymax else ymax


    # Add a line for FDR = 0.05
    fig.add_shape(type='line', x0=-10, x1=10,
                  y0=-np.log10(0.05), y1=-np.log10(0.05),
                  line=dict(color='violet', width=3, dash='dash'))

    fig.update_traces( 
        textposition='top center',
        hovertemplate =
        '<b>%{text}</b>' + 
        '<br>LogFC: %{x}'+
        '<br>FDR: %{y}<br>')
    
    # Define ranges to avoid autoresizing when hiding data
    if fig_fixed_range:
        # Center the data by taking the bigger value between xmin and xmax 
        xmax = max(abs(xmin),abs(xmax))
        fig.update_xaxes(range=[-xmax, xmax])
        fig.update_yaxes(range=[-1, ymax])

    if fig_dimension is None:
        fig.update_layout(
            title=fig_title,
            autosize=True,
            template='simple_white'
        )
    else: 
        fig.update_layout(
            title=fig_title,
            width=fig_dimension[0],
            height=fig_dimension[1],
            template='simple_white'
        )
    
    
    return fig 


def reduced_dimension_plot(count_df: pd.DataFrame(), 
                           reduction_method = 'umap',
                           pca_n_components = 2):
    
    import umap
    
    assert (reduction_method == 'umap') | (reduction_method == 'pca'), print('reduction_method needs to be either \'pca\' or \'umap\'')
    
    if reduction_method.lower() == 'pca':
        pca = PCA(n_components=pca_n_components)  # You can adjust the number of components
        result = pca.fit_transform(count_df.transpose())
        col = [f"pca_{i}"  for i in range(1,pca_n_components+1)]
    elif reduction_method.lower() == 'umap':
        result = umap.UMAP().fit_transform(count_df.transpose())
        col = ['umap_1', 'umap_2']
    
    result = pd.DataFrame(result, columns = col)
    result['sample_name'] = count_df.columns
    result['sample'] = result['sample_name'].str.split('_').str[0]
    result['odor'] = result['sample_name'].str.split('_').str[1]

    fig = px.scatter(result, 
                    x = col[0], 
                    y = col[1], 
                    hover_name = 'sample_name', 
                    color = 'sample', 
                    symbol = 'odor')
    
    fig.update_traces(marker={'size': 15})
    
    fig.update_layout(
        title=reduction_method,
        autosize=True,
        template='simple_white'
    )
    
    return fig
def add_p_value_annotation(fig, 
                           array_columns, 
                           just_annotate = None,
                           test_type = 'ranksums', 
                           popmean = None, 
                           y_padding = True, 
                           subplot=None, 
                           include_tstat=None, 
                           p_round=3,
                           font_size=14, 
                           show=None, 
                           select_datatype=None,
                           _format=dict(interline=0.07, text_height=1.07, color='black')):
    ''' Adds notations giving the p-value between two box plot data (t-test two-sided comparison)
    
    Parameters:
    ----------
    fig: figure
        plotly boxplot figure
    array_columns: np.array
        array of which columns to compare 
        e.g.: [[0,1], [1,2]] compares column 0 with 1 and 1 with 2
    subplot: None or int
        specifies if the figures has subplots and what subplot to add the notation to
    _format: dict
        format characteristics for the lines

    Returns:
    -------
    fig: figure
        figure with the added notation
    '''
    
    assert test_type in ['ranksums', 'ttest_ind', 'ttest_rel', 'ttest_1samp'], \
        "Please specify test_type to be either ranksums or ttest"
    if test_type == 'ttest_1samp': 
        assert popmean is not None, "ttest_1samp requires popmean value"
    
    if just_annotate is not None: 
        assert len(just_annotate) == len(array_columns), "'just_annotate' and 'array_columns' len must be identical "


    # Filter data for subplots
    if subplot:
        subplot_str = '' if subplot == 1 else str(subplot)
        selected_data = [
            (i, data) for i, data in enumerate(fig.data) 
            if data['xaxis'] == 'x' + subplot_str
        ]
    else:
        subplot_str = ''
        selected_data = list(enumerate(fig.data))

    # Filter by datatype if specified
    if select_datatype:
        selected_data = [(i, d) for i, d in selected_data if d['type'] == select_datatype]
    
    # Extract only y-data and indices after filtering
    filtered_indices = [i for i, _ in selected_data]
    filtered_y_data = [d['y'] for _, d in selected_data]

    # Prepare annotation y-positions
    y_range = np.zeros([len(array_columns), 2])
    for i in range(len(array_columns)):
        base = 1.01 + (i * _format['interline'] if y_padding else _format['interline'])
        y_range[i] = [base, base + 0.01]

    # Main loop for annotation
    for idx, column_pair in enumerate(array_columns):
        idx0, idx1 = column_pair
        y0 = filtered_y_data[idx0]
        y1 = filtered_y_data[idx1]
        
        if test_type == 'ttest_ind':
            tstat, pvalue = stats.ttest_ind(y0, y1, equal_var=False)
        elif test_type == 'ttest_rel':
            tstat, pvalue = stats.ttest_rel(y0, y1)
        elif test_type == 'ranksums':
            tstat, pvalue = stats.ranksums(y0, y1)
        elif test_type == 'ttest_1samp':
            tstat, pvalue = stats.ttest_1samp(y0, popmean=popmean)

        # Symbol formatting
        symbol = just_annotate[idx] if just_annotate else format_pvalue(pvalue, p_round, tstat if include_tstat else None, show)

        # Draw lines only if different
        if idx0 != idx1:
            for x in [idx0, idx1]:
                fig.add_shape(type="line",
                    xref="x"+subplot_str, yref="y"+subplot_str+" domain",
                    x0=x, y0=y_range[idx][0]*_format['text_height'],
                    x1=x, y1=y_range[idx][1]*_format['text_height'],
                    line=dict(color=_format['color'], width=2)
                )
            fig.add_shape(type="line",
                xref="x"+subplot_str, yref="y"+subplot_str+" domain",
                x0=idx0, y0=y_range[idx][1]*_format['text_height'],
                x1=idx1, y1=y_range[idx][1]*_format['text_height'],
                line=dict(color=_format['color'], width=2)
            )

        fig.add_annotation(dict(
            font=dict(color=_format['color'], size=font_size),
            x=(idx0 + idx1)/2,
            y=y_range[idx][1]*(_format['text_height']+.07),
            showarrow=False,
            text=symbol,
            textangle=0,
            xref="x"+subplot_str,
            yref="y"+subplot_str+" domain"
        ))
    return fig

def format_pvalue(pvalue, p_round=3, t=None, show=None):
    """
    Format a p-value as a string with significance symbols.

    Parameters:
    pvalue (float): The p-value to be formatted.
    p_round (int): Number of digits to show in scientific notation.
    t (float, optional): Optional t-statistic to include.
    show (str or None): 
        If 'symbol', returns only the significance label (e.g. '*').
        If 'pvalue', returns only the formatted p-value (e.g. 'p=1.23e-03').
        If None, returns full string with significance and p-value (and t if provided).

    Returns:
    str: The formatted p-value string.
    """
    
    # Choose significance symbol
    if np.isnan(pvalue): 
        symbol_part= ''
    elif pvalue >= 0.05:
        symbol_part = 'ns'
    elif pvalue >= 0.01:
        symbol_part = '*'
    elif pvalue >= 0.001:
        symbol_part = '**'
    else:
        symbol_part = '***'
    
    # Format p-value using scientific notation
    pval_part = f'p={pvalue:.{p_round}e}'
    
    # Handle return based on return_part
    if show == 'symbol':
        return symbol_part
    elif show == 'pvalue':
        return pval_part
    
    # Default: full output
    symbol = f'{symbol_part} <br>{pval_part}'
    if t is not None:
        symbol += f'<br>t={round(t, 3)}'
    
    return symbol



from scipy.spatial import distance

def umap_euclidean_distance(umap_df, 
                            by,
                            include_shuffle = True, 
                            shuffled_fraction = 0.3):
    # Calculate pairwise distances
    distances = []

    groups = list(umap_df[by].sort_values().unique())
    if include_shuffle:
        groups += ['shuffled']
    
    for group in groups:
        if (group == 'shuffled'): 
            unique_top_Olfr = umap_df.sample(frac = shuffled_fraction)['top_Olfr'].unique()
        else: 
            unique_top_Olfr = umap_df[umap_df[by] == group]['top_Olfr'].unique()
        
        for olfr in unique_top_Olfr:
            olfr_data = umap_df[umap_df['top_Olfr'] == olfr]
            open_coords = olfr_data[olfr_data['nostril'] == 'open'][['umap_x', 'umap_y']].values
            close_coords = olfr_data[olfr_data['nostril'] == 'close'][['umap_x', 'umap_y']].values

            for open_point in open_coords:
                for close_point in close_coords:
                    dist = distance.euclidean(open_point, close_point)
                    distances.append({'top_Olfr': olfr, 'distance': dist, 'group': group})

    # Create a DataFrame with pairwise distances
    return pd.DataFrame(distances)