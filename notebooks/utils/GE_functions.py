import pandas as pd 
import numpy as np 
import plotly.express as px
import plotly.graph_objects as go
from matplotlib import cm, colors as mcolors
from sklearn.metrics import mean_squared_error
from sklearn.decomposition import PCA
import copy 
import random

def get_ge_table(adata, gene_names, ge_across = 'leiden', sort_by = 'Rtp1', dpt = True):
# construct a pd dataframe for average gene expression of interest across dpt_bins / selected_dpt
    ge_table = pd.DataFrame(columns=gene_names, index = adata.obs[ge_across].unique()).sort_index()
    if dpt:
        dpt_average = pd.DataFrame(columns = ['dpt_average'], index = adata.obs[ge_across].unique()).sort_index()
    for clust in ge_table.index: 
    #     following line subsets the adata by leiden AND subset from genes with gene_names. Then find mean across index(leiden)
        ge_table.loc[clust] = adata[adata.obs[ge_across].isin([clust]),:][:,gene_names].X.mean(0) #mean(0) performs average by the columns
        if dpt:
            dpt_average.loc[clust] = adata.obs[adata.obs[ge_across] == clust].dpt_pseudotime.mean()
    if dpt:
        ge_table = ge_table.join(dpt_average)
    # drop genes with 0 expression/counts 
    ge_table = ge_table.loc[:, (ge_table != 0).any(axis=0)]
    # assigns index as an column 
    ge_table.reset_index(inplace=True)
    # sort the rows by ascending value of Rtp1 expression 
    if not sort_by == None: 
        ge_table.sort_values(by=sort_by, inplace=True)
    return ge_table

# quick plotting function 
def expression_line_plot(ge_table, mse_table):
    ge_data = get_ge_data(ge_table, mse_table)
    # plotting with plolty express
    fig = px.line(ge_data, x='index', y= 'expression', color='gene')
    fig.update_layout(xaxis_type = 'category') # update x-axis to category so that it doesn't sort the numbers
    return fig

def GO_expression_line_plot(ge_table, mse_table, opacity = 0.3, style = 'dot', width = 5):
    ge_data = get_ge_data(ge_table, mse_table)
    fig = go.Figure()
    for g in ge_data['gene'].unique():
        fig.add_trace(go.Scatter(x= ge_data[ge_data['gene'] == g]['index'], 
                                 y= ge_data[ge_data['gene'] == g]['expression'],
                                opacity = opacity,
                                mode ='lines',
                                name = g,
                                line=dict(
#                                     color='grey',
                                    width=width,
                                    dash = style
                                )
                            ))
    return fig


def get_ge_data(ge_table, genes):
    temp_genes = genes.copy()
    temp_genes.append('index')
        
    ge_data = pd.melt(ge_table[temp_genes], 
                      id_vars='index', value_vars=ge_table[temp_genes],
                      var_name='gene', value_name='expression')
    print('ge_data constructed')
    return ge_data

# Deprecated as calculating across numpy is faster than iterating columns
# def get_mse_table(ge_table, gene_names, mse_gene = 'Rtp1'):
#     gene_names = np.intersect1d(ge_table.columns, gene_names)
#     mse_table = pd.DataFrame(columns=gene_names, index=['mse'])  

#     # Calculate mean squared error of dpt_average in respect to Rtp1
#     for gene in gene_names:
#         mse_table[gene] = mean_squared_error(ge_table[gene],ge_table[mse_gene]) 
#     print('mse_table constructed')
#     return mse_table

def get_rmse_table(ge_table, subject_gene = 'Rtp1', no_filter=False):
    
    if not no_filter: 
        genes = [col for col in ge_table.columns if not any(exc_str in col for exc_str in ['dpt_average', '_sum'])]
    else: 
        genes = ge_table.columns
    # Calculating mean_square_error for all the columns in a dictionary form 
    mse_values = (ge_table[genes].sub(ge_table[subject_gene], axis=0)**2).mean()
    return mse_values

def filter_mse_table(mse_table, excluded_genes = [], mse_gene = 'Rtp2'):
    drop_gene = []
    # Filter for genes that have smaller mse than Rtp2 - Rtp1. These are defined as associated to Rtp1
    for gene in mse_table.columns:
        if not gene in excluded_genes:
            if (mse_table[gene] > mse_table[mse_gene]).bool():
    #         if (mse_table[gene] > 0.1).bool():
                drop_gene += [gene]
    mse_table = mse_table.drop(columns=drop_gene)
    print(len(drop_gene), 'genes dropped from mse_table')
    return mse_table

def get_gene_names(adata):
    """
    construct a gene_names of list of genes to include analysis
    exclude mt- , Rik ... etc 
    """
    rm_genes = adata.var_names.str.startswith('mt-')
    rm_genes = np.add(rm_genes, adata.var_names.str.startswith('n-'))
#     rm_genes = np.add(rm_genes, adata.var_names.str.startswith('Rik'))
#     rm_genes = np.add(rm_genes, adata.var_names.str.startswith('Olfr'))
#     rm_genes = np.add(rm_genes, adata.var_names.str.startswith('Gm'))
    keep = np.invert(rm_genes)
    gene_names = adata[:,keep].var_names
    # Filters out the genes with double capitalization unknown genes 
    double_upper = []
    for gene in gene_names: 
        if len(gene) > 1:
            if gene[1].isupper():
                double_upper.append(gene)
    gene_names = [i for i in gene_names if i not in double_upper]
    return gene_names

def get_excluded_genes(gene_names = None):
    """
    list of genes to exclude exemption from being eliminated 
    Finds the manually exclude genes in the gene_names list 
    These are mature OSN markers 
    """
    excluded_genes = ['Omp','Cnga4','Ano2', 'Cnga2', 'Cngb1',
                     'Stoml3', 'Stom', 'Rtp1', 'Rtp2', 'Clgn', 'Gfy', 'Gnal',
                     'Gng13','Adcy3','Ric8b']
    if gene_names != None:
        excluded_genes = [x for x in gene_names if x in excluded_genes]

    return excluded_genes

def add_dpt_bins(adata):
    # Creates a column of dpt_bins to group cells by increments of 0.1 across 0-1
    adata.obs['dpt_bins'] = None 
    for dpt_cutoff in np.linspace(1,0, num = 11):
        adata.obs.loc[adata.obs.dpt_pseudotime < dpt_cutoff, 'dpt_bins'] = str(abs(round(
            dpt_cutoff-0.1,1))) +"-"+ str(round(dpt_cutoff,1))
    return adata
    
def add_selected_dpt(adata, dpt_cutoff = [0,0.5,0.9,0.98,1]):
    """
    Creates a column of dpt_bins to group cells by dpt_bins 
    dpt_bins step size can be tuned by adjusting the num slices across 1-0 in np.linespace. 
    """
    adata.obs['selected_dpt'] = None 
    for cutoff in enumerate(dpt_cutoff):
        if not cutoff[1] == 1:
            adata.obs.loc[adata.obs.dpt_pseudotime >= cutoff[1], 'selected_dpt'] = str((cutoff[1])) +"-"+ str(dpt_cutoff[cutoff[0]+1])
    return adata


def x_y_scatter(mse_x_y, title="", label_genes=['Rtp1','Rtp2']):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x = mse_x_y[mse_x_y.columns[0]], 
                             y = mse_x_y[mse_x_y.columns[1]], 
                             mode = 'markers', 
                             hovertemplate = mse_x_y.index))
    fig.add_trace(go.Scatter(x = mse_x_y.loc[label_genes][mse_x_y.columns[0]], 
                             y = mse_x_y.loc[label_genes][mse_x_y.columns[1]],
                             mode = 'markers',
                             hovertemplate = mse_x_y.loc[label_genes].index,
                             marker = dict(size = 10)))
    fig.update_layout(
        title=title,
        xaxis_title=mse_x_y.columns[0],
        yaxis_title=mse_x_y.columns[1],
        showlegend=False,
        width = 800,
        height = 800,
        plot_bgcolor = "white"
    )
    fig.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1,
      )
    return fig 

def rank_plot(rank_x_y, title="", excluded_genes = ['Rtp1', 'Rtp2']):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x = rank_x_y[rank_x_y.columns[0]], 
                             y = rank_x_y[rank_x_y.columns[1]], 
                             mode = 'markers', 
                             marker = dict(size = 4, opacity = 0.5),
                             hovertemplate = rank_x_y.index))
    fig.add_trace(go.Scatter(x = rank_x_y.loc[excluded_genes][rank_x_y.columns[0]], 
                             y = rank_x_y.loc[excluded_genes][rank_x_y.columns[1]],
                             mode = 'markers',
                             marker = dict(size = 8, color = 'pink'),
                            hovertemplate = rank_x_y.loc[excluded_genes].index))
    fig.update_layout(
        title=title,
        xaxis_title=rank_x_y.columns[0],
        yaxis_title=rank_x_y.columns[1],
        showlegend=False,
        width = 800,
        height = 800,
        plot_bgcolor = "white"
    )
    return fig

# transforms cell_top_markers table into ranks matching the rank_genes_table
def gene_to_rank(cell_top_markers, gene_names, reverse_rank = True):
    cell_genes_rank = pd.DataFrame(index = gene_names)
    for column in cell_top_markers.columns:
        temp = cell_top_markers[[column]].set_index(column)
        temp = temp[temp.index.isin(cell_genes_rank.index)]
        if reverse_rank:
            temp[column] = range(len(temp),0,-1) 
        else: 
            temp[column] = range(0,len(temp))
        cell_genes_rank = cell_genes_rank.join(temp)
    return cell_genes_rank

def gene_to_relative_expression(ge_table):
    # Normalizes gene expression to it's maximum expression value 
    for gene in ge_table.columns:
        if gene not in ['index', 'dpt_average']:
            ge_table[gene] = ge_table[gene]/max(ge_table[gene])
    return ge_table

def mse_to_rank(mse_combined, col_names):
    rank_index = []
    for i in range(1,len(mse_combined)+1):
        rank_index.append(i)
        
    normalized = mse_combined[mse_combined.columns[0]]
    log1p = mse_combined[mse_combined.columns[1]]
    
    rank_normalized = pd.DataFrame(index = normalized.sort_values().index)
    rank_normalized[col_names[0]] = rank_index
    rank_log1p = pd.DataFrame(index = log1p.sort_values().index)
    rank_log1p[col_names[1]] = rank_index
    rank_combined = rank_normalized.join(rank_log1p)
    return rank_combined

def get_pca_df(data, dimensions=2):
    pca = PCA(n_components=dimensions)
    principalComponents = pca.fit_transform(data)
    columns = ['PC_%i' % i for i in range(dimensions)]
    PCA_df = pd.DataFrame(data = principalComponents
                 , columns = columns)
    PCA_df = PCA_df.set_index(data.index)
    return PCA_df


def get_color(colorscale_name, loc):
    """
# This function allows you to retrieve colors from a continuous color scale
# by providing the name of the color scale, and the normalized location between 0 and 1
# Reference: https://stackoverflow.com/questions/62710057/access-color-from-plotly-color-scale
    """


    from _plotly_utils.basevalidators import ColorscaleValidator
    # first parameter: Name of the property being validated
    # second parameter: a string, doesn't really matter in our use case
    cv = ColorscaleValidator("colorscale", "")
    # colorscale will be a list of lists: [[loc1, "rgb1"], [loc2, "rgb2"], ...] 
    colorscale = cv.validate_coerce(colorscale_name)
    
    if hasattr(loc, "__iter__"):
        return [get_continuous_color(colorscale, x) for x in loc]
    return get_continuous_color(colorscale, loc)
        

# Identical to Adam's answer
import plotly.colors
from PIL import ImageColor

def get_continuous_color(colorscale, intermed):
    """
    Plotly continuous colorscales assign colors to the range [0, 1]. This function computes the intermediate
    color for any value in that range.

    Plotly doesn't make the colorscales directly accessible in a common format.
    Some are ready to use:
    
        colorscale = plotly.colors.PLOTLY_SCALES["Greens"]

    Others are just swatches that need to be constructed into a colorscale:

        viridis_colors, scale = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Viridis)
        colorscale = plotly.colors.make_colorscale(viridis_colors, scale=scale)

    :param colorscale: A plotly continuous colorscale defined with RGB string colors.
    :param intermed: value in the range [0, 1]
    :return: color in rgb string format
    :rtype: str
    """
    if len(colorscale) < 1:
        raise ValueError("colorscale must have at least one color")

    hex_to_rgb = lambda c: "rgb" + str(ImageColor.getcolor(c, "RGB"))

    if intermed <= 0 or len(colorscale) == 1:
        c = colorscale[0][1]
        return c if c[0] != "#" else hex_to_rgb(c)
    if intermed >= 1:
        c = colorscale[-1][1]
        return c if c[0] != "#" else hex_to_rgb(c)

    for cutoff, color in colorscale:
        if intermed > cutoff:
            low_cutoff, low_color = cutoff, color
        else:
            high_cutoff, high_color = cutoff, color
            break

    if (low_color[0] == "#") or (high_color[0] == "#"):
        # some color scale names (such as cividis) returns:
        # [[loc1, "hex1"], [loc2, "hex2"], ...]
        low_color = hex_to_rgb(low_color)
        high_color = hex_to_rgb(high_color)

    return plotly.colors.find_intermediate_color(
        lowcolor=low_color,
        highcolor=high_color,
        intermed=((intermed - low_cutoff) / (high_cutoff - low_cutoff)),
        colortype="rgb",
    )

def normalize_list(lst):
    max_val = max(lst)
    min_val = min(lst)
    normalized_lst = [(x - min_val) / (max_val - min_val) for x in lst]
    return normalized_lst

def reverse_scale(values, factor = 1):
    """
    Reverses the scale of a list of values such that the smallest value becomes 1 and the largest value becomes 0.
    """
    min_val = min(values)
    max_val = max(values)
    scaled_values = [(val - min_val) / (max_val - min_val)*factor for val in values]
    reversed_values = [1 - val for val in scaled_values]
    return reversed_values




def distinct_colors(label_list, category='tab10', custom_color=None, random_state=0):
    """
    Generate distinct colors for a list of labels.

    Parameters:
    label_list (list): A list of labels for which you want to generate distinct colors.
    category (str): Category of distinct colors. Options are 'warm', 'floral', 'rainbow', 'pastel',
                    matplotlib color palettes (e.g., 'tab10', 'Set2'), or None for random. Default is None.
    custom_color (list): A custom list of colors to use.
    random_state (int): Seed for random color generation. Default is 0.

    Returns:
    dict: A dictionary where labels are keys and distinct colors (in hexadecimal format) are values.

    Example:
    >>> labels = ['A', 'B', 'C']
    >>> color_mapping = distinct_colors(labels, category='tab10')
    >>> print(color_mapping)
    {'A': '#1f77b4', 'B': '#ff7f0e', 'C': '#2ca02c'}
    """
    random.seed(random_state)
    
    warm_colors = ['#fabebe', '#ffd8b1', '#fffac8', '#ffe119', '#ff7f00', '#e6194B']
    floral_colors = ['#bfef45', '#fabed4', '#aaffc3', '#ffd8b1', '#dcbeff', '#a9a9a9']
    rainbow_colors = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4']
    pastel_colors = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', 
                     '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928', 
                     '#8DD3C7', '#BEBADA', '#FFED6F']
    
    color_dict = {}

    if custom_color is not None: 
        assert len(custom_color) >= len(label_list), "Provided label_list needs to be shorter than provided custom_color"
        for i, _label in enumerate(label_list): 
            color_dict[_label] = custom_color[i]
        return color_dict

    color_palette = None

    # Handle custom categories
    if category in ['warm', 'floral', 'rainbow', 'pastel']: 
        if category == 'warm':
            color_palette = random.sample(warm_colors, len(warm_colors))
        elif category == 'floral':
            color_palette = random.sample(floral_colors, len(floral_colors))
        elif category == 'rainbow':
            color_palette = random.sample(rainbow_colors, len(rainbow_colors))
        elif category == 'pastel': 
            color_palette = random.sample(pastel_colors, len(pastel_colors))
        # else:
        #     color_palette = random.sample(warm_colors + floral_colors + rainbow_colors + pastel_colors, len(label_list))
    # Handle matplotlib color palettes
    elif category in mcolors.TABLEAU_COLORS or category in cm.cmaps_listed or hasattr(cm, category):
        cmap = cm.get_cmap(category) if hasattr(cm, category) else cm.get_cmap('tab10')
        num_colors = len(label_list)
        color_palette = [mcolors.to_hex(cmap(i / num_colors)) for i in range(num_colors)]
    
    for i, label in enumerate(label_list):
        color_dict[label] = color_palette[i % len(color_palette)]
    
    return color_dict

import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

def continuous_colors(label_list, colormap='viridis', custom_color=None, orders=None):
    """
    Generate continuous colors for a list of labels.

    Parameters:
    label_list (list): A list of labels for which you want to generate continuous colors.
    colormap (str or matplotlib colormap, optional): The colormap to use for color scaling. Default is 'viridis'.
    custom_color (list, optional): A list of color tuples defining the custom colormap.
                                   Default is None.
    orders (list, optional): A list defining the hierarchy of label_list. Default is None.

    Returns:
    dict: A dictionary where labels are keys and continuous colors (in hexadecimal format) are values.

    Example:
    >>> labels = ['A', 'B']
    >>> custom_color = [(0, '#DBE5EB'), (0.5, '#67879B'), (1, '#073763')]
    >>> color_mapping = continuous_colors(labels, custom_color=custom_color)
    >>> print(color_mapping)
    {'A': '#DBE5EB', 'B': '#073763'}
    """
    color_dict = {}

    # Choose colormap
    if isinstance(colormap, str):
        cmap = cm.get_cmap(colormap)
    else:
        cmap = colormap

    # Generate custom colormap
    if custom_color is not None:
        custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', custom_color)
    else:
        custom_cmap = None

    # Generate colors
    num_labels = len(label_list)
    for i, label in enumerate(label_list):
        if custom_cmap is not None:
            norm_color = i / (num_labels - 1)  # Normalize color index
            color = cm.colors.rgb2hex(custom_cmap(norm_color))
        else:
            color = cm.colors.rgb2hex(cmap(i / (num_labels - 1)))  # Normalize color index
        color_dict[label] = color

    # Reorder color_dict based on orders if provided
    if orders is not None:
        color_dict = {label: color_dict[label] for label in orders if label in color_dict}

    return color_dict