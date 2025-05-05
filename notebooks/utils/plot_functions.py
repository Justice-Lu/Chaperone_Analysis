import numpy as np 
import pandas as pd 
import plotly.graph_objects as go 
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, linregress

import utils.color_function as cf 

# import logomaker

def _plotly_blank_style(fig): 
    """
    Simply update the plotly go figures to transparent background for better 3D visualization
    """
    
    fig.update_layout(
        scene=dict(
                xaxis=dict(visible=False, showbackground=False),
                yaxis=dict(visible=False, showbackground=False),
                zaxis=dict(visible=False, showbackground=False)
                ),
        margin=dict(r=10, l=10, b=10, t=10)
    )
    return fig

def _plotly_fixed_axes_ranges(fig: go.Figure, padding=0, percentage_padding=0.1, fixedrange=False) -> go.Figure:
    """
    Fix the axis ranges of a Plotly go.Figure so that hiding traces 
    does not resize the plot axes, with options for absolute or percentage-based padding.

    Parameters:
    - fig (go.Figure): A Plotly figure.
    - padding (float): Absolute padding to add/subtract from min/max axis limits.
    - percentage_padding (float): Percentage padding applied to min/max as a multiplier.

    Returns:
    - go.Figure: The modified figure with fixed axis ranges.
    """
    # Get all x and y data from traces
    x_data, y_data = [], []

    for trace in fig.data:
        if 'x' in trace and trace.x is not None:
            x_data.extend(trace.x)
        if 'y' in trace and trace.y is not None:
            y_data.extend(trace.y)

    # Determine the axis ranges
    if x_data:
        x_min, x_max = min(x_data), max(x_data)
        x_range_span = x_max - x_min  # Total span of x values
        x_padding = (x_range_span * percentage_padding) / 2  # Apply percentage padding to both sides
        x_range = [x_min - padding - x_padding, x_max + padding + x_padding]
    else:
        x_range = None

    if y_data:
        y_min, y_max = min(y_data), max(y_data)
        y_range_span = y_max - y_min  # Total span of y values
        y_padding = (y_range_span * percentage_padding) / 2  # Apply percentage padding to both sides
        y_range = [y_min - padding - y_padding, y_max + padding + y_padding]
    else:
        y_range = None

    # Update figure layout with fixed ranges
    fig.update_layout(
        xaxis=dict(range=x_range, fixedrange=fixedrange) if x_range else {},
        yaxis=dict(range=y_range, fixedrange=fixedrange) if y_range else {}
    )

    return fig

import numpy as np
import plotly.graph_objects as go

def plot_coordinates(coordinate_sets, labels=None, colors=None, opacity=0.8, size=5, 
                     mode='markers'):
    """
    Plots one or multiple sets of 3D coordinates in an interactive 3D scatter plot.

    :param coordinate_sets: List of numpy arrays, where each array has shape (N, 3) representing (x, y, z) points.
    :param labels: Optional list of labels for each coordinate set.
    :param colors: Optional list of colors corresponding to each coordinate set.
    :param opacity: Opacity of the markers (default 0.8).
    :param size: Size of the markers (default 5).
    """
    fig = go.Figure()
    
    if not isinstance(coordinate_sets, list):
        coordinate_sets = [coordinate_sets]  # Convert single input into a list
    
    num_sets = len(coordinate_sets)
    
    if labels is None:
        labels = [f"Set {i+1}" for i in range(num_sets)]
    
    if colors is None:
        colors = ["blue", "red", "green", "orange", "purple"] * (num_sets // 5 + 1)  # Cycle colors
    elif type(colors) == dict: 
        colors = [colors[_key] for _key in colors]

    for i, coords in enumerate(coordinate_sets):
        coords = np.array(coords)  # Ensure it's a numpy array
        if coords.shape[1] != 3:
            raise ValueError(f"Each coordinate set must have shape (N, 3), but got {coords.shape}")

        x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]

        fig.add_trace(go.Scatter3d(
            x=x, y=y, z=z,
            mode=mode,
            name=labels[i],
            marker=dict(
                size    = size[i]    if isinstance(size, list)    else size,
                color   = colors[i]  if isinstance(colors, list)  else colors,
                opacity = opacity[i] if isinstance(opacity, list) else opacity
            )
        ))

    # Set layout
    fig = _plotly_blank_style(fig)

    return fig

    
# def visualize_voxel_grid(voxel_data, 
#                          labels, 
#                          colormap, 
#                          size=5,
#                          opacity=0.1,
#                          highlight_labels=None, 
#                          highlight_opacity=0.5):
#     """
#     Visualizes 3D voxel grids using Plotly Scatter3D.

#     :param voxel_data: 
#         List of voxel grids. Each grid is a numpy array of shape (X, Y, Z).
#     :param coordinate_labels: 
#         Dictionary where keys are labels (e.g., "BACKBONE", "BOUND_POINT") and values are voxel grids.
#     :param color_map: 
#         Dictionary mapping labels to colors for visualization.
#     :param highlight_labels: 
#         List of labels to highlight with higher opacity. Default is None.
#     """
#     fig = go.Figure()

#     for i, (label, voxel_grid) in enumerate(zip(labels, voxel_data)):
#         # Get occupied voxels
#         occupied_voxels = np.array(np.where(voxel_grid != 0)).T
        
#         if len(occupied_voxels) == 0:
#             continue  # Skip empty grids
        
#         # Extract x, y, z coordinates
#         x = occupied_voxels[:, 0]
#         y = occupied_voxels[:, 1]
#         z = occupied_voxels[:, 2]
        
#         # Add scatter plot for the current label
#         fig.add_trace(go.Scatter3d(
#             x=x, y=y, z=z,
#             mode='markers',
#             name=label,
#             marker=dict(
#                 size=size[i] if isinstance(size, list) else size,
#                 color=colormap.get(label, 'gray') if type(colormap) == dict else colormap[i],  # Default color is gray if label is not in the color_map
#                 opacity= highlight_opacity if highlight_labels and label in highlight_labels else opacity
#             )
#         ))
    
#     # Apply a blank style 
#     fig = _plotly_blank_style(fig)
    
#     return fig

def visualize_voxel_grid(voxel_data, 
                         labels=None, 
                         colormap=None, 
                         size=5,
                         opacity=0.1,
                         voxel_type='4D', 
                         highlight_labels=None, 
                         highlight_opacity=0.5,
                         property_indices=None):
    """
    Visualizes 3D voxel grids using Plotly Scatter3D. Works with both 3D and 4D voxel representations.

    Parameters:
    -----------
    voxel_data : list or array
        For voxel_type='3D': List of voxel grids, each with shape (X, Y, Z)
        For voxel_type='4D': List of voxel grids, each with shape (X, Y, Z, P) where P is properties
        For single voxel visualization, you can pass a single array
    labels : list
        List of labels corresponding to each voxel grid for legend display
        If None, will generate labels like "Voxel 1", "Voxel 2", etc.
    colormap : dict or list
        Dictionary mapping labels to colors, or a list of colors for each voxel grid
        If None, will use default Plotly colors
    size : int or list
        Size of markers. Can be a single value or a list for different sizes per voxel grid
    opacity : float
        Opacity of markers
    voxel_type : str
        '3D' - Properties encoded as integers in Z layer (shape: X, Y, Z)
        '4D' - One-hot encoded properties (shape: X, Y, Z, P)
    highlight_labels : list
        List of labels to highlight with higher opacity
    highlight_opacity : float
        Opacity for highlighted voxels
    property_indices : list or None
        For 4D voxels: Which property indices to visualize (e.g., [0, 3, 5])
        For 3D voxels: Which property values to visualize (e.g., [1, 3, 5])
        If None:
            - Will visualize all properties combined in a single trace per voxel (without property-specific labeling)
            - For 3D voxels: Will include all values except -1 (negative space)
            - For 4D voxels: Will include any non-zero property values

    Returns:
    --------
    fig : plotly.graph_objects.Figure
        The Plotly figure object for the visualization
    """
    fig = go.Figure()
    
    # Ensure voxel_data is a np.array
    if not isinstance(voxel_data, type(np.array)):
        voxel_data = np.array(voxel_data)
    
    # Generate default labels if not provided
    if labels is None:
        labels = [f"Voxel {i+1}" for i in range(len(voxel_data))]
    
    # Generate default colormap if not provided
    if colormap is None:
        # Default plotly colorscale
        colormap = cf.distinct_colors(labels)
        # colormap = {label: colors[i % len(colors)] for i, label in enumerate(labels)}
    
    # Process each voxel grid
    for i, voxel_grid in enumerate(voxel_data):
        label = labels[i] if i < len(labels) else f"Voxel {i+1}"
        
        if voxel_type == '3D':
            if property_indices is None:
                # Combined mode: plot all valid properties together
                # Get all voxels except negative space (-1)
                occupied_voxels = np.array(np.where(voxel_grid > -1)).T
                
                if len(occupied_voxels) == 0:
                    continue  # Skip if no valid voxels
                
                # Extract coordinates
                x, y, z = occupied_voxels[:, 0], occupied_voxels[:, 1], occupied_voxels[:, 2]
                
                # Determine color
                if isinstance(colormap, dict):
                    color = colormap.get(label, 'gray')
                else:
                    color = colormap[i % len(colormap)]
                
                # Determine if this should be highlighted
                is_highlighted = highlight_labels and label in highlight_labels
                
                fig.add_trace(go.Scatter3d(
                    x=x, y=y, z=z,
                    mode='markers',
                    name=label,
                    marker=dict(
                        size=size[i] if isinstance(size, list) else size,
                        color=color,
                        opacity=highlight_opacity if is_highlighted else opacity
                    )
                ))
            else:
                # Original mode: plot each property separately
                for prop_idx in property_indices:
                    # Get voxels with this property value
                    occupied_voxels = np.array(np.where(voxel_grid == prop_idx)).T
                    
                    if len(occupied_voxels) == 0:
                        continue  # Skip if no voxels have this property
                    
                    # Extract coordinates
                    x, y, z = occupied_voxels[:, 0], occupied_voxels[:, 1], occupied_voxels[:, 2]
                    
                    # Create label for this property
                    prop_label = f"{label} - Property {prop_idx}"
                    
                    # Determine color - if colormap is a dict, try to use label, then prop_label
                    if isinstance(colormap, dict):
                        color = colormap.get(prop_label, 
                               colormap.get(label, 
                               colormap.get(prop_idx, 'gray')))
                    else:
                        # If colormap is a list, use the corresponding color
                        color = colormap[i % len(colormap)]
                    
                    # Determine if this should be highlighted
                    is_highlighted = highlight_labels and (label in highlight_labels or prop_label in highlight_labels)
                    
                    fig.add_trace(go.Scatter3d(
                        x=x, y=y, z=z,
                        mode='markers',
                        name=prop_label,
                        marker=dict(
                            size=size[i] if isinstance(size, list) else size,
                            color=color,
                            opacity=highlight_opacity if is_highlighted else opacity
                        )
                    ))
                
        elif voxel_type == '4D':
            # 4D voxels - one-hot encoded properties
            if len(voxel_grid.shape) != 4:
                raise ValueError(f"Expected 4D voxel with shape (X,Y,Z,P), got {voxel_grid.shape}")
            
            if property_indices is None:
                # Combined mode: Create a single trace with all non-zero properties
                # First identify all occupied voxels (where any property is non-zero)
                # Sum across all properties
                any_property = np.sum(voxel_grid, axis=3) > 0
                occupied_voxels = np.array(np.where(any_property)).T
                
                if len(occupied_voxels) == 0:
                    continue  # Skip if no valid voxels
                
                # Extract coordinates
                x, y, z = occupied_voxels[:, 0], occupied_voxels[:, 1], occupied_voxels[:, 2]
                
                # Determine color
                if isinstance(colormap, dict):
                    color = colormap.get(label, 'gray')
                else:
                    color = colormap[i % len(colormap)]
                
                # Determine if this should be highlighted
                is_highlighted = highlight_labels and label in highlight_labels
                
                fig.add_trace(go.Scatter3d(
                    x=x, y=y, z=z,
                    mode='markers',
                    name=label,
                    marker=dict(
                        size=size[i] if isinstance(size, list) else size,
                        color=color,
                        opacity=highlight_opacity if is_highlighted else opacity
                    )
                ))
            else:
                # Original mode: plot each property separately
                for p_idx in property_indices:
                    # Extract the 3D grid for this property
                    prop_grid = voxel_grid[:, :, :, p_idx]
                    
                    # Get occupied voxels (non-zero values)
                    occupied_voxels = np.array(np.where(prop_grid > 0)).T
                    
                    if len(occupied_voxels) == 0:
                        continue  # Skip if no voxels have this property
                    
                    # Extract coordinates
                    x, y, z = occupied_voxels[:, 0], occupied_voxels[:, 1], occupied_voxels[:, 2]
                    
                    # Create label for this property
                    prop_label = f"{label} - Property {p_idx}"
                    
                    # Determine color
                    if isinstance(colormap, dict):
                        color = colormap.get(prop_label, 
                               colormap.get(label, 
                               colormap.get(p_idx, 'gray')))
                    else:
                        color = colormap[i % len(colormap)]
                    
                    # Determine if this should be highlighted
                    is_highlighted = highlight_labels and (label in highlight_labels or prop_label in highlight_labels)
                    
                    fig.add_trace(go.Scatter3d(
                        x=x, y=y, z=z,
                        mode='markers',
                        name=prop_label,
                        marker=dict(
                            size=size[i] if isinstance(size, list) else size,
                            color=color,
                            opacity=highlight_opacity if is_highlighted else opacity
                        )
                    ))
        else:
            raise ValueError(f"Unknown voxel_type: {voxel_type}. Must be '3D' or '4D'")
    
    # Apply a blank style
    fig = _plotly_blank_style(fig)
    return fig


def plot_correlation(df, x_by, y_by, 
                     label_by=None, 
                     xlabel='', ylabel='', title='', 
                     color_by=None, 
                     size=10, showlegend=True, 
                     plot_pearson_line=True, 
                     linestyle='dash', linecolor=None,
                     opacity=0.7, text_xy=[0.05, 0.95], 
                     figsize=(800, 800), 
                     padding=0.05):
    """
    Creates an interactive scatter plot with data grouped by a categorical column (e.g., `cid`).

    :param df: Pandas DataFrame containing the data.
    :param x_by: Column name for x-axis values.
    :param y_by: Column name for y-axis values.
    :param label_by: Column for hover text (optional).
    :param xlabel: String, label for the x-axis.
    :param ylabel: String, label for the y-axis.
    :param title: String, title of the plot.
    :param color_by: Column to color traces by (optional). Defaults to gray if not provided.
    :param size: Size of scatter points.
    :param showlegend: Whether to show legend entries.
    :param plot_pearson_line: Whether to show Pearson correlation line.
    :return: A Plotly figure.
    """

    # Extract relevant columns
    values1 = df[x_by].values
    values2 = df[y_by].values
    labels = df[label_by].values if label_by else [f"({x:.2f}, {y:.2f})" for x, y in zip(values1, values2)]  # Default hover text
    groups = df[color_by].unique() if color_by else ['All Data']  # Single group if color_by is missing

    # Generate distinct colors or default to gray
    if color_by:
        color_map = cf.distinct_colors(groups, category= 'tab10' if len(groups) < 10 else 'tab20')
    else:
        color_map = {'All Data': 'gray'}

    fig = go.Figure()

    # Loop over each unique group and add separate traces
    for group in groups:
        subset = df[df[color_by] == group] if color_by else df  # Use full data if no color_by
        
        hover_text = [
            f"{x_by}: {x:.3f}<br>{y_by}: {y:.3f}<br>{label}" 
            for x, y, label in zip(subset[x_by], subset[y_by], subset[label_by] if label_by else [''])
        ]
        
        fig.add_trace(go.Scatter(
            x=subset[x_by], 
            y=subset[y_by], 
            mode='markers',
            marker=dict(color=color_map[group], size=size, opacity=opacity),
            name=f"{group}" if color_by else "Data",  # Legend entry
            # text=subset[label_by] if label_by else [f"({x:.2f}, {y:.2f})" for x, y in zip(subset[x_by], subset[y_by])],
            text=hover_text,
            hoverinfo="text",
            legendgroup=f"{group}" if color_by else "Data",  # Group legend items
            showlegend=showlegend if color_by else False  # Hide legend if no grouping
        ))

    # Calculate Pearson correlation and best-fit line
    r, p_value = pearsonr(values1, values2)
    slope, intercept, _, _, _ = linregress(values1, values2)

    # Add Pearson correlation line
    if plot_pearson_line:
        sorted_idx = np.argsort(values1)
        sorted_x = values1[sorted_idx]
        sorted_y = slope * sorted_x + intercept

        fig.add_trace(go.Scatter(
            x=sorted_x, y=sorted_y, mode='lines',
            line=dict(color=linecolor if linecolor else 'black' if p_value < 0.05 else 'gray', 
                      dash=linestyle),
            text=[f"x: {x:.2f}<br>y: {y:.2f}<br>Pearson r = {r:.3f}<br>r² = {r**2:.3f}<br>p-value = {p_value:.3e}" for x, y in zip(sorted_x, sorted_y)],
            hoverinfo="text",
            name="Pearson Line",
            showlegend=showlegend
        ))

    # Add Pearson correlation text annotation
    fig.add_annotation(
        x=text_xy[0], y=text_xy[1], 
        text=f"Pearson r = {r:.3f}<br>r² = {r**2:.3f}<br>p-value = {p_value:.3e}",
        showarrow=False,
        xref="paper", yref="paper",
        # bgcolor="white", 
        # bordercolor="gray",
        opacity=0.8
    )

    # Compute axis limits with padding
    x_min, x_max = values1.min(), values1.max()
    y_min, y_max = values2.min(), values2.max()
    x_pad = padding * (x_max - x_min)
    y_pad = padding * (y_max - y_min)

    fig.update_layout(
        title=title,
        xaxis_title=xlabel or x_by,
        yaxis_title=ylabel or y_by,
        xaxis=dict(range=[x_min - x_pad, x_max + x_pad]),  # Fixed range with padding
        yaxis=dict(range=[y_min - y_pad, y_max + y_pad]),  # Fixed range with padding
        template="simple_white"
    )

    if figsize: 
        fig.update_layout(width=figsize[0], height=figsize[1])

    return fig


def plt_correlation(values1, values2,
                     xlabel='', 
                     ylabel='',
                     title='', 
                     plot_pearson_line=True, 
                     linestyle='dotted',
                     linecolor='black',
                     edgecolor='gray',
                     edgesize=10, 
                     opacity=0.5, 
                     text_xy=[0.05, 0.95], 
                     figsize=[8,8]):
    """
    Creates a scatter plot with a line of best fit and Pearson correlation annotation.

    :param values1: Array-like, the first set of values (e.g., Grantham distances).
    :param values2: Array-like, the second set of values (e.g., response correlations).
    :param xlabel: String, label for the x-axis.
    :param ylabel: String, label for the y-axis.
    :param title: String, title of the plot.
    :return: None, displays the plot.
    """
    # Ensure inputs are numpy arrays for consistency
    values1 = np.array(values1)
    values2 = np.array(values2)

    # Calculate Pearson correlation and p-value
    r, p_value = pearsonr(values1, values2)

    # Calculate the line of best fit
    slope, intercept, _, _, _ = linregress(values1, values2)

    # Create the scatter plot
    plt.figure(figsize=(figsize[0], figsize[1]))
    plt.scatter(values1, values2, color=edgecolor, alpha=opacity, s=edgesize)
    
    # Add labels, title, and grid
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title(title, fontsize=14)

    # Add text annotation for Pearson correlation, r^2, and p-value

    if plot_pearson_line: 
        sorted_idx = np.argsort(values1)  # Get indices that would sort values1
        sorted_x = values1[sorted_idx]    # Sort values1
        sorted_y = (slope * sorted_x + intercept)  # Compute y-values for sorted x

        plt.plot(sorted_x, sorted_y, color=linecolor, linestyle=linestyle, label="Pearson Line", alpha=0.6)
        plt.text(
            text_xy[0], text_xy[1],  # Adjust position as needed
            f"Pearson r = {r:.3f}\n$r^2$ = {r**2:.3f}\np-value = {p_value:.3e}",
            fontsize=12,
            ha="left", va="top",
            transform=plt.gca().transAxes,
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.5),
        )

    # Add tight layout
    plt.tight_layout()

    # Show the plot
    return plt


# def plot_weblogo(frequency_matrix, 
#                  alphabet=list("ACDEFGHIKLMNPQRSTVWY-"), 
#                  highlight_index=None, 
#                  highlight_color="gold", 
#                  highlight_text_color="gold", 
#                  highlight_text_size=10, 
#                  bw_index=None,
#                  figsize=[12, 4], 
#                  row_length=50, 
#                  plot_title='',
#                  show_axticks=True):
#     """
#     Plots a long WebLogo in multiple rows for better visualization, with optional highlighted positions.
    
#     Parameters:
#         frequency_matrix (np.ndarray): Frequency matrix (AA x Positions).
#         alphabet (list): List of amino acids including gap ('-').
#         highlight_index (dict): Mapping of highlight labels (e.g., "H104") to positions in the alignment.
#         figsize (list): Size of the figure [width, height per row].
#         row_length (int): Number of positions per row in the WebLogo.
    
#     Returns:
#         matplotlib.figure.Figure: The figure containing the WebLogo plot.
#     """
#     import matplotlib.pyplot as plt
#     import logomaker

#     # Identify and remove gaps ('-') from the frequency matrix and alphabet
#     if '-' in alphabet:
#         gap_index = alphabet.index('-')
#         frequency_matrix = np.delete(frequency_matrix, gap_index, axis=0)
#         alphabet = [aa for aa in alphabet if aa != '-']

#     # Split frequency matrix into chunks
#     num_positions = frequency_matrix.shape[1]
#     num_rows = (num_positions + row_length - 1) // row_length  # Compute the number of rows

#     # Initialize the figure
#     fig = plt.figure(figsize=(figsize[0], figsize[1] * num_rows))

#     for i in range(num_rows):
#         start = i * row_length
#         end = min((i + 1) * row_length, num_positions)

#         # Extract the chunk for the current row
#         chunk = frequency_matrix[:, start:end]
#         frequency_df = pd.DataFrame(chunk, index=alphabet)

#         # Create subplot for each row
#         ax = plt.subplot(num_rows, 1, i + 1)
#         logo = logomaker.Logo(
#             frequency_df.T,
#             ax=ax,
#             shade_below=0.5,
#             fade_below=0.5
#         )

#         # Highlight specified positions within the current chunk
#         if highlight_index:
#             for label, global_position in highlight_index.items():
#                 global_position -= 1
#                 if start <= global_position < end:
#                     local_position = global_position - start
#                     logo.highlight_position(p=local_position, color=highlight_color, alpha=0.5)
#                     ax.text(local_position, -0.1, label, color=highlight_text_color, 
#                             ha="center", fontsize=highlight_text_size, rotation=0)

#         if bw_index:
#             for label, global_position in bw_index.items():
#                 global_position -= 1
#                 if start <= global_position < end:
#                     local_position = global_position - start
#                     logo.highlight_position(p=local_position, color='gray', alpha=0.5)
#                     ax.text(local_position, -0.05, label, color="gray", ha="center", fontsize=10, rotation=0)

#         # Style and label the subplot
#         logo.style_spines(visible=False)
#         logo.style_spines(spines=["left", "bottom"], visible=True)
#         ax.set_ylabel("Frequency")
#         # ax.set_xlabel(f"Position {start + 1}-{end}")
#         # Set ticks at the first and last position of the chunk
#         if show_axticks: 
#             ax.set_xticks([0, chunk.shape[1] - 1])  # Positions relative to the chunk (0-based index)
#             ax.set_xticklabels([start + 1, end])    # Labels corresponding to the actual sequence positions
#         else: # Add blank to x ticks 
#             ax.set_xticks([0])
#             ax.set_xticklabels([''])
#         ax.set_title(plot_title) if i == 0 else None # Only print label in the first plot 

#     plt.tight_layout()

#     return fig