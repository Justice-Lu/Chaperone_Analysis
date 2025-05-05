import numpy as np 
import plotly.colors as pc 
import random
from matplotlib import cm, colors as mcolors

from PIL import ImageColor


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
        


def get_continuous_colors(values, colormap="RdBu_r", midpoint=None):
    """
    Assigns colors to a list or dictionary of values using a continuous colorscale.

    Parameters:
    - values (list, dict, or Series): A dictionary {index: value}, Pandas Series, or a list of values.
    - colormap (str): A Matplotlib colormap name (e.g., 'RdBu_r', 'viridis', 'coolwarm').
    - midpoint (float, optional): The value to center the colormap around. If None, uses the median.

    Returns:
    - dict: A dictionary mapping indices (or positions if input is a list) to colors in hex format.
    """
    # Convert Pandas Series to dictionary if applicable
    if hasattr(values, "to_dict"):
        values = values.to_dict()

    # Extract indices and values
    if isinstance(values, dict):
        indices, vals = list(values.keys()), np.array(list(values.values()))
    else:
        indices, vals = range(len(values)), np.array(values)

    # Handle edge case where all values are the same
    if np.all(vals == vals[0]):
        return {idx: mcolors.to_hex(cm.get_cmap(colormap)(0.5)) for idx in indices}

    # Define midpoint if not provided (default: median)
    if midpoint is None:
        midpoint = np.median(vals)

    # Determine min and max values
    min_val, max_val = np.min(vals), np.max(vals)

    # Adjust the color scale if all values are on one side of the midpoint
    if min_val >= midpoint:  # All values are positive
        max_val = max(max_val, abs(midpoint))  # Ensure symmetry
        min_val = -max_val
    elif max_val <= midpoint:  # All values are negative
        min_val = min(min_val, -abs(midpoint))  # Ensure symmetry
        max_val = -min_val

    # Normalize values between [0, 1] while keeping midpoint centered
    if min_val < midpoint < max_val:
        norm = mcolors.TwoSlopeNorm(vmin=min_val, vcenter=midpoint, vmax=max_val)
    else:
        norm = mcolors.Normalize(vmin=min_val, vmax=max_val)

    # Get the colormap and apply normalization
    cmap = cm.get_cmap(colormap)
    colors = {idx: mcolors.to_hex(cmap(norm(value))) for idx, value in zip(indices, vals)}

    return colors

def distinct_colors(label_list=None, category='tab10', custom_color=None, random_state=0, form='dict', num_colors=None):
    """
    Generate distinct colors for a list of labels or a specified number of colors.

    Parameters:
    label_list (list, optional): A list of labels for which you want to generate distinct colors.
                                If None, will generate num_colors distinct colors.
    category (str): Category of distinct colors. Options are 'warm', 'floral', 'rainbow', 'pastel',
                    matplotlib color palettes (e.g., 'tab10', 'Set2'), or 'random'. Default is 'tab10'.
    custom_color (list, optional): A custom list of colors to use.
    random_state (int): Seed for random color generation. Default is 0.
    form (str): Output format - 'dict' for a label-to-color dictionary, 'list' for a list of colors. Default is 'dict'.
    num_colors (int, optional): Number of colors to generate if label_list is None.

    Returns:
    dict or list: Colors in the requested format (dictionary mapping labels to colors or list of colors).
    """
    
    random.seed(random_state)
    np.random.seed(random_state)
    
    # Predefined color palettes
    warm_colors = ['#fabebe', '#ffd8b1', '#fffac8', '#ffe119', '#ff7f00', '#e6194B']
    floral_colors = ['#bfef45', '#fabed4', '#aaffc3', '#ffd8b1', '#dcbeff', '#a9a9a9']
    rainbow_colors = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4']
    pastel_colors = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', 
                     '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928', 
                     '#8DD3C7', '#BEBADA', '#FFED6F']
    
    # Determine the number of colors needed
    if label_list is not None:
        required_colors = len(label_list)
    elif num_colors is not None:
        required_colors = num_colors
    else:
        raise ValueError("Either label_list or num_colors must be provided")
    
    # Handle custom colors
    if custom_color is not None:
        if len(custom_color) < required_colors:
            raise ValueError(f"Not enough custom colors ({len(custom_color)}) for the required number ({required_colors})")
        color_palette = custom_color[:required_colors]
    
    # Generate color palette based on category
    else:
        if category == 'warm':
            base_palette = warm_colors
        elif category == 'floral':
            base_palette = floral_colors
        elif category == 'rainbow':
            base_palette = rainbow_colors
        elif category == 'pastel':
            base_palette = pastel_colors
        elif category == 'random':
            # Generate completely random colors
            color_palette = ["#{:06x}".format(random.randint(0, 0xFFFFFF)) for _ in range(required_colors)]
        else:
            # Try to use matplotlib colormap
            try:
                cmap = cm.get_cmap(category)
                base_palette = [mcolors.to_hex(cmap(i)) for i in np.linspace(0, 1, min(required_colors, 20))]
            except (ValueError, AttributeError):
                # Fallback to tab10 if the specified category is not available
                cmap = cm.get_cmap('tab10')
                base_palette = [mcolors.to_hex(cmap(i)) for i in np.linspace(0, 1, min(required_colors, 10))]
        
        # If category wasn't 'random', ensure we have enough unique colors
        if category != 'random':
            color_palette = []
            
            # If we need more colors than in the base palette, use interpolation and HSV manipulation
            if required_colors > len(base_palette):
                # Start with all colors from the base palette
                color_palette = base_palette.copy()
                
                # Convert to HSV for better interpolation and manipulation
                hsv_colors = [mcolors.rgb_to_hsv(mcolors.to_rgb(color)) for color in base_palette]
                
                # Generate additional colors by manipulating hue and saturation
                while len(color_palette) < required_colors:
                    new_hsv = hsv_colors[len(color_palette) % len(hsv_colors)].copy()
                    # Modify hue and saturation slightly
                    new_hsv[0] = (new_hsv[0] + 0.1 * (len(color_palette) // len(hsv_colors))) % 1.0
                    new_hsv[1] = max(0.4, min(1.0, new_hsv[1] + 0.05 * ((len(color_palette) // len(hsv_colors)) % 3 - 1)))
                    
                    # Convert back to RGB, then hex
                    new_rgb = mcolors.hsv_to_rgb(new_hsv)
                    new_hex = mcolors.to_hex(new_rgb)
                    
                    # Only add if the color is visually distinct enough (simple check)
                    if all(mcolors.rgb_to_hsv(mcolors.to_rgb(new_hex))[0] != 
                           mcolors.rgb_to_hsv(mcolors.to_rgb(existing))[0] for existing in color_palette[-10:]):
                        color_palette.append(new_hex)
                    else:
                        # If too similar, add some randomness to the hue
                        new_hsv[0] = (new_hsv[0] + random.random() * 0.2) % 1.0
                        new_rgb = mcolors.hsv_to_rgb(new_hsv)
                        color_palette.append(mcolors.to_hex(new_rgb))
            else:
                # If we have enough colors in the base palette, just use those
                color_palette = base_palette[:required_colors]
    
    # Return results in the requested format
    if label_list is None:
        return color_palette if form == 'list' else {i: color for i, color in enumerate(color_palette)}
    else:
        if form == 'list':
            return color_palette
        else:
            return {label: color_palette[i] for i, label in enumerate(label_list)}

from plotly.validators.scatter.marker import SymbolValidator
import random

from plotly.validators.scatter.marker import SymbolValidator
import random

def distinct_shapes(label_list=None, random_state=0, form='dict', num_shapes=None):
    """
    Generate distinct shapes for a list of labels or a specified number of shapes.

    Parameters:
    label_list (list, optional): A list of labels for which you want to generate distinct shapes.
                                 If None, will generate num_shapes distinct shapes.
    random_state (int): Seed for reproducibility. Default is 0.
    form (str): Output format - 'dict' for a label-to-shape dictionary, 'list' for a list of shapes. Default is 'dict'.
    num_shapes (int, optional): Number of shapes to generate if label_list is None.

    Returns:
    dict or list: Shapes in the requested format (dictionary mapping labels to shapes or a list of shapes).
    """
    
    random.seed(random_state)

    # Get all available marker symbols in Plotly
    all_shapes = [s for i, s in enumerate(SymbolValidator().values) if i % 3 == 2]  

    # Separate main shapes (without '-') and sub-shapes (with '-')
    main_shapes = [s for s in all_shapes if '-' not in s]
    sub_shapes = [s for s in all_shapes if '-' in s]

    # Determine number of required shapes
    if label_list is not None:
        required_shapes = len(label_list)
    elif num_shapes is not None:
        required_shapes = num_shapes
    else:
        raise ValueError("Either label_list or num_shapes must be provided.")

    # Prioritize main shapes, then use sub-shapes if needed
    shape_list = main_shapes[:required_shapes]  # Take as many main shapes as available
    if len(shape_list) < required_shapes:
        remaining = required_shapes - len(shape_list)
        shape_list.extend(sub_shapes[:remaining])  # Fill remaining slots with sub-shapes

    # Ensure cycling if still not enough
    if len(shape_list) < required_shapes:
        combined_shapes = main_shapes + sub_shapes
        shape_list = [combined_shapes[i % len(combined_shapes)] for i in range(required_shapes)]

    # Return the requested format
    if label_list is None:
        return shape_list if form == 'list' else {i: shape for i, shape in enumerate(shape_list)}
    else:
        return {label: shape_list[i] for i, label in enumerate(label_list)} if form == 'dict' else shape_list

def scale(values, reverse=False, factor = 1, scale_between = [1,0]):
    """
    Reverses the scale of a list of values such that the smallest value becomes 1 and the largest value becomes 0.
    """
    
    min_val = min(values)
    max_val = max(values)
    new_max = scale_between[0]
    new_min = scale_between[1]

    if min_val != max_val: #if all the min and max is the same value. assign  max size 
        scaled_values = [(value - min_val) * (new_max - new_min) / (max_val - min_val) + new_min for value in values]
    else: 
        scaled_values = [value*new_max for value in values]
#     scaled_values = [(val - min_val) / (max_val - min_val)*factor for val in values]
    if reverse:
        scaled_values = [1 - val for val in scaled_values]
    return scaled_values
