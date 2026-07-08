import matplotlib.pyplot as plt

import statistics as stats
from typing import List, Tuple
import math

# calculate the radius of a regular polygon
polygon_radius = lambda num_sides, side_length: side_length / (2 * math.sin(math.pi / num_sides))

def get_info_from_contig(contig:str) -> Tuple[List[int], List[int]]:
    """Split up contig to extract domain and linkage information."""
    contig_split = contig.split('-')
    
    # either even iteration describes a polygon (or domain)
    polygons = [int(n) for n in contig_split[::2]]
    
    # and every odd iteration describes its connectivity
    shared_helices = [int(n) for n in contig_split[1::2]]
    
    return (polygons, shared_helices)

def combine_domain_angles(angles:List[List[int]], shared_helices:List[int],
                          windings:List[str]=None) -> List[float]:
    """Combine per-domain polygon angles into one turtle path.

    windings: optional list of 'CW'/'ACW' per domain. The winding is seeded from the
    FIRST domain's handedness ('ACW' = anticlockwise/positive, 'CW' = clockwise/negative)
    and then alternates between stacked domains so they don't overlap. For chains whose
    handedness alternates (all even-size modules) this reproduces each domain's true
    handedness. If None, the winding just alternates starting anticlockwise (original).
    """
    combined_angles = []

    # need to add a dummy value here to simplify the code
    shared_helices = [0] + shared_helices
    seed = 1 if (windings is None or windings[0] == 'ACW') else -1

    for i, (angle_set, shared) in enumerate(zip(angles, shared_helices)):

        sign = seed * (1 if i % 2 == 0 else -1)
        corrected_angles = [sign * angle for angle in angle_set]

        # take one angle away for the first domain (otherwise we'll make a circle)
        if i == 0:
            combined_angles = combined_angles + corrected_angles

        elif shared == 1:
            combined_angles = combined_angles + [sign * 26.7] + corrected_angles[2:]

        elif shared == 2:
            combined_angles = combined_angles + corrected_angles[2:]

    return combined_angles

def calc_cart_coords(combined_angles:List[int]) -> Tuple[List[int], List[int]]:
    """Using trigonometry, convert polar coordinates to cartesian coordinates."""
    current_x = 0
    current_y = 0
    current_angle = 0
    
    # This is the side length. Keep at 11 for consistency between figures
    side_lengths = [11] * len(combined_angles)

    # initialise empty lists for x and y cart coords
    x = []
    y = []
    
    for side_length, angle in zip(side_lengths, combined_angles):
        
        # do the maths
        current_angle += math.radians(angle)
        current_x += side_length * math.cos(current_angle)
        current_y += side_length * math.sin(current_angle)
        x.append(current_x)
        y.append(current_y)
    
    return x, y

def create_plot(x:List[int], y:List[int], add_labels:bool=False, order:List[int]=None, outpath:str=None, ax=None, N_colour:str = "#0077bb", C_colour:str = "#cc3311",loops=True) -> None:
    """Draws a plot from cartesian coordinates representing helices, loops and colours N/C-terms.
    
    Args:
        x (List[int]): x-coordinate for helix.
        y (List[int]): y-coordinate for helix.
        add_labels (bool): to add labels to circles for relooping, default = False
        order (List[int]): list of helix numbers for looping, default = None
        outpath (str): where to save the file, if not provided, doesn't save, default = None
        ax: plot object to draw on, useful for iterations in multiplots, default = None
        N_colour (str): colour of the N-terminal helix
        C_colour (str): colour of the C-terminal helix
        
    Returns:
        None: but will save a graphic using the outpath argument
    """
    
    def add_circle(plot_obj, x, y, size, fc, ec, lw, label=None):
        """Adds circle at x, y coordinates with a few customisable parameters"""
        circle = plt.Circle((x, y), size, color=fc, fill=True, ec=ec, lw=lw)
        plot_obj.add_artist(circle)
        if label is not None:
            # Adding the label to the center of the circle
            plot_obj.text(x, y, str(label), ha='center', va='center', size=20, weight='bold', color='w')

    # reorder x, y coordinates for correct connections
    if order is None:
        order = range(len(x))
    x_ordered = [x[ord] for ord in order]
    y_ordered = [y[ord] for ord in order]   
                
    # if not plot_obj is provide, makes one
    if ax is None:
        ax = plt.subplot()

    # draws the "loops" (sides) connecting the "helices" (vertices)
    if loops:
        for i, _ in enumerate(x_ordered):
            if i % 2==0:
                ax.plot(x_ordered[i:i+2], y_ordered[i:i+2], lw=6, c='k', zorder=2)
            else:
                ax.plot(x_ordered[i:i+2], y_ordered[i:i+2], lw=6, c='k', zorder=0)

    # Add a circle at each vertex
    for i, (x_coord, y_coord) in enumerate(zip(x_ordered, y_ordered)):
        
        # specifies label number if needed
        if add_labels:
            label = i
        else:
            label = None
        
        # each cirlce has different colouring combinations depending which
        # terminus we can see from the forced perspective. By default the N-term
        # of the sequence starts away from us.
        
        # if N-term of sequence, add N_colour outer ring to denote
        if loops and i == 0:
            add_circle(ax, x_coord, y_coord, 4.5, N_colour, None, 1, label=label)
            add_circle(ax, x_coord, y_coord, 3, C_colour, 'k', 4)

        # if helix iteration is even i.e, we see C-term of helix
        # from our perspective, colour circle C_colour
        elif i % 2 == 0:
            
            # if C-term of sequence, add C_colour outer ring to denote
            if loops and i == len(x)-1:
                add_circle(ax, x_coord, y_coord, 4.5, C_colour, None, 1, label=label)
                add_circle(ax, x_coord, y_coord, 3, C_colour, 'k', 4)

            # else make outer ring white
            else:
                    add_circle(ax, x_coord, y_coord, 4.5, "w", None, 1, label=label)
                    add_circle(ax, x_coord, y_coord, 3, C_colour, 'k', 4)

        # else we can see N-term of helix from our perspective, colour circle N_colour
        else:
            
            # if C-term of sequence, add C_colour outer ring to denote
            if loops and i == len(x)-1:
                add_circle(ax, x_coord, y_coord, 4.5, C_colour, None, 1, label=label)
                add_circle(ax, x_coord, y_coord, 3, N_colour, 'k', 4)

            # else make outer ring white
            else:
                add_circle(ax, x_coord, y_coord, 4.5, "w", None, 1, label=label)
                add_circle(ax, x_coord, y_coord, 3, N_colour, 'k', 4)
    
    # now we try to make all plots consistent
    # first we get the geometric centre of our coordinates
    
    x_bar = stats.mean(x)
    y_bar = stats.mean(y)

    # then using an arbitrary scale, set x,y limits consistent for all plots
    # play with this is if you're making some crazy big toplogies
    arbitrary_radius = 37.5
    
    # set the axis limits from this value
    ax.set_xlim(x_bar-arbitrary_radius, x_bar+arbitrary_radius)
    ax.set_ylim(y_bar-arbitrary_radius, y_bar+arbitrary_radius)
    
    # make even aspect
    ax.set_aspect('equal', adjustable='box')
    
    # remove axes
    ax.axis('off')
    
    # save file
    plt.tight_layout()
    if outpath != None:
        fig_file_format = outpath.split('.')[-1]
        plt.savefig(outpath, transparent=True, format=fig_file_format, dpi=400)
    plt.show()
