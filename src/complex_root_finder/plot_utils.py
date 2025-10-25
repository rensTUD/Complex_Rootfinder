# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 11:03:14 2025

@author: rensv
"""

# %% IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString, MultiLineString
    
# %% FUNCTIONS

def debug_plot_contours(
    contours: list,
    color: str = 'r',
    alpha: float = 0.5,
    label: str = None,
    pause: float = 0.01,
    ax=None,
):
    """
    Plot a list of ContourBase instances for debugging.

    Parameters
    ----------
    contours : list
        List of contours (RectangleContour or CircleContour).
    color : str, optional
        Color of the contour lines.
    alpha : float, optional
        Transparency.
    label : str, optional
        Optional label for the contour (added only once).
    pause : float, optional
        Pause duration for interactive plotting.
    ax : matplotlib axis, optional
        Axis to plot on. If None, uses current axis.
    """


    ax = ax or plt.gca()

    for i, c in enumerate(contours):
        Z = c.Z
        ax.plot(np.real(Z), np.imag(Z), color=color, alpha=alpha, label=label if i == 0 else None)

    if label:
        ax.legend()

    plt.pause(pause)

def plot_boundary_orientation(contour, num_arrows=10):
    """
    Plot the .boundary of a contour with arrows showing orientation.

    Parameters
    ----------
    contour : object
        Must have a `.boundary` attribute (shapely Polygon).
    num_arrows : int
        Number of arrows to show along the boundary.
    """
    boundary = contour.boundary
    
    if not isinstance(boundary, Polygon):
        raise TypeError(f"Expected boundary to be a Polygon, got {type(boundary)}")
    
    coords = list(boundary.exterior.coords)
    x, y = zip(*coords)

    fig, ax = plt.subplots()
    ax.plot(x, y, '.', label='Boundary')

    # Plot arrows
    total_points = len(x) - 1  # last point repeats
    step = max(1, total_points // num_arrows)

    for i in range(0, total_points, step):
        ax.annotate(
            '', 
            xy=(x[i+1], y[i+1]), 
            xytext=(x[i], y[i]),
            arrowprops=dict(arrowstyle='->', color='red', lw=1),
            size=15
        )

    # ax.set_aspect('equal')
    ax.legend()
    plt.title(f'Boundary Orientation ({type(contour).__name__})')
    plt.show()