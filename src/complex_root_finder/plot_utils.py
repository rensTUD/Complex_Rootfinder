# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 11:03:14 2025

@author: rensv
"""

# %% IMPORTS
import numpy as np
import matplotlib.pyplot as plt
    
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
