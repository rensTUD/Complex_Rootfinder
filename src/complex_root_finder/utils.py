# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 17:27:51 2025

@author: RensvanLeijden
"""

# %% IMPORTS

import numpy as np
from typing import Callable


# %% FUNCTIONS

def finite_difference_first_derivative_4th_order(
        f: Callable[[np.ndarray], np.ndarray], 
        Z
    ):
        """
        Compute derivative using 4th order accurate finite difference.
        
        Parameters
        ----------
            f : callable 
                the function to differentiate
            Z : array-like
                points at which to evaluate the derivative
        
        Returns:
        -------
            array-like
                The approximated derivative values
        """
        h = 1e-5
        return (-f(Z + 2*h) + 8*f(Z + h) - 8*f(Z - h) + f(Z - 2*h)) / (12*h)