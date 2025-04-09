# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 20:47:51 2024

@author: RensvanLeijden
"""

# %%

import numpy as np
from scipy.integrate import trapezoid
from scipy.fft import fft, ifft
from typing import Callable

# %%

def approx_deriv(
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

def count_roots_numerical(
    f: Callable[[np.ndarray], np.ndarray], 
    Z, 
    df = None
):
    """
    Find number of roots using numerical derivative.
    
    Parameters
    ----------
        f : callable 
            the function to analyze
        Z : array-like
            points along the path
        df : callable, optional
            the derivative of f
    
    Returns:
    -------
        float 
            The number of roots
    """
    if df is None:
        df = approx_deriv(f, Z)
    else:
        df = df(Z)
    
    y = df / f(Z)

    n_roots = np.trapezoid(y, x=Z) / (2 * np.pi * 1j)
        
    return process_root_count_result(n_roots)

def count_roots_unity(
    f: Callable[[np.ndarray], np.ndarray],
    z0: complex,
    r0: float,
    step_size: float = 1e-6,
    max_points: int = 2**14,
) -> float:
    """
    Estimate the number of roots of an analytic function inside a circular contour
    using the roots of unity method and FFT.

    Parameters
    ----------
    f : Callable[[np.ndarray], np.ndarray]
        Analytic function to analyze, vectorized over complex inputs.
    z0 : complex
        Center of the circular contour.
    r0 : float
        Radius of the contour.
    step_size : float, optional
        Step size for arc length sampling (default is 1e-6).
    max_points : int, optional
        Maximum number of FFT points (default is 2**14).

    Returns
    -------
    float
        Estimated number of roots inside the contour.
    """
    # Estimate number of points based on arc length / step_size
    arc_length = 2 * np.pi * r0
    n_est = int(np.ceil(arc_length / step_size))

    # Use next power of 2, capped by max_points
    Npoints = 2 ** int(np.floor(np.log2(min(n_est, max_points))))

    # Discrete roots of unity (equally spaced points on circle)
    k = np.arange(Npoints)
    theta = 2 * np.pi * k / Npoints
    Z = z0 + r0 * np.exp(1j * theta)

    fk = f(Z)
    c = fft(fk) / Npoints

    # Derivative of FFT terms (equivalent to analytic derivative)
    cp = np.arange(1, Npoints) * c[1:]
    ppzk = Npoints * ifft(np.concatenate((cp, [0]))) / r0

    # Apply formula: mean of Z * f'(Z) / f(Z) to get number of roots
    n_roots =  np.real(np.mean(Z * ppzk / fk))
    
    return process_root_count_result(n_roots)

def process_root_count_result(Nroots):
    """
    Process and validate the root count result.
    
    Parameters
    ----------
        Nroots : complex or float
            raw root count result
    
    Returns:
    -------
        int or complex
            Processed number of roots
    """
    if np.isnan(Nroots) or np.isinf(Nroots):
        return 0

    # Check for imaginary value, indicating crossing a branch cut
    if round(np.imag(Nroots) / (1 / np.pi)) != 0:
        N_BC = round(np.imag(Nroots) / (1 / np.pi))
        
        if round(np.real(Nroots)) == 0:
            return 0
        
        if np.mod(np.real(Nroots), 1) != 0:
            # When values are 1.5, 2.5, etc.
            Nroots = round(np.real(Nroots), 2) + 1j * N_BC
        else:
            # When value is 0.5 (indicating one branch point)
            Nroots = round(np.real(Nroots), 1) + 1j * N_BC
        
        return Nroots

    # Round Nroots to the nearest integer
    return int(round(np.real(Nroots)))


