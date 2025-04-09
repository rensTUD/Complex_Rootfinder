# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 15:20:28 2025

@author: RensvanLeijden
"""

# %% IMPORTS

# imports of python libraries
import numpy as np
from typing import Callable
from scipy.fft import fft, ifft
from scipy.linalg import hankel, eig

# imports from curent library
from .utils import finite_difference_first_derivative_4th_order
from .contours import ContourBase

# %% FUNCTION DEFNITIONS

def newtons_identities(s_N, n_roots):
    """
    Compute the coefficients of a monic polynomial from power sums using Newton's identities.

    Parameters
    ----------
    s_N : array_like
        Power sums (moments) s_1, s_2, ..., s_N.
    n_roots : int
        Degree of the polynomial (number of roots).

    Returns
    -------
    p : ndarray
        Polynomial coefficients in descending order (highest degree first).
    """
    e_N = np.ones(n_roots + 1, dtype=complex)
    for iN in range(1, n_roots + 1):
        e = 0
        for ii in range(1, iN + 1):
            e += (-(-1)**ii) * e_N[iN - ii] * s_N[ii - 1]
        e_N[iN] = e / iN

    # Construct polynomial coefficients
    p = np.zeros(n_roots + 1, dtype=complex)
    for i in range(n_roots + 1):
        p[i] = (-1)**i * e_N[i]
    
    return p


def find_roots_delves_lynes(
    f: Callable[[np.ndarray], np.ndarray],
    contour: ContourBase,
    n_roots: int,
    df: Callable[[np.ndarray], np.ndarray] | None = None,
    previous_roots: np.ndarray | None = None
):
    """
    Finds the roots based on the polynomial created by newtons identities using the argument principle.
    Local deflation is applied when previously found roots are located within the domain

    Parameters
    ----------
    f : callable 
            the function to analyze
    contour: ContourBase
            Contour class       
    n_roots : int
        Number of roots to be found, determines the order of the polynomial
    df : Callable[[np.ndarray], np.ndarray] | None, optional
        the derivative of f, by default None
    previous_roots : np.array | None, optional
        numpy array of previously found roots in the domain within Z, by default None

    Returns
    -------
    roots
        list of the found roots
    """
    # get Z
    Z = contour.Z
    
    # Get the derivative
    if df is None:
        df = finite_difference_first_derivative_4th_order(f, Z)
    else:
        df = df(Z)
    
    # Calculate s_N
    s_N = np.zeros(n_roots, dtype=complex)

    for iN in range(1, n_roots + 1):
            s_N[iN - 1] = np.trapz((Z ** iN) * df / f(Z), Z) / (2 * np.pi * 1j)
            
    # If previous roots are found then apply local deflation
    if previous_roots:
        n_roots = n_roots - len(previous_roots)
        for iN in range(1,n_roots+1):
            s_N[iN] = s_N[iN] - np.sum(previous_roots**iN)
    
    # calculate the roots based on the roots of the polynomial
    p = newtons_identities(s_N, n_roots)
    
    roots = np.roots(p)
    
    return roots

def find_roots_austin_kravanja(
    f: Callable[[np.ndarray], np.ndarray],
    contour: ContourBase,
    n_roots: int,
    step_size: float = 1e-6,
    max_points: int = 2**14,
    df: Callable[[np.ndarray], np.ndarray] | None = None,
    previous_roots: np.ndarray | None = None
):
    """
    Finds the roots based on the cauchy integral on a unit disk

    Parameters
    ----------
    f : callable 
            the function to analyze
    contour: ContourBase
            Contour class  
    n_roots : int
        Number of roots to be found, determines the order of the polynomial
    df : Callable[[np.ndarray], np.ndarray] | None, optional
        the derivative of f, by default None
    previous_roots : np.array | None, optional
        numpy array of previously found roots in the domain within Z, by default None

    Returns
    -------
    roots
        list of the found roots
    """
    # get radius and center from contour
    radius = contour.radius
    center = contour.center
    
    # Estimate number of points based on arc length / step_size
    arc_length = 2 * np.pi * radius
    n_est = int(np.ceil(arc_length / step_size))

    # Use next power of 2, capped by max_points
    Npoints = 2 ** int(np.floor(np.log2(min(n_est, max_points))))

    # Discrete roots of unity (equally spaced points on circle)
    k = np.arange(Npoints)
    theta = 2 * np.pi * k / Npoints
    Z = center + radius * np.exp(1j * theta)

    fk = f(Z)
    c = fft(fk) / Npoints

    # Derivative of FFT terms (equivalent to analytic derivative)
    cp = np.arange(1, Npoints) * c[1:]
    ppzk = Npoints * ifft(np.concatenate((cp, [0]))) / radius
    
    s = ifft(ppzk/fk)
    
    # Hankel matrices
    col1_H = s[1:n_roots+1]            
    row1_H = s[n_roots:2*n_roots]            
    H = hankel(col1_H, row1_H)

    col1_H2 = s[2:n_roots+2]      
    row1_H2 = s[n_roots+1:2*n_roots+1]  
    H2 = hankel(col1_H2, row1_H2)
    
    # get eigenvalues as roots and scale with r0 and translate by z0
    eigvals = eig(H2, H, right=False)
    roots = radius * eigvals + center
    
    return roots

