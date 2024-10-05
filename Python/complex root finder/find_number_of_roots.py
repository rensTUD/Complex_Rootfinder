# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 20:47:51 2024

@author: RensvanLeijden
"""

# %%

import numpy as np
from scipy.integrate import simps
from scipy.fft import fft, ifft

# %%

def find_number_of_roots(f, Z, option_der='numerical', df=None, **kwargs):
    """
    Estimates the number of roots of a complex-valued function `f` in a given domain.
    
    Parameters:
    f : callable
        The complex-valued function for which roots need to be found.
    Z : numpy.ndarray
        A path along which the function is evaluated.
    option_der : str, optional
        Method to use for calculating the derivative. Options are 'symbolic', 'numerical', 'argp', or 'rootsofunity'.
    df : callable, optional
        Derivative function of `f`, used if `option_der` is 'symbolic'.
    kwargs : dict, optional
        Additional arguments required for specific options like 'argp' or 'rootsofunity'.
    
    Returns:
    Nroots : int
        The estimated number of roots.
    """

    # Calculate the number of roots based on the chosen derivative method
    if option_der == 'symbolic' and df is not None:
        # Use the provided symbolic derivative
        Nroots = simps(df(Z) / f(Z), Z) / (2 * np.pi * 1j)
    
    elif option_der == 'numerical':
        # Numerically calculate the derivative
        dfdz = np.diff(f(Z)) / np.diff(Z)
        Nroots = simps(dfdz / f(Z[:-1]), Z[:-1]) / (2 * np.pi * 1j)
    
    elif option_der == 'argp':
        # Derivation without explicit derivative (Peter & Kravanja method)
        Npoints = kwargs.get('Npoints', len(Z))
        Args = np.zeros(Npoints, dtype=np.complex_)
        for ii in range(Npoints - 1):
            Args[ii] = np.angle(f(Z[ii + 1]) / f(Z[ii]))
        Args[-1] = np.angle(f(Z[0]) / f(Z[-1]))
        Nroots = np.sum(Args) / (2 * np.pi)
    
    elif option_der == 'rootsofunity':
        z0 = kwargs['z0']
        r0 = kwargs['r0']
        step_size = 1e-6
        Npoints = int(np.round(2 * np.pi * r0 / step_size))
        Npoints = min(Npoints, 10000)  # Keep within max limit
        k = np.arange(Npoints)
        Z = z0 + r0 * np.exp(2 * np.pi * 1j * k / Npoints)
        
        n = Npoints
        fk = f(Z)
        c = fft(fk) / n
        cp = np.arange(1, n) * c[1:]
        ppzk = n * ifft(np.concatenate(([0], cp))) / r0

        # By Austin / Kravanja method
        Nroots = np.real(np.mean(Z * ppzk / fk))
    
    else:
        raise ValueError(f"Unknown option_der: {option_der}")

    # Handle branch cuts and adjust Nroots accordingly
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
