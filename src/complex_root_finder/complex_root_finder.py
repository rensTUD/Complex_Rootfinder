# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 15:51:28 2025

@author: RensvanLeijden
"""
# %% IMPORTS

# imports from python libraries
from typing import Callable, List, Tuple
import numpy as np
import matplotlib.pyplot as plt
from collections import deque
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, squareform

# imports from current library
from .utils import finite_difference_first_derivative_4th_order
from .contours import ContourBase, RectangleContour, CircleContour
from .count_roots import count_roots_numerical, count_roots_unity
from .find_roots import find_roots_delves_lynes, find_roots_austin_kravanja

# %% CLASS DEFINITION


class ComplexRootFinder:
    def __init__(
        self,
        f: Callable[[np.ndarray], np.ndarray],
        df: Callable[[np.ndarray], np.ndarray] | None = None,
    ):
        """
        Initialize a ComplexRootFinder instance.

        Parameters
        ----------
        f : Callable
            Function to analyze. Can be scalar-valued or matrix-valued.
        df : Callable, optional
            Derivative of f. Can be scalar- or matrix-valued. If None, will be approximated.
        """
        # assign values to self
        self.original_f = f
        self.original_df = df

        # process input
        self._process_function_inputs()

    def _process_function_inputs(
        self
    ):
        """
        Determines whether f and df are scalar- or matrix-valued and processes them accordingly.
        Sets self.f and self.df as callables returning scalar values.
        """
        # random test value
        test_z = 1+1j

        # Handle f
        f_output = self.original_f(test_z)
        is_matrix_f = isinstance(f_output, np.ndarray) and f_output.ndim == 2

        if is_matrix_f:
            self.f = np.vectorize(lambda z: np.linalg.det(self.original_f(z)))
        else:
            self.f = self.original_f

        # Handle df 
        if self.original_df is None:
            self.df = lambda z: finite_difference_first_derivative_4th_order(self.f, z)  # Use finite difference if no derivative provided
        else:
            df_output = self.original_df(test_z)

            if isinstance(df_output, np.ndarray):
                if df_output.ndim == 2:  # Matrix-valued -> compute derivative of determinant 
                    self.df = self._wrap_det_derivative(self.original_f, self.original_df)
                elif df_output.ndim == 0 or df_output.shape == ():  # Scalar array
                    self.df = self.original_df
                else:
                    raise ValueError("df appears to be a vector or tensor — expected scalar or matrix.")
            elif np.isscalar(df_output):  # float, complex, etc.
                self.df = self.original_df
            else:
                raise TypeError(f"Unsupported df return type: {type(df_output)}")
                
    def _wrap_det_derivative(
        self, 
        f_matrix, 
        df_matrix
    ):
        def df_scalar(z):
            z = np.atleast_1d(z)
            A = f_matrix(z)        # shape (N, m, m)
            A_prime = df_matrix(z) # shape (N, m, m)

            # Compute inverse of each matrix in A: shape (N, m, m)
            A_inv = np.linalg.inv(A)

            # Compute trace of A_inv @ A_prime for each i
            traces = np.einsum('nij,njk->nik', A_inv, A_prime)  # matrix product
            traces = np.trace(traces, axis1=1, axis2=2)         # trace of each matrix

            # Compute determinant of each A
            dets = np.linalg.det(A)

            return traces * dets  # final result is shape (N,)
        
        return df_scalar
    
    def count_roots_domains_rectangle(
        self,
        real_min: float,
        real_max: float,
        imag_min: float,
        imag_max: float
    ):
        
        contour = RectangleContour(
            real_min, 
            real_max, 
            imag_min,
            imag_max
            )
        
        domains_number_of_roots = count_root_containing_domains(
            self.f, 
            self.df, 
            contour,
            count_roots_fn = count_roots_numerical,
            debug=True)
        
        return domains_number_of_roots
    
    def find_roots_domains_rectangle(
        self,
        real_min: float,
        real_max: float,
        imag_min: float,
        imag_max: float
    ):
        
        contour = RectangleContour(
            real_min, 
            real_max, 
            imag_min,
            imag_max
            )
        
        domains_number_of_roots = count_root_containing_domains(
            self.f, 
            self.df, 
            contour,
            count_roots_fn = count_roots_numerical,
            debug=True)
        
        all_roots = find_all_roots_in_domains(
            self.f,
            self.df,
            domains_number_of_roots,
            find_roots_fn=find_roots_delves_lynes,
            debug=True
        )
        
        return all_roots


def count_root_containing_domains(
    f: Callable[[np.ndarray], np.ndarray],
    df: Callable[[np.ndarray], np.ndarray],
    contour: ContourBase,
    count_roots_fn: Callable[..., float],
    n_divide: int = None,
    max_roots_per_domain: int = 3,
    debug: bool = False,
) -> List[Tuple[ContourBase, int]]:
    """
    Recursively subdivide a contour until each subdomain contains
    at most `max_roots_per_domain` roots, using a root-counting function.

    Parameters
    ----------
    f : Callable
        Scalar-valued function or determinant function.
    df : Callable
        Derivative of f. Will be passed to `count_roots_fn`.
    contour : ContourBase
        Initial domain (rectangle or circle).
    count_roots_fn : Callable
        Function that estimates the number of roots in a contour,
        such as `count_roots_numerical` or `count_roots_unity`.
    max_roots_per_domain : int, optional
        Maximum number of roots allowed in a final subdomain (default is 3).
    n_divide : int, optional
        How many times to divide along each axis/direction (default is governed by the Contour).
        For rectangles: creates n_divide^2 sub-rectangles.
    debug : bool, optional
        Whether to plot subdomains during recursion (default is False).

    Returns
    -------
    List[Tuple[ContourBase, int]]
        List of (contour, estimated_root_count) pairs.
    """
    # add original contour to deque (list that has a fast .pop() method)
    queue = deque([contour])
    final_domains: List[Tuple[ContourBase, int]] = []

    while queue:
        current_contour = queue.popleft()

        try:
            n_roots = count_roots_fn(f, current_contour, df=df)
        except Exception as e:
            print(f"[Warning] Failed to count roots in domain {current_contour}: {e}")
            continue

        if np.iscomplex(n_roots):
            # Handle branch cut effects heuristically
            if np.imag(n_roots) != 0:
                n_roots = np.real(n_roots) - 0.5 * np.imag(n_roots)

        n_roots = float(np.real(n_roots))  # ensure it's a float

        if np.abs(n_roots) < 1e-8:
            continue

        if n_roots <= max_roots_per_domain and not np.isclose(n_roots % 1, 0.5):
            final_domains.append((current_contour, int(round(n_roots))))
        else:
            children = current_contour.subdivide()
            queue.extend(children)

            if debug:
                try:
                    import matplotlib.pyplot as plt
                    for c in children:
                        Z = c.Z
                        plt.plot(np.real(Z), np.imag(Z), 'r-', alpha=0.5)
                    plt.pause(0.01)
                except ImportError:
                    print("[Debug] matplotlib not available for plotting.")

    return final_domains

def find_all_roots_in_domains(
    f: Callable[[np.ndarray], np.ndarray],
    df: Callable[[np.ndarray], np.ndarray],
    domains: List[Tuple[ContourBase, int]],
    find_roots_fn: Callable,
    polish: bool = True,
    polish_tol: float = 1e-10,
    merge_tol: float = 1e-8,
    debug: bool = False,
) -> np.ndarray:
    """
    Find and polish roots in multiple subdomains using a specified root-finding function.

    Parameters
    ----------
    f : Callable
        Function whose roots to find.
    df : Callable
        Derivative of f.
    domains : list of tuples (ContourBase, int)
        Each entry is a subdomain and the number of roots expected in it.
    find_roots_fn : Callable
        Function that estimates roots in a subdomain: 
        find_roots_fn(f, contour, n_roots, df, previous_roots=None)
    polish : bool
        Whether to apply Nelder-Mead to refine each root (default True).
    polish_tol : float
        Tolerance for polishing step (default 1e-10).
    merge_tol : float
        Tolerance for merging close roots (default 1e-8).
    debug : bool
        Whether to plot the found roots during the process.

    Returns
    -------
    np.ndarray
        Array of all unique (polished) roots found.
    """
    all_roots = []

    for contour, n_roots in domains:
        if n_roots <= 0:
            continue

        try:
            roots = find_roots_fn(f, contour, n_roots, df=df)
        except Exception as e:
            print(f"[Warning] Root finding failed in domain {contour}: {e}")
            continue

        roots = np.atleast_1d(roots)

        if polish:
            polished = []
            for z0 in roots:
                result = minimize(
                    lambda z: np.abs(f(z[0] + 1j * z[1]))**2,
                    x0=[np.real(z0), np.imag(z0)],
                    method="Nelder-Mead",
                    tol=polish_tol,
                    options={"disp": False}
                )
                z_star = result.x[0] + 1j * result.x[1]
                polished.append(z_star)
            roots = np.array(polished)

        all_roots.extend(roots)

        if debug:
            plt.scatter(np.real(roots), np.imag(roots), marker='x', color='black', s=50)

    if debug:
        plt.pause(0.01)  # allow figure to render interactively

    if len(all_roots) == 0:
        return np.array([])

    all_roots = np.array(all_roots)

    # Merge close roots (deduplication)
    dist_matrix = squareform(pdist(all_roots[:, np.newaxis].view(np.float64).reshape(-1, 2)))
    keep = np.ones(len(all_roots), dtype=bool)

    for i in range(len(all_roots)):
        if not keep[i]:
            continue
        keep[(dist_matrix[i] < merge_tol) & (np.arange(len(all_roots)) > i)] = False

    return all_roots[keep]


    
    