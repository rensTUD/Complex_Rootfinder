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
from matplotlib.path import Path

# imports from current library
from .utils import finite_difference_first_derivative_4th_order, partition_contour_along_cut_general, is_point_inside_contour, wrap_contour_around_branch_point, plot_contour
from .contours import ContourBase, RectangleContour, CircleContour, BranchCut, CompositeContour
from .count_roots import count_roots_numerical, count_roots_unity
from .find_roots import find_roots_delves_lynes, find_roots_austin_kravanja
from .plot_utils import debug_plot_contours

# %% CLASS DEFINITION


class ComplexRootFinder:
    def __init__(
        self,
        f: Callable[[np.ndarray], np.ndarray],
        df: Callable[[np.ndarray], np.ndarray] | None = None,
        branch_cuts: List[BranchCut] | None = None
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
        self.f = f
        
        if df is None:
            self.df = lambda z: finite_difference_first_derivative_4th_order(self.f, z)  # Use finite difference if no derivative provided
        else:
            self.df = df
            
        # assign branch cuts
        if len(branch_cuts) > 1:
            self.branch_cuts = branch_cuts
        else:
            self.branch_cuts = [branch_cuts]
                    
    def find_roots_domain(
        self,
        contour: ContourBase,
        n_divide: int = None,
        max_roots_per_domain: int = 3,
        max_depth: int = 5,
        polish: bool = True,
        polish_tol: float = 1e-10,
        merge_tol: float = 1e-8,
        debug: bool = False,
    ):
        
        # count roots
        domains_number_of_roots = self.count_root_containing_domains(
            contour=contour,
            count_roots_fn = count_roots_numerical,
            n_divide = n_divide,
            max_roots_per_domain = max_roots_per_domain,
            max_depth = max_depth,
            debug=debug)
        
        # find roots
        all_roots = self.find_all_roots_in_domains(
            original_contour = contour,
            domains = domains_number_of_roots,
            find_roots_fn = find_roots_delves_lynes,
            polish = polish,
            polish_tol = polish_tol,
            merge_tol = merge_tol,
            debug = debug,
        )
        
        return all_roots
        

    def count_root_containing_domains(
        self,
        contour: ContourBase,
        count_roots_fn: Callable[..., float] = count_roots_numerical,
        n_divide: int = None,
        max_roots_per_domain: int = 3,
        max_depth: int = 5,
        debug: bool = False,
    ) -> List[Tuple[ContourBase, int]]:
        """
        Recursively subdivide a contour until each subdomain contains
        at most `max_roots_per_domain` roots, using a root-counting function.
    
        Parameters
        ----------
        contour : ContourBase
            Initial domain (rectangle or circle).
        count_roots_fn : Callable
            Function that estimates the number of roots in a contour,
            such as `count_roots_numerical` or `count_roots_unity`.
        max_roots_per_domain : int, optional
            Maximum number of roots allowed in a final subdomain (default is 3).
        n_divide : int, optional
            How many times to divide along each axis/direction (default is governed by the Contour).
        max_depth : int, optional
            Maximum number of recursive subdivisions allowed (default is 10).
        debug : bool, optional
            Whether to plot subdomains during recursion (default is False).
    
        Returns
        -------
        List[Tuple[ContourBase, int]]
            List of (contour, estimated_root_count) pairs.
        """
        
        queue = deque([(contour, 0)])  # Each entry is (contour, depth)
        final_domains: List[Tuple[ContourBase, int]] = []
    
        while queue:
            current_contour, depth = queue.popleft()
    
            if depth >= max_depth:
                if debug:
                    print(f"[Max Depth Reached] Skipping subdivision beyond depth {max_depth}")
                continue
    
            try:
                n_roots = count_roots_fn(self.f, current_contour, df=self.df)
            except Exception as e:
                print(f"[Warning] Failed to count roots in domain {current_contour}: {e}")
                continue
    
            # Subdivide if branch cut is detected
            if np.iscomplex(n_roots) and not np.isclose(np.imag(n_roots), 0):
                children = current_contour.subdivide(n_divide=n_divide) if n_divide else current_contour.subdivide()
                queue.extend((child, depth + 1) for child in children)
    
                if debug:
                    print(f"[Branch Cut] Subdividing due to Im(n_roots) = {np.imag(n_roots):.3f}")
                    debug_plot_contours(children, color='orange')
    
                continue
    
            n_roots = float(np.real(n_roots))
    
            if np.abs(n_roots) < 1e-8:
                continue
    
            if n_roots <= max_roots_per_domain and not np.isclose(n_roots % 1, 0.5):
                final_domains.append((current_contour, int(round(n_roots))))
                if debug:
                    print(f"{int(round(n_roots))} number of roots found!")
            else:
                children = current_contour.subdivide(n_divide=n_divide) if n_divide else current_contour.subdivide()
                queue.extend((child, depth + 1) for child in children)
    
                if debug:
                    debug_plot_contours(children, color='red')
    
        return final_domains
    
    def count_root_containing_domains_branch_cut(
        self,
        contour,
        count_roots_fn: Callable[..., float] | None = None,
        n_divide: int = None,
        max_roots_per_domain: int = 3,
        max_depth: int = 5,
        debug: bool = False,
    ) -> List[Tuple["CompositeContour", int]]:
        """
        Recursively subdivide a contour until each subdomain contains
        at most `max_roots_per_domain` roots, taking into account
        branch cuts and branch points.
    
        Parameters
        ----------
        contour : ContourBase or CompositeContour
            The initial contour (rectangular, circular, or custom).
        count_roots_fn : Callable
            Root-counting function to use. Must be provided.
        n_divide : int
            Number of divisions when subdividing.
        max_roots_per_domain : int
            Maximum roots allowed in a single domain.
        max_depth : int
            Maximum recursion depth.
        debug : bool
            Enable debug printing and plotting.
    
        Returns
        -------
        List[Tuple[CompositeContour, int]]
            List of final (contour, root_count) pairs.
        """
        if count_roots_fn is None:
            raise ValueError("A root-counting function must be provided.")
    
        queue = deque([(contour, 0)])
        final_domains = []
    
        while queue:
            current_contour, depth = queue.popleft()
            Z = current_contour.Z 
    
            if depth >= max_depth:
                if debug:
                    print(f"[Max Depth Reached] Skipping subdivision beyond depth {max_depth}")
                continue
                    
            # Check for branch cut intersections or branch points
            split_due_to_cut = False
            if self.branch_cuts:
                for cut in self.branch_cuts:
                    # check for:
                    branch_cut_intersects = cut.intersects(current_contour)
                    if cut.branch_point:
                        branch_point_in_contour = is_point_inside_contour(Z, cut.branch_point)
                    else:
                        branch_point_in_contour = None
                    if branch_cut_intersects and not branch_point_in_contour:
                        if debug:
                            print("[Branch Cut] Splitting domain along cut")
                            plot_contour(Z, color='orange')
                        subpaths = partition_contour_along_cut_general(Z, cut.points)
                        for path in subpaths:
                            queue.append((CompositeContour([path]), depth + 1))
                        split_due_to_cut = True
                        break
                    elif branch_cut_intersects and branch_point_in_contour:
                        if debug:
                            print("[Branch Point] Wrapping around branch point")
                            plot_contour(Z, color='blue')
                        wrapped = wrap_contour_around_branch_point(current_contour, cut)
                        queue.append((wrapped, depth + 1))
                        split_due_to_cut = True
                        break
            if split_due_to_cut:
                continue
    
            # Root counting
            try:
                n_roots = count_roots_fn(self.f, current_contour, df=self.df)
            except Exception as e:
                print(f"[Warning] Failed to count roots in domain {current_contour}: {e}")
                continue
            
            # complex root count, meaning we are on a branch cut, or something else went wrong: subdivide
            if np.iscomplex(n_roots) and not np.isclose(np.imag(n_roots), 0):
                if debug:
                    print("[Complex Root Count] Subdividing due to ambiguity")
                    plot_contour(Z, color='gray')
                children = current_contour.subdivide() if hasattr(current_contour, "subdivide") else []
                queue.extend((child, depth + 1) for child in children)
                continue
            
            # if 0 roots are found, continue
            n_roots = float(np.real(n_roots))
            if np.abs(n_roots) < 1e-8:
                continue
            
            #
            if n_roots <= max_roots_per_domain and not np.isclose(n_roots % 1, 0.5):
                if not isinstance(current_contour, CompositeContour):
                    current_contour = CompositeContour([Z])
                final_domains.append((current_contour, int(round(n_roots))))
                if debug:
                    print(f"[Accepted] {int(round(n_roots))} root(s) in contour.")
                    plot_contour(Z, color='green')
            else:
                if debug:
                    print("[Subdivision] Too many or ambiguous roots")
                    plot_contour(Z, color='red')
                children = current_contour.subdivide() if hasattr(current_contour, "subdivide") else []
                queue.extend((child, depth + 1) for child in children)
    
        return final_domains
    
    def find_all_roots_in_domains(
        self,
        original_contour: ContourBase,
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
        original_contour : ContourBase
            original contour that we are looking in
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
                roots = find_roots_fn(self.f, contour, n_roots, df=self.df)
            except Exception as e:
                print(f"[Warning] Root finding failed in domain {contour}: {e}")
                continue
    
            roots = np.atleast_1d(roots)
            accepted = []
    
            if polish:
                # Build path for domain check
                Z_poly = np.column_stack((np.real(contour.Z), np.imag(contour.Z)))
                path = Path(Z_poly)
                for z0 in roots:
                    if not path.contains_point((np.real(z0), np.imag(z0)), radius=-1e-5):
                        z0 = contour.center
                    result = minimize(
                        lambda z: np.abs(self.f(z[0] + 1j * z[1])),
                        x0=[np.real(z0), np.imag(z0)],
                        method='L-BFGS-B',
                        bounds = contour.bounds,
                        tol=polish_tol,
                        options={"disp": False}
                    )
                    x, y = result.x
                    z_star = x + 1j * y
                    
                    if path.contains_point((x, y), radius=-1e-5):
                        accepted.append(z_star)
                    else:
                        if debug:
                            print(f"[Polish Rejected] Root {z_star:.4f} outside domain — skipped.")
            else:
                accepted = roots
    
            if debug and accepted:
                plt.scatter(
                    np.real(accepted),
                    np.imag(accepted),
                    s=50,
                    c='green',
                    marker='x',
                    linewidths=1.5,
                    label='Polished roots'
                )
    
            all_roots.extend(accepted)
    
        if debug:
            plt.pause(0.01)
    
        if len(all_roots) == 0:
            return np.array([])
    
        all_roots = np.array(all_roots)
    
        # Deduplicate roots using pairwise distance
        dist_matrix = squareform(pdist(all_roots[:, np.newaxis].view(np.float64).reshape(-1, 2)))
        keep = np.ones(len(all_roots), dtype=bool)
    
        for i in range(len(all_roots)):
            if not keep[i]:
                continue
            keep[(dist_matrix[i] < merge_tol) & (np.arange(len(all_roots)) > i)] = False
    
        return all_roots[keep]

    
    