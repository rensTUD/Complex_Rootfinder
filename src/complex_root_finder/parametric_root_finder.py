# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 15:28:47 2025

@author: rensv
"""

# %% IMPORTS

from typing import Callable, List, Dict
import numpy as np
from joblib import Parallel, delayed
from multiprocessing import cpu_count
from tqdm import tqdm

from .complex_root_finder import ComplexRootFinder
from .count_roots import count_roots_numerical
from .find_roots import find_roots_delves_lynes
from .contours import ContourBase



# %% DEFINITIONS

class ParametricRootFinder:
    def __init__(
        self,
        param_list: List[float],
        f_generator: Callable[[float], Callable[[np.ndarray], np.ndarray]],
        df_generator: Callable[[float], Callable[[np.ndarray], np.ndarray]] | None,
        contour_generator: Callable[[float], ContourBase],
        count_method: Callable = count_roots_numerical,
        root_method: Callable = find_roots_delves_lynes,
    ):
        """
        Initialize a ParametricRootFinder.

        Parameters
        ----------
        param_list : List[float]
            List of parameter values (e.g. omega).
        f_generator : Callable
            Function f_generator(param) -> f(z), the function to analyze.
        df_generator : Callable or None
            Function df_generator(param) -> df(z), the derivative of f(z).
        contour_generator : Callable
            Function contour_generator(param) -> ContourBase
        count_method : Callable
            Root-counting function (default: count_roots_numerical).
        root_method : Callable
            Root-finding function (default: find_roots_delves_lynes).
        """
        self.param_list = param_list
        self.f_generator = f_generator
        self.df_generator = df_generator
        self.contour_generator = contour_generator
        self.count_method = count_method
        self.root_method = root_method

    def run_all(
        self,
        max_depth: int = 5,
        n_divide: int | None = None,
        max_roots_per_domain: int = 3,
        polish: bool = True,
        polish_tol: float = 1e-10,
        merge_tol: float = 1e-8,
        n_jobs: int | str = 1,
        debug: bool = False,
    ) -> Dict[float, np.ndarray]:
        """
        Run root-finding for all parameters (parallel if desired).
    
        Parameters
        ----------
        n_jobs : int or str
            Number of parallel jobs:
            - "all" to use all physical cores
            - negative to subtract from total cores (e.g. -2 uses total - 2)
            - default is 1 (serial)
    
        Other parameters: see class docstring.
    
        Returns
        -------
        Dict[float, np.ndarray]
            Dictionary mapping each param to its root array.
        """
        # Determine number of jobs
        total_cores = cpu_count()
        if n_jobs == "all":
            n_jobs_resolved = total_cores
        elif isinstance(n_jobs, int) and n_jobs < 0:
            n_jobs_resolved = max(1, total_cores + n_jobs)
        else:
            n_jobs_resolved = int(n_jobs)
    
        def process_param(iN):
            param = self.param_list[iN]
            f = self.f_generator(param)
            df = self.df_generator(param) if self.df_generator else None
            contour = self.contour_generator(param)
    
            finder = ComplexRootFinder(f, df)
            roots = finder.find_roots_domain(
                contour=contour,
                n_divide=n_divide,
                max_roots_per_domain=max_roots_per_domain,
                max_depth=max_depth,
                polish=polish,
                polish_tol=polish_tol,
                merge_tol=merge_tol,
                debug=debug,
            )
    
            return param, roots
    
        # with Parallel(n_jobs=n_jobs_resolved) as parallel:
        #     results = parallel(
        #         delayed(process_param)(param) for param in tqdm(self.param_list, desc="Root Finding")
        #     )
        
        with Parallel(n_jobs=n_jobs_resolved) as parallel:
            results = parallel(delayed(process_param)(iN) for iN in tqdm(range(len(self.param_list))))
            
            
        return dict(results)