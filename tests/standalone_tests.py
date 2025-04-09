# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 15:20:28 2025

@author: RensvanLeijden
"""

# %% IMPORTS

import numpy as np
from typing import Callable, Dict, Tuple
from complex_root_finder.argument_principle import argument_principle
from complex_root_finder import count_roots_numerical, count_roots_unity, find_roots_delves_lynes, find_roots_austin_kravanja, RectangleContour, CircleContour

import time


# %%
# Test functions and their known roots
TEST_FUNCTIONS: Dict[str, Dict] = {
    'quadratic': {
        'func': lambda z: z**2 - 1,
        'derivative': lambda z: 2*z,
        'roots': [1, -1],
        'description': 'f(z) = z² - 1'
    },
    'cubic': {
        'func': lambda z: z**3 - 1,
        'derivative': lambda z: 3*z**2,
        'roots': [1, -0.5 + 0.866j, -0.5 - 0.866j],
        'description': 'f(z) = z³ - 1'
    },
    'quartic': {
        'func': lambda z: z**4 + 2*z**2 + 1,
        'derivative': lambda z: 4*z**3 + 4*z,
        'roots': [1j, - 1j, 1j, - 1j],
        'description': 'f(z) = z⁴ + 2z² + 1'
    },
    'quartic_minus_16': {
        'func': lambda z: z**4 - 16,
        'derivative': lambda z: 4*z**3,
        'roots': [2, -2, 2j, -2j],
        'description': 'f(z) = z⁴ - 16'
    }
}

def test_root_counting_method(
    method: Callable,
    func_name: str,
    bounds: Tuple[float, float, float, float],
    method_name: str
) -> Tuple[bool, float]:
    """
    Test a specific root-finding method on a given function.
    
    Parameters
    ----------
    method : Callable
        The root-finding method to test
    func_name : str
        Name of the test function to use
    bounds : Tuple[float, float, float, float]
        (real_min, real_max, imag_min, imag_max)
    method_name : str
        Name of the method for reporting
        
    Returns
    -------
    Tuple[bool, float]
        (passed/failed, execution time in seconds)
    """
    test_func = TEST_FUNCTIONS[func_name]
    expected_roots = len(test_func['roots'])
    
    print(f"\nTesting {method_name} on {test_func['description']}")
    
    # Start timing
    start_time = time.perf_counter()
    
    if method_name == 'argument_principle':
        n_roots = method(
            f=test_func['func'],
            real_min=bounds[0],
            real_max=bounds[1],
            imag_min=bounds[2],
            imag_max=bounds[3],
            step_size=0.1
        )
    elif method_name == 'symbolic':
        contour = RectangleContour(*bounds)
        n_roots = method(
            f=test_func['func'],
            contour=contour,
            df=test_func['derivative']
        )
    elif method_name == 'numerical':
        contour = RectangleContour(*bounds)
        n_roots = method(
            f=test_func['func'],
            contour=contour
        )
    elif method_name == 'unity':
        # For unity method, use a circle that contains all roots
        center = (bounds[0] + bounds[1])/2 + 1j*(bounds[2] + bounds[3])/2
        radius = max(bounds[1] - bounds[0], bounds[3] - bounds[2])/2
        contour = CircleContour(center,radius)
        n_roots = method(
            f=test_func['func'],
            contour=contour
        )
    
    # End timing
    execution_time = time.perf_counter() - start_time
    
    print(f"Expected roots: {expected_roots}")
    print(f"Found roots: {n_roots}")
    print(f"Execution time: {execution_time:.4f} seconds")
    passed = abs(n_roots - expected_roots) < 0.5
    print(f"Test {'passed' if passed else 'failed'}")
    
    return passed, execution_time

def test_root_finding_method(
    method: Callable,
    func_name: str,
    bounds: Tuple[float, float, float, float],
    method_name: str,
    tolerance: float = 1e-4
) -> Tuple[bool, float]:
    """
    Test a specific root-finding method on a given function.
    
    Parameters
    ----------
    method : Callable
        The root-finding method to test
    func_name : str
        Name of the test function to use
    bounds : Tuple[float, float, float, float]
        (real_min, real_max, imag_min, imag_max)
    method_name : str
        Name of the method for reporting
    tolerance : float
        Maximum allowed distance between found and expected roots
        
    Returns
    -------
    Tuple[bool, float]
        (passed/failed, execution time in seconds)
    """
    test_func = TEST_FUNCTIONS[func_name]
    expected_roots = np.array(test_func['roots'])
    
    print(f"\nTesting {method_name} on {test_func['description']}")
    
    # Start timing
    start_time = time.perf_counter()
    
    if method_name == 'delves_lynes_rectangle_df_known':
        contour = RectangleContour(*bounds)
        found_roots = method(
            f=test_func['func'],
            contour=contour,
            n_roots=len(test_func['roots']),
            df=test_func['derivative']
        )
    elif method_name == 'delves_lynes_circle_df_known':
        center = (bounds[0] + bounds[1])/2 + 1j*(bounds[2] + bounds[3])/2
        radius = max(bounds[1] - bounds[0], bounds[3] - bounds[2])/2
        contour = CircleContour(center,radius)
        found_roots = method(
            f=test_func['func'],
            contour=contour,
            n_roots=len(test_func['roots']),
            df=test_func['derivative']
        )
    elif method_name == 'delves_lynes_rectangle':
        contour = RectangleContour(*bounds)
        found_roots = method(
            f=test_func['func'],
            contour=contour,
            n_roots=len(test_func['roots']),
        )
    elif method_name == 'delves_lynes_circle':
        center = (bounds[0] + bounds[1])/2 + 1j*(bounds[2] + bounds[3])/2
        radius = max(bounds[1] - bounds[0], bounds[3] - bounds[2])/2
        contour = CircleContour(center,radius)
        found_roots = method(
            f=test_func['func'],
            contour=contour,
            n_roots=len(test_func['roots']),
        )
    elif method_name == 'austin_kravanja':
        # For unity method, use a circle that contains all roots
        center = (bounds[0] + bounds[1])/2 + 1j*(bounds[2] + bounds[3])/2
        radius = max(bounds[1] - bounds[0], bounds[3] - bounds[2])/2
        contour = CircleContour(center,radius)
        found_roots = method(
            f=test_func['func'],
            contour=contour,
            n_roots=len(test_func['roots']),
        )

    # End timing
    execution_time = time.perf_counter() - start_time

    if found_roots is None or len(found_roots) == 0:
        print("No roots found!")
        print(f"Execution time: {execution_time:.4f} seconds")
        return False, execution_time
    
    found_roots = np.array(found_roots)
    
    # Check if we found the correct number of roots
    if len(found_roots) != len(expected_roots):
        print(f"Wrong number of roots! Expected {len(expected_roots)}, found {len(found_roots)}")
        return False, execution_time
    
    # For each expected root, find the closest found root and check distance
    max_error = 0 
    for expected in expected_roots:
        distances = np.abs(found_roots - expected)
        min_distance = np.min(distances)
        max_error = max(max_error, min_distance)
    
    passed = max_error < tolerance
    
    print(f"Expected roots: {expected_roots}")
    print(f"Found roots: {found_roots}")
    print(f"Maximum error: {max_error}")
    print(f"Execution time: {execution_time:.4f} seconds")
    print(f"Test {'passed' if passed else 'failed'}")
    
    return passed, execution_time

def run_all_tests(configurations, methods, test_function):
    """Run all tests for all methods on all test functions."""
     
    # Run all tests
    results = {}
    timing_results = {}
    for method_name, method in methods.items():
        results[method_name] = {}
        timing_results[method_name] = {}
        for func_name, bounds in configurations.items():
            passed, execution_time = test_function(
                method, func_name, bounds, method_name
            )
            results[method_name][func_name] = passed
            timing_results[method_name][func_name] = execution_time
    
    # Print summary
    print("\nTest Summary:")
    print("-" * 50)
    for method_name, method_results in results.items():
        passed = sum(method_results.values())
        total = len(method_results)
        avg_time = np.mean(list(timing_results[method_name].values()))
        max_time = max(timing_results[method_name].values())
        print(f"{method_name}:")
        print(f"  Passed: {passed}/{total} tests")
        print(f"  Average time: {avg_time:.4f} seconds")
        print(f"  Maximum time: {max_time:.4f} seconds")
        print()

if __name__ == "__main__":
    
    # Test configurations
    test_configs = {
        'quadratic': (-2, 2, -2, 2),
        'cubic': (-2, 2, -2, 2),
        'quartic': (-2, 2, -2, 2),
        'quartic_minus_16': (-4, 4, -4, 4)
    }
    # test_configs = {
    #     'quartic_minus_16': (-4, 4, -4, 4)
    # }
    
    # Methods to test root finding
    methods_root_counting = {
        'argument_principle': argument_principle,
        'symbolic': count_roots_numerical,
        'numerical': count_roots_numerical,
        'unity': count_roots_unity
    }
    
    methods_root_finding = {
        'delves_lynes_rectangle_df_known': find_roots_delves_lynes,
        'delves_lynes_circle_df_known': find_roots_delves_lynes,
        'delves_lynes_rectangle': find_roots_delves_lynes,
        'delves_lynes_circle': find_roots_delves_lynes,
        'austin_kravanja': find_roots_austin_kravanja,
    }
    
    methods_root_finding = {
        'delves_lynes_rectangle': find_roots_delves_lynes,
        'austin_kravanja': find_roots_austin_kravanja,
    }

    # Run root counting tests
    print("Testing root counting methods...")
    print("=" * 50)
    run_all_tests(test_configs, methods_root_counting, test_root_counting_method)
    
    # Run root finding tests
    print("\nTesting root finding methods...")
    print("=" * 50)
    run_all_tests(test_configs, methods_root_finding, test_root_finding_method)
    
    