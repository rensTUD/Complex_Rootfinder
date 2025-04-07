from typing import Callable, Tuple
import numpy as np
import sys
import os

from complex_root_finder.argument_principle import argument_principle


def test_polynomial_roots():
    """Test the argument principle with different polynomial functions."""
    
    def quadratic(z: complex) -> complex:
        """Test function: f(z) = z² - 1
        Has roots at z = 1 and z = -1
        """
        return z**2 - 1
    
    def cubic(z: complex) -> complex:
        """Test function: f(z) = z³ - 1
        Has roots at z = 1, -0.5 ± 0.866i
        """
        return z**3 - 1
    
    def quartic(z: complex) -> complex:
        """Test function: f(z) = z⁴ + 2z² + 1
        Has roots at z = ±(1±i)
        """
        return z**4 + 2*z**2 + 1

    def quartic_minus_16(z: complex) -> complex:
        """Test function: f(z) = z⁴ - 16
        Has roots at z = ±2 and z = ±2i
        """
        return z**4 - 16

    # Test cases setup
    test_cases = [
        {
            'name': 'Quadratic z² - 1',
            'func': quadratic,
            'bounds': (-2, 2, -2, 2),
            'expected_roots': 2
        },
        {
            'name': 'Cubic z³ - 1',
            'func': cubic,
            'bounds': (-2, 2, -2, 2),
            'expected_roots': 3
        },
        {
            'name': 'Quartic z⁴ + 2z² + 1',
            'func': quartic,
            'bounds': (-2, 2, -2, 2),
            'expected_roots': 4
        },
        {
            'name': 'Quartic z⁴ - 16',
            'func': quartic_minus_16,
            'bounds': (-4, 4, -4, 4),  # Increased bounds to capture all roots at ±2 and ±2i
            'expected_roots': 4
        }
    ]

    # Run tests
    for test in test_cases:
        print(f"\nTesting {test['name']}")
        real_min, real_max, imag_min, imag_max = test['bounds']
        n_roots = argument_principle(
            real_min=real_min,
            real_max=real_max,
            imag_min=imag_min,
            imag_max=imag_max,
            step_size=0.1   ,
            det_func=test['func']
        )
        print(f"Expected roots: {test['expected_roots']}")
        print(f"Found roots: {n_roots}")
        print(f"Test {'passed' if n_roots == test['expected_roots'] else 'failed'}")

if __name__ == "__main__":
    test_polynomial_roots() 