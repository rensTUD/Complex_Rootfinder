# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 19:42:28 2024

@author: RensvanLeijden
"""

# %% import packages
from typing import Callable

# %% functions

def argument_principle(
    f: Callable[[complex], complex],
    real_min: float, 
    real_max: float, 
    imag_min: float, 
    imag_max: float, 
    step_size: float
) -> int:
    """
    Find the number of roots of a complex analytic function using the argument principle.

    The function traverses a rectangular boundary in the complex plane and counts the number
    of times the function's argument changes, which corresponds to the number of roots inside
    the boundary.

    Parameters
    ----------
    f : Callable[[complex], complex]
        Complex analytic function to analyze.
    real_min : float
        Minimum real value of the boundary rectangle.
    real_max : float
        Maximum real value of the boundary rectangle.
    imag_min : float
        Minimum imaginary value of the boundary rectangle.
    imag_max : float
        Maximum imaginary value of the boundary rectangle.
    step_size : float
        Initial step size for traversal along the boundary.

    Returns
    -------
    int
        Number of roots inside the given boundary. Returns special values:
        - 2001: Step limit reached on top boundary
        - 2002: Step limit reached on right boundary
        - 2003: Step limit reached on bottom boundary
        - 2004: Step limit reached on left boundary

    Notes
    -----
    The function uses the argument principle to count roots by tracking changes in the
    function's argument as it traverses the boundary. The step size is automatically
    reduced when necessary to ensure accurate counting.
    """
    # Constants and Limits
    step_limit = 1e-12  # Limit to prevent getting stuck in step reduction

    # Initialize variables
    calc_counter = 0
    n_singularities_unknown = True
    n_singularities_previous = 0
    n_count = 0
    n_roots = 0
    step = step_size

    # Start Argument Principle Calculation
    while n_singularities_unknown:
        calc_counter += 1
        n_count = 0
        original_step = step

        # First Part: Horizontal Path from real_min to real_max (top boundary)
        current_real_value = f(complex(real_min, imag_max)).real
        current_imag_value = f(complex(real_min, imag_max)).imag

        real_pos = real_min

        while real_pos < real_max:
            if real_pos + step < real_max:
                next_real_pos = real_pos + step
                next_wc = complex(next_real_pos, imag_max)
            else:
                next_real_pos = real_max
                step = real_max - real_pos
                next_wc = complex(real_max, imag_max)

            next_value = f(next_wc)
            next_real_value = next_value.real
            next_imag_value = next_value.imag

            if current_real_value * next_real_value < 0 and current_imag_value * next_imag_value < 0:
                step = 0.3 * step
                if step < step_limit:
                    print(f"Warning: Step limit reached without converging, might be stuck at: ({next_wc.real}, {next_wc.imag})")
                    n_roots = 2001  # Indicator for a potential issue
                    return n_roots
            else:
                n_count = update_n_count(n_count, current_real_value, current_imag_value, next_real_value, next_imag_value)
                current_real_value = next_real_value
                current_imag_value = next_imag_value
                step = original_step
                real_pos = next_real_pos

        # Second Part: Vertical Path from imag_max to imag_min (right boundary)
        current_real_value = f(complex(real_max, imag_max)).real
        current_imag_value = f(complex(real_max, imag_max)).imag

        imag_pos = imag_max

        while imag_pos > imag_min:
            if imag_pos - step > imag_min:
                next_imag_pos = imag_pos - step
                next_wc = complex(real_max, next_imag_pos)
            else:
                next_imag_pos = imag_min
                step = imag_pos - imag_min
                next_wc = complex(real_max, imag_min)

            next_value = f(next_wc)
            next_real_value = next_value.real
            next_imag_value = next_value.imag

            if current_real_value * next_real_value < 0 and current_imag_value * next_imag_value < 0:
                step = 0.3 * step
                if step < step_limit:
                    print(f"Warning: Step limit reached without converging, might be stuck at: ({next_wc.real}, {next_wc.imag})")
                    n_roots = 2002  # Indicator for a potential issue
                    return n_roots
            else:
                n_count = update_n_count(n_count, current_real_value, current_imag_value, next_real_value, next_imag_value)
                current_real_value = next_real_value
                current_imag_value = next_imag_value
                step = original_step
                imag_pos = next_imag_pos

        # Third Part: Horizontal Path from real_max to real_min (bottom boundary)
        current_real_value = f(complex(real_max, imag_min)).real
        current_imag_value = f(complex(real_max, imag_min)).imag

        real_pos = real_max

        while real_pos > real_min:
            if real_pos - step > real_min:
                next_real_pos = real_pos - step
                next_wc = complex(next_real_pos, imag_min)
            else:
                next_real_pos = real_min
                step = real_pos - real_min
                next_wc = complex(real_min, imag_min)

            next_value = f(next_wc)
            next_real_value = next_value.real
            next_imag_value = next_value.imag

            if current_real_value * next_real_value < 0 and current_imag_value * next_imag_value < 0:
                step = 0.3 * step
                if step < step_limit:
                    print(f"Warning: Step limit reached without converging, might be stuck at: ({next_wc.real}, {next_wc.imag})")
                    n_roots = 2003  # Indicator for a potential issue
                    return n_roots
            else:
                n_count = update_n_count(n_count, current_real_value, current_imag_value, next_real_value, next_imag_value)
                current_real_value = next_real_value
                current_imag_value = next_imag_value
                step = original_step
                real_pos = next_real_pos

        # Fourth Part: Vertical Path from imag_min to imag_max (left boundary)
        current_real_value = f(complex(real_min, imag_min)).real
        current_imag_value = f(complex(real_min, imag_min)).imag

        imag_pos = imag_min

        while imag_pos < imag_max:
            if imag_pos + step < imag_max:
                next_imag_pos = imag_pos + step
                next_wc = complex(real_min, next_imag_pos)
            else:
                next_imag_pos = imag_max
                step = imag_max - imag_pos
                next_wc = complex(real_min, imag_max)

            next_value = f(next_wc)
            next_real_value = next_value.real
            next_imag_value = next_value.imag

            if current_real_value * next_real_value < 0 and current_imag_value * next_imag_value < 0:
                step = 0.3 * step
                if step < step_limit:
                    print(f"Warning: Step limit reached without converging, might be stuck at: ({next_wc.real}, {next_wc.imag})")
                    n_roots = 2004  # Indicator for a potential issue
                    return n_roots
            else:
                n_count = update_n_count(n_count, current_real_value, current_imag_value, next_real_value, next_imag_value)
                current_real_value = next_real_value
                current_imag_value = next_imag_value
                step = original_step
                imag_pos = next_imag_pos

        # Determine the number of roots based on the count of argument changes
        n_roots = abs(n_count // 4)
        if calc_counter == 1:
            step = step * 0.5
            n_singularities_previous = n_roots
        elif n_roots != n_singularities_previous:
            step = step * 0.5
            n_singularities_previous = n_roots
        else:
            n_singularities_unknown = False

    return n_roots

def update_n_count(n_count, real_A, imag_A, real_B, imag_B):
    """
    update_n_count: Updates the count of changes in argument for the root-finding algorithm.
    
    Parameters:
        n_count : int
            Current count of argument changes.
        real_A, imag_A : float
            Real and imaginary parts of the function value at the current position.
        real_B, imag_B : float
            Real and imaginary parts of the function value at the next position.
    
    Returns:
        n_count : int
            Updated count of argument changes.
    """
    if (real_A > 0 and imag_A < 0 and real_B > 0 and imag_B > 0):
        n_count += 1
    elif (real_A > 0 and imag_A > 0 and real_B < 0 and imag_B > 0):
        n_count += 1
    elif (real_A < 0 and imag_A > 0 and real_B < 0 and imag_B < 0):
        n_count += 1
    elif (real_A < 0 and imag_A < 0 and real_B > 0 and imag_B < 0):
        n_count += 1
    elif (real_A > 0 and imag_A > 0 and real_B > 0 and imag_B < 0):
        n_count -= 1
    elif (real_A < 0 and imag_A > 0 and real_B > 0 and imag_B > 0):
        n_count -= 1
    elif (real_A < 0 and imag_A < 0 and real_B < 0 and imag_B > 0):
        n_count -= 1
    elif (real_A > 0 and imag_A < 0 and real_B < 0 and imag_B < 0):
        n_count -= 1

    return n_count
 

