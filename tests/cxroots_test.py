# %% IMPORTS
import numpy as np
from cxroots import Rectangle, Circle
from complex_root_finder.argument_principle import argument_principle
# %% 
# Material properties
c_p = 1.719084289205639e3 + 3.322560356224935e1j
c_s = 3.704985415873700e2 + 1.875674989161015e1j
lambda_ = 5.114035573342751e9 + 1.649233839281022e8j
mu = 2.612383107677712e8 + 2.651871379876653e7j

# Geometry and excitation
omega = 10 * 2 * np.pi
H = 10

# Define the function
def f(k_x):
    # Helper terms
    A = np.sqrt((omega**2 - (k_x**2) * c_p**2) / c_p**2)
    B = np.sqrt((omega**2 - (k_x**2) * c_s**2) / c_s**2)

    exp_AB = np.exp(-2j * H * (A + B))
    exp_2A = np.exp(-2j * A * H)
    exp_2B = np.exp(-2j * B * H)
    exp_sum = np.exp(-1j * H * (A + B))

    term1 = ((-1/16 - exp_AB / 16 - exp_2A / 16 - exp_2B / 16) * mu +
             (-1/32 - exp_AB / 32 - exp_2A / 32 - exp_2B / 32) * lambda_) * omega**4

    term2 = (-1/4 - exp_2A / 4 - exp_2B / 4 - exp_AB / 4 + exp_sum) \
            * (mu * c_p**2 + 2 * c_s**2 * (mu + lambda_ / 2)) * k_x**2 * omega**2 / 4

    term3 = (-1/4 - exp_2A / 4 - exp_2B / 4 - exp_AB / 4 + exp_sum) \
            * c_p**2 * c_s**2 * mu * k_x**4

    term4 = (exp_2A + exp_2B - exp_AB - 1) * (
        omega**4 * (3/4 * mu + lambda_ / 8)
        - 3/4 * (mu * c_p**2 + 4/3 * c_s**2 * (mu + lambda_ / 4)) * k_x**2 * omega**2
        + c_p**2 * c_s**2 * k_x**4 * mu
    ) * k_x**2 / 4

    return 32 * mu / (c_p**2 * c_s**2) * (term1 + term2 + term3 + term4)



# %%
# rect = Rectangle(x_range=(-5, 5), y_range=(-5, 1))

# rect.roots(f)




# %%
# circ = Circle(0, 5)

# circ.roots(f)

# %%
real_min= -5
real_max= 5
imag_min= -5
imag_max= 1

n_roots = argument_principle(
            f=f,
            real_min=real_min,
            real_max=real_max,
            imag_min=imag_min,
            imag_max=imag_max,
            step_size=0.1 
        )
# %%
