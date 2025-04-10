# -*- coding: utf-8 -*-
"""
Parametric root finding over omega using ComplexRootFinder
"""

# %% Imports
import numpy as np
from complex_root_finder import ComplexRootFinder, RectangleContour
from complex_root_finder.parametric_root_finder import ParametricRootFinder

# %% Define physical constants and param sweep

c_p = 1.719084289205639e+03 + 3.322560356224935e+01j
c_s = 3.704985415873700e+02 + 1.875674989161015e+01j
Lambda = 5.114035573342751e+09 + 1.649233839281022e+08j
mu = 2.612383107677712e+08 + 2.651871379876653e+07j

# Sweep over omega
omega_list = np.linspace(10, 300, 20) * 2 * np.pi  # 20 frequencies from 200–800 Hz

# %% Define f(z, omega)

# # Set elastic medium parameters
# E = 7e7
# nu = 0.4
# rho = 1700

# # Lamé parameters
# Lambda = E*nu/((1+nu)*(1-2*nu))
# mu = E/(2*(1+nu))

# # Wave speeds
# c_p = np.sqrt((Lambda+2*mu)/rho)
# c_s = np.sqrt(mu/rho)

# # Attenuation
# alpha_p = 1
# alpha_s = 2
# eta = (40*np.pi*np.log10(np.e))**(-1)

# c_p = np.sqrt((Lambda+2*mu)/rho)
# c_s = np.sqrt(mu/rho)
# c_p = c_p/(1-1j*eta*alpha_p)
# c_s = c_s/(1-1j*eta*alpha_s)


def f_omega(k_x: np.ndarray, omega: float) -> np.ndarray:

    
    H = 10

    f = 32 * ((((-(0.1e1 / 0.16e2) - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 16 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 16 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 16) * mu + (-(0.1e1 / 0.32e2) - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 32 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 32 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 32) * Lambda) * (omega ** 4) - (-(0.1e1 / 0.4e1) - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 4 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 4 - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 4 + np.exp(-1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2)))) * (mu * (c_p ** 2) + 2 * (c_s ** 2) * (mu + Lambda / 2)) * (k_x ** 2) * (omega ** 2) / 4 + (-(0.1e1 / 0.4e1) - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 4 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 4 - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 4 + np.exp(-1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2)))) * (c_p ** 2) * (c_s ** 2) * mu * (k_x ** 4)) * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + (np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) + np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) - 1) * ((omega ** 4) * ((0.3e1 / 0.4e1) * mu + Lambda / 8) - (0.3e1 / 0.4e1) * (mu * (c_p ** 2) + (0.4e1 / 0.3e1) * (c_s ** 2) * (mu + Lambda / 4)) * (k_x ** 2) * (omega ** 2) + (c_p ** 2) * (c_s ** 2) * (k_x ** 4) * mu) * (k_x ** 2) / 4) * mu / (c_p ** 2) / (c_s ** 2)

    return f

# %% Function generator

def f_generator(omega):
    return lambda z: f_omega(z, omega)

def df_generator(omega):
    return None  # Use automatic finite difference

# %% Contour generator

def contour_generator(omega):
    k_s = omega / c_s
    real_max = np.max(np.real(k_s)) / 0.2
    real_min = -real_max
    imag_min = -3
    imag_max = 5e-5
    return RectangleContour(real_min, real_max, imag_min, imag_max)

# %% Run parallel root finder

sweep = ParametricRootFinder(
    param_list=omega_list,
    f_generator=f_generator,
    df_generator=df_generator,
    contour_generator=contour_generator,
)

# Use all CPU cores except 2
roots_dict = sweep.run_all(
    max_roots_per_domain=3,
    max_depth=5,
    polish=True,
    n_jobs=-2,
    debug=False,
)

# %% Print results

for omega, roots in roots_dict.items():
    print(f"\nω = {omega / (2*np.pi):.2f} Hz → {len(roots)} roots:")
    for r in roots:
        print(f"  {r:.6f}")
