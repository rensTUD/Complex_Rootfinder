# -*- coding: utf-8 -*-
"""
Created on Mon Apr 9 21:45:28 2025

@author: RensvanLeijden
"""

# %% IMPORTS

import numpy as np
import matplotlib.pyplot as plt
from complex_root_finder import ComplexRootFinder, RectangleContour, CircleContour

# %% INPUTS

# Set elastic medium parameters
E = 7e7
nu = 0.4
rho = 1700

# Lamé parameters
Lambda = E*nu/((1+nu)*(1-2*nu))
mu = E/(2*(1+nu))

# Wave speeds
c_L = np.sqrt((Lambda+2*mu)/rho)
c_T = np.sqrt(mu/rho)

# Attenuation
alpha_L = 1
alpha_T = 2
eta = (40*np.pi*np.log10(np.e))**(-1)
c_L = c_L/(1-1j*eta*alpha_L)
c_T = c_T/(1-1j*eta*alpha_T)

# Set frequency and search domain
freq = 100
omega = freq*2*np.pi

# medium wavenumbers
k_p = omega/c_L
k_s = omega/c_T

# Depth of the elastic layer
Z_1 = 20
Nz = 200
z = np.linspace(0, Z_1, Nz)

# search domain
dx = 5e-3
real_min = -np.max(np.real(k_s)) / 0.2
real_max = np.max(np.real(k_s)) / 0.2
imag_min = -1
imag_max = 50e-6  

def f(k_r):
    """
    Calculate determinant using matrix method
    
    Parameters:
    -----------
    k_r : complex
        Radial wave number
    omega : float
        Angular frequency
    c_L, c_T : complex
        Longitudinal and transverse wave speeds
    Lambda, mu : float
        Lamé parameters
    Z_1 : float
        Depth of the elastic layer
        
    Returns:
    --------
    det_value : complex
        Determinant value
    """
    k_r = np.atleast_1d(k_r)
    n = len(k_r)
    
    # Calculate wave numbers and wave parameters
    i = 1j
    
    k_L = (omega/c_L)
    k_T = (omega/c_T)
    gamma_L = np.sqrt((k_L)**2 - k_r**2)
    gamma_T = np.sqrt((k_T)**2 - k_r**2)
    exp_GL = np.exp(-i * gamma_L * Z_1)
    exp_GT = np.exp(-i * gamma_T * Z_1)

    # Create N matrices of shape (4, 4) → full array of shape (N, 4, 4)
    M = np.empty((n, 4, 4), dtype=complex)

    M[:, 0, 0] = -(2*mu*gamma_L**2 + Lambda*k_L**2) * exp_GL
    M[:, 0, 1] = -2*mu*gamma_L**2 - Lambda*k_L**2
    M[:, 0, 2] = 2*i*mu*k_r*gamma_T * exp_GT
    M[:, 0, 3] = -2*i*mu*k_r*gamma_L

    M[:, 1, 0] = -2*i*k_r*gamma_L * exp_GL
    M[:, 1, 1] = 2*i*k_r*gamma_L
    M[:, 1, 2] = (gamma_T**2 - k_r**2) * exp_GT
    M[:, 1, 3] = gamma_T**2 - k_r**2

    M[:, 2, 0] = i*gamma_L
    M[:, 2, 1] = -i*gamma_L * exp_GL
    M[:, 2, 2] = k_r
    M[:, 2, 3] = k_r * exp_GT

    M[:, 3, 0] = -k_r
    M[:, 3, 1] = -k_r * exp_GL
    M[:, 3, 2] = -i*gamma_T
    M[:, 3, 3] = i*gamma_T * exp_GT

    # Now compute all determinants at once
    dets = np.array([np.linalg.det(M[j]) for j in range(n)])
    
    # magnitudes = np.abs(dets)
    # powers = np.floor(np.log10(magnitudes))
    # dets = dets * np.exp(-powers)

    return dets
        
# test

k = 1+1j
f(k)


# check for branch point
branch_point = k_p
theta = np.linspace(0, 2*np.pi, 100)
loop = branch_point + 1e-3 * np.exp(1j * theta)
values = f(loop)
plt.figure()
plt.plot(np.real(loop), np.imag(values))
plt.title("Looping around suspected branch point")


# %%% PLOTTING

# Create meshgrid and complex matrix Z
x = np.arange(real_min, real_max + dx, dx)
y = np.arange(imag_min, imag_max + dx, dx)
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y

# Evaluate f(z) on the grid
fz = np.vectorize(f)(Z)  # assuming f can be vectorized

# --- Optional plotting curves ---
del_ = 0
a = np.arange(1e-6, np.real(k_p), 1e-5)
b = np.real(k_p) * np.imag(k_p) / a
b0 = (np.real(k_p) * np.imag(k_p)) * (1 + 1e-10) / a

c = np.arange(1e-6, np.real(k_s), 1e-5)
d = np.real(k_s) * np.imag(k_s) / c
d0 = (-del_ + np.real(k_s) * np.imag(k_s)) * (1 - 1e-10) / c

# --- Plot ---
fig, ax = plt.subplots(figsize=(10, 6))

# Contours
cont1 = ax.contour(np.real(Z), np.imag(Z), np.real(fz), colors='b', linestyles='--')
cont2 = ax.contour(np.real(Z), np.imag(Z), np.imag(fz), colors='r', linestyles='-.')

# Characteristic curves
ax.plot(a, b, '--k', linewidth=0.5)
ax.plot(c, d, '--k', linewidth=0.5)

# Special points
ax.scatter(np.real(omega / c_L), np.imag(omega / c_L), color='magenta')
ax.scatter(np.real(omega / c_T), np.imag(omega / c_T), color='cyan')

# Optional surface plot (magnitude)
surf = ax.contourf(X, Y, np.abs(fz), levels=100, cmap='viridis')  # use contourf for 2D density
# To use a full 3D surface (like MATLAB's `surf`), use plot_surface with Axes3D

# Limits
ax.set_ylim([imag_min, imag_max])
ax.set_title('Contours of Re(f), Im(f), and |f(z)|')
ax.set_xlabel('Re(z)')
ax.set_ylabel('Im(z)')

# Colorbar
plt.colorbar(surf, ax=ax, label='|f(z)|')

plt.tight_layout()
plt.show()

# %% COUNT ROOTS AND PLOT

# initialise rootfinder
rootfinder = ComplexRootFinder(f)

# create rectangular domain
rectangle = RectangleContour(real_min, real_max, imag_min, imag_max)

# count roots in chosen domain
# rootfinder.count_root_containing_domains(rectangle)

# find all roots in chosen domain
all_roots = rootfinder.find_roots_domain(rectangle, max_roots_per_domain=3, max_depth=5, debug=True)

print(f"{all_roots}")
