# -*- coding: utf-8 -*-
"""
Created on Mon Apr 9 20:02:28 2025

@author: RensvanLeijden
"""

# %% IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from complex_root_finder import ComplexRootFinder, RectangleContour, CircleContour, BranchCut


# %% TEST 1 - single soil layer?

# %%% INPUTS
c_p = 1.719084289205639e+03 + 3.322560356224935e+01j
c_s = 3.704985415873700e+02 + 1.875674989161015e+01j
Lambda = 5.114035573342751e+09 + 1.649233839281022e+08j
mu = 2.612383107677712e+08 + 2.651871379876653e+07j

rho_f = 1000
c_f = 1500


omega = 1000*2*np.pi
H = 10
H0 = H
D0 = H
k_p = omega/c_p
k_s = omega/c_s
k_f = omega/c_f
Z_1 = H

def f(k_x): 
    return 32 * ((((-(0.1e1 / 0.16e2) - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 16 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 16 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 16) * mu + (-(0.1e1 / 0.32e2) - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 32 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 32 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 32) * Lambda) * (omega ** 4) - (-(0.1e1 / 0.4e1) - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 4 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 4 - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 4 + np.exp(-1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2)))) * (mu * (c_p ** 2) + 2 * (c_s ** 2) * (mu + Lambda / 2)) * (k_x ** 2) * (omega ** 2) / 4 + (-(0.1e1 / 0.4e1) - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 4 - np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 4 - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 4 + np.exp(-1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2)))) * (c_p ** 2) * (c_s ** 2) * mu * (k_x ** 4)) * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + (np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) + np.exp(-2*1j * np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) - np.exp(-2*1j * H * (np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) - 1) * ((omega ** 4) * ((0.3e1 / 0.4e1) * mu + Lambda / 8) - (0.3e1 / 0.4e1) * (mu * (c_p ** 2) + (0.4e1 / 0.3e1) * (c_s ** 2) * (mu + Lambda / 4)) * (k_x ** 2) * (omega ** 2) + (c_p ** 2) * (c_s ** 2) * (k_x ** 4) * mu) * (k_x ** 2) / 4) * mu / (c_p ** 2) / (c_s ** 2)
# def f2(k_x): 
#     return 32 * ((((-(0.1e1 / 0.16e2) - np.exp(-2*1j * H * (-np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 16 - np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 16 - np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 16) * mu + (-(0.1e1 / 0.32e2) - np.exp(-2*1j * H * (-np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 32 - np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 32 - np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 32) * Lambda) * (omega ** 4) - (-(0.1e1 / 0.4e1) - np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 4 - np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 4 - np.exp(-2*1j * H * (-np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 4 + np.exp(-1j * H * (-np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2)))) * (mu * (c_p ** 2) + 2 * (c_s ** 2) * (mu + Lambda / 2)) * (k_x ** 2) * (omega ** 2) / 4 + (-(0.1e1 / 0.4e1) - np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) / 4 - np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) / 4 - np.exp(-2*1j * H * (-np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) / 4 + np.exp(-1j * H * (-np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2)))) * (c_p ** 2) * (c_s ** 2) * mu * (k_x ** 4)) * -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * -np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + (np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) * H) + np.exp(-2*1j * -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2) * H) - np.exp(-2*1j * H * (-np.sqrt((-k_x ** 2 * c_p ** 2 + omega ** 2) / c_p ** 2) + -np.sqrt((-k_x ** 2 * c_s ** 2 + omega ** 2) / c_s ** 2))) - 1) * ((omega ** 4) * ((0.3e1 / 0.4e1) * mu + Lambda / 8) - (0.3e1 / 0.4e1) * (mu * (c_p ** 2) + (0.4e1 / 0.3e1) * (c_s ** 2) * (mu + Lambda / 4)) * (k_x ** 2) * (omega ** 2) + (c_p ** 2) * (c_s ** 2) * (k_x ** 4) * mu) * (k_x ** 2) / 4) * mu / (c_p ** 2) / (c_s ** 2)
# def f(k_x): 
#     return f1(k_x)*f2(k_x)
def f(k_r):
    k_zT = np.sqrt((k_s)**2 - k_r**2)
    k_zL = np.sqrt((k_p)**2 - k_r**2)
    cg = Lambda
    return complex(0, -2) * mu * (((k_zT + k_r) * (mu + cg / 2) * (k_r - k_zT) * k_zL ** 2 - 2 * k_r ** 2 * k_zL * k_zT * mu + k_r ** 4 * cg / 2 - k_r ** 2 * k_zT ** 2 * cg / 2) * (k_r ** 2 + k_zT * k_zL) * np.exp(complex(0, 2) * Z_1 * (k_zL + k_zT)) + 4 * k_r ** 2 * k_zL * ((-2 * mu - cg) * k_zL ** 2 + (mu - cg) * k_r ** 2 - k_zT ** 2 * mu) * k_zT * np.exp(complex(0, 1) * Z_1 * (k_zL + k_zT)) - ((k_zT + k_r) * (mu + cg / 2) * (k_r - k_zT) * k_zL ** 2 + 2 * k_r ** 2 * k_zL * k_zT * mu + k_r ** 4 * cg / 2 - k_r ** 2 * k_zT ** 2 * cg / 2) * (k_r ** 2 - k_zT * k_zL) * np.exp(complex(0, 2) * k_zL * Z_1) - ((k_zT + k_r) * (mu + cg / 2) * (k_r - k_zT) * k_zL ** 2 + 2 * k_r ** 2 * k_zL * k_zT * mu + k_r ** 4 * cg / 2 - k_r ** 2 * k_zT ** 2 * cg / 2) * (k_r ** 2 - k_zT * k_zL) * np.exp(complex(0, 2) * k_zT * Z_1) + ((k_zT + k_r) * (mu + cg / 2) * (k_r - k_zT) * k_zL ** 2 - 2 * k_r ** 2 * k_zL * k_zT * mu + k_r ** 4 * cg / 2 - k_r ** 2 * k_zT ** 2 * cg / 2) * (k_r ** 2 + k_zT * k_zL))



# check for branch point
branch_point = k_p
theta = np.linspace(0, 2*np.pi, 100)
loop = branch_point + 1e-3 * np.exp(1j * theta)
values = f(loop)
plt.figure()
plt.plot(np.real(loop), np.imag(values))
plt.title("Looping around suspected branch point")

# --- Parameters ---
dx = 5e-3
real_min = -np.max(np.real(k_s)) / 0.5
real_min = -1
real_max = np.max(np.real(k_s)) / 0.8
imag_min = -3
imag_max = 50e-6  

real_min = -np.max(np.real(k_s)) / 0.5
real_min = -1
real_max = np.max(np.real(k_s)) / 0.8
imag_min = -3
imag_max = -5e-2 
 



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
a = np.arange(1e-6, np.real(k_p), 1e-4)
a = np.arange(np.real(k_p), 1.02*real_max, 1e-4)
b = np.real(k_p) * np.imag(k_p) / a
b0 = (np.real(k_p) * np.imag(k_p)) * (1 + 1e-10) / a

c = np.arange(1e-6, np.real(k_s), 1e-4)
c = np.arange(np.real(k_s), 1.02*real_max, 1e-4)
d = np.real(k_s) * np.imag(k_s) / c
d0 = (-del_ + np.real(k_s) * np.imag(k_s)) * (1 - 1e-10) / c

# %%%%
# --- Plot ---
fig, ax = plt.subplots(figsize=(10, 6))

# Contours
cont1 = ax.contour(np.real(Z), np.imag(Z), np.real(fz), levels=[0], colors='b', linestyles='--')
cont2 = ax.contour(np.real(Z), np.imag(Z), np.imag(fz), levels=[0], colors='r', linestyles='-.')

# Characteristic curves
ax.plot(a, b, '--k', linewidth=0.5)
ax.plot(c, d, '--k', linewidth=0.5)

# Special points
ax.scatter(np.real(omega / c_p), np.imag(omega / c_p), color='magenta')
ax.scatter(np.real(omega / c_s), np.imag(omega / c_s), color='cyan')

# Optional surface plot (magnitude)
# surf = ax.contourf(X, Y, np.abs(fz), levels=100, cmap='viridis')  # use contourf for 2D density


# Limits
ax.set_ylim([imag_min, imag_max])
ax.set_title('Contours of Re(f), Im(f), and |f(z)|')
ax.set_xlabel('Re(z)')
ax.set_ylabel('Im(z)')



# Colorbar
# plt.colorbar(surf, ax=ax, label='|f(z)|')

plt.tight_layout()
plt.show()




# %%% root counting

# define branch cuts
branch_cut_kp = BranchCut(a+b*1j, k_p)
branch_cut_ks = BranchCut(c+d*1j, k_s)
branch_cuts = [branch_cut_kp, branch_cut_ks]

# initialise rootfinder
rootfinder = ComplexRootFinder(f, branch_cuts = branch_cuts)

# create rectangular domain
rectangle = RectangleContour(real_min, real_max, imag_min, imag_max)

# count roots in chosen domain
# rootfinder.count_root_containing_domains(rectangle)


# find all roots in chosen domain
# all_roots = rootfinder.find_roots_domain(rectangle, max_roots_per_domain=4, max_depth=5, debug=True)

# test version of above
all_roots = rootfinder.find_roots_domain_test(rectangle, max_roots_per_domain=1, max_depth=15, debug=True)



print(f"{all_roots}")