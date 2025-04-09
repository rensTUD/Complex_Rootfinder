

# %% IMPORTS

import numpy as np
from abc import ABC, abstractmethod

# %% CLASSES

class ContourBase(ABC):
    @abstractmethod
    def Z(self) -> np.ndarray:
        """Return array of complex numbers defining the contour."""
        pass

    @abstractmethod
    def subdivide(self) -> list:
        """Return list of new ContourBase-derived instances (subdomains)."""
        pass

class RectangleContour(ContourBase):
    def __init__(
        self,
        real_min: float,
        real_max: float,
        imag_min: float,
        imag_max: float,
        n_points: int = 1000,
    ):
        """
        Creates a rectangular contour object, which has functions to return it's Z values in the complex domain
        and a method to subdivide itself into other rectangles

        Parameters
        ----------
        real_min : float
            left boundary of rectangle
        real_max : float
            rigth boundary of rectangle
        imag_min : float
            bottom boundary of rectangle
        imag_max : float
            top boundary of rectangle
        n_points : int, optional
            discretisation points, by default 1000
        """
        self.real_min = real_min
        self.real_max = real_max
        self.imag_min = imag_min
        self.imag_max = imag_max
        self.n_points = int(n_points)

    @property
    def Z(self) -> np.ndarray:
        """
        Generates the points Z in the complex domain

        Returns
        -------
        np.ndarray
            Z points in the complex domain
        """
        # set width (real valued) and heigth (imaginary valued)
        width = np.linspace(self.real_min, self.real_max, self.n_points)
        heigth = np.linspace(self.imag_min, self.imag_max, self.n_points)
        
        return np.concatenate([
            width + 1j * self.imag_min,                 # bottom
            self.real_max + 1j * heigth,                 # right
            width[::-1] + 1j * self.imag_max,           # top
            self.real_min + 1j * heigth[::-1],           # left
        ])

    def subdivide(
        self, 
        n_divide: int = 3
    ) -> list:
        """
        Subdivide rectangle into n_divide x n_divide smaller rectangles.

        Returns
        -------
        list of RectangleContour instances
        """
        new_width = self.width / n_divide
        new_height = self.height / n_divide

        subdomains = []
        for i in range(n_divide):
            for j in range(n_divide):
                re_min = self.real_min + i * new_width
                re_max = re_min + new_width
                im_min = self.imag_min + j * new_height
                im_max = im_min + new_height
                subdomains.append(RectangleContour(re_min, re_max, im_min, im_max, self.n_points))
        return subdomains
    
    @property
    def width(self) -> float:
        return self.real_max - self.real_min
    
    @property
    def height(self) -> float:
        return self.imag_max - self.imag_min
    
    @property
    def center(self) -> float:
        """
        center of the rectangle

        Returns
        -------
        float
            center
        """
        real = self.width / 2
        imag = self.height / 2
        return real + 1j * imag
    
    @property
    def radius(self) -> float:
        """
        returns the radius of a circle with the rectangle inscribed

        Returns
        -------
        float
            radius
        """
        return np.sqrt(self.width**2 + self.height**2) / 2

class CircleContour(ContourBase):
    def __init__(
        self,
        center: complex,
        radius: float,
        n_points: int = 1000,
    ):
        """
        Creates a circular contour object for root finding.

        Parameters
        ----------
        center : complex
            Center of the circle.
        radius : float
            Radius of the circle.
        n_points : int, optional
            Number of discretization points along the contour (default is 1000).
        """
        self.center = center
        self.radius = radius
        self.n_points = int(n_points)

    @property
    def Z(self) -> np.ndarray:
        """
        Generates the points Z in the complex domain

        Returns
        -------
        np.ndarray
            Z points in the complex domain
        """
        theta = np.linspace(0, 2 * np.pi, self.n_points, endpoint=True)
        return self.center + self.radius * np.exp(1j * theta)

    def subdivide(self, n_circle: int = 8) -> list:
        """
        Subdivide circle into n_circle overlapping smaller circles and 1 center circle.

        Returns
        -------
        list of CircleContour instances
        """
        subdomains = [CircleContour(self.center, self.radius / 2, self.n_points)]

        r0_n = 5 * self.radius / 12
        angles = 2 * np.pi * np.arange(n_circle) / n_circle
        centers = self.center + (3 / 4) * self.radius * np.exp(1j * angles)

        for c in centers:
            subdomains.append(CircleContour(c, r0_n, self.n_points))

        return subdomains

