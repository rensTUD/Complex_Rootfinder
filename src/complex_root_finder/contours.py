

# %% IMPORTS

import numpy as np
from abc import ABC, abstractmethod
from typing import List, Union
from shapely.geometry import Polygon, box as shapely_box, Point, LineString, Polygon
from shapely.geometry.polygon import orient

# %% CLASSES

class ContourBase(ABC):
    # @abstractmethod
    def Z(self) -> np.ndarray:
        """Return array of complex numbers defining the contour."""
        pass

    @abstractmethod
    def subdivide(self) -> list:
        """Return list of new ContourBase-derived instances (subdomains)."""
        pass
    
    def resample_from_boundary(boundary: Polygon, dz: float = None, n_points: int = None) -> np.ndarray:
        """
        Resample the exterior of a polygon boundary into complex points.
    
        Parameters
        ----------
        boundary : shapely.geometry.Polygon
            The polygon to resample.
        dz : float, optional
            Desired arc-length spacing between points.
        n_points : int, optional
            Number of points to sample.
    
        Returns
        -------
        np.ndarray
            Complex-valued array of points tracing the boundary.
        """
        line = LineString(boundary.exterior.coords)
        length = line.length
    
        if dz is not None:
            n_points = max(int(np.ceil(length / dz)), 3)
        elif n_points is not None:
            n_points = max(int(n_points), 3)
        else:
            raise ValueError("You must specify either dz or n_points")
    
        distances = np.linspace(0, length, n_points)
        points = [line.interpolate(d) for d in distances]
        return np.array([complex(p.x, p.y) for p in points])

class RectangleContour(ContourBase):
    def __init__(
        self,
        real_min: float,
        real_max: float,
        imag_min: float,
        imag_max: float,
        n_points: int = 1000,
        dz: float | int = 1e-3
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
        self.boundary = shapely_box(real_min, imag_min, real_max, imag_max)
        self.n_points = int(n_points)
        self.dz = dz
        
    def domain(self):
        """
        Prints the current domain

        Returns
        -------
        None.

        """
        print(
            f"real_min = {self.real_min}\n"
            f"real_max = {self.real_max}\n"
            f"imag_min = {self.imag_min}\n"
            f"imag_max = {self.imag_max}\n"
        )
        
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
        # width = np.linspace(self.real_min, self.real_max, self.n_points)
        # heigth = np.linspace(self.imag_min, self.imag_max, self.n_points)
        width = np.arange(self.real_min, self.real_max, self.dz)
        heigth = np.arange(self.imag_min, self.imag_max, self.dz)
        
        return np.concatenate([
            width + 1j * self.imag_min,                 # bottom
            self.real_max + 1j * heigth,                 # right
            width[::-1] + 1j * self.imag_max,           # top
            self.real_min + 1j * heigth[::-1],           # left
        ])

    def subdivide(
        self, 
        n_divide: int = 3,
        overlap_percent: float = 0.3
    ) -> list:
        """
        Subdivide rectangle into n_divide x n_divide overlapping rectangles.
    
        Parameters
        ----------
        n_divide : int
            Number of divisions along each axis.
        overlap_percent : float, optional
            Fraction (e.g. 0.1 for 10%) to enlarge each subrectangle in all directions.
    
        Returns
        -------
        list of RectangleContour instances
        """
        new_width = self.width / n_divide
        new_height = self.height / n_divide
    
        # How much each side is extended (in real and imag axes)
        overlap_re = new_width * overlap_percent / 2
        overlap_im = new_height * overlap_percent / 2
    
        subdomains = []
        for i in range(n_divide):
            for j in range(n_divide):
                # Start from unmodified bounds
                re_min = self.real_min + i * new_width
                re_max = re_min + new_width
                im_min = self.imag_min + j * new_height
                im_max = im_min + new_height
    
                # Expand each rectangle, but clip to the parent bounds
                re_min = max(self.real_min, re_min - overlap_re)
                re_max = min(self.real_max, re_max + overlap_re)
                im_min = max(self.imag_min, im_min - overlap_im)
                im_max = min(self.imag_max, im_max + overlap_im)
    
                subdomains.append(RectangleContour(re_min, re_max, im_min, im_max, self.n_points, self.dz))
        
        return subdomains
    
    @property
    def bounds(self) -> list[tuple[float, float]]:
        return [
            (self.real_min, self.real_max),
            (self.imag_min, self.imag_max)
        ]
    
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
        real_center = (self.real_min + self.real_max) / 2
        imag_center = (self.imag_min + self.imag_max) / 2
        return real_center + 1j * imag_center
    
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
        dz: float = 1e-3
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
        self.boundary = Point(self.center.real, self.center.imag).buffer(self.radius, resolution=self.n_points)
        self.n_points = int(n_points)
        self.dz = dz

    def domain(self):
        """
        Prints the current domain

        Returns
        -------
        None.

        """
        print(
            f"center = {self.center}\n"
            f"radius = {self.radius}\n"
        )

    @property
    def Z(self) -> np.ndarray:
        """
        Generates the points Z in the complex domain

        Returns
        -------
        np.ndarray
            Z points in the complex domain
        """
        # theta = np.linspace(0, 2 * np.pi, self.n_points, endpoint=True)
        theta = np.arange(0, 2 * np.pi, self.dz, endpoint=True)
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
            subdomains.append(CircleContour(c, r0_n, self.n_points, self.dz))

        return subdomains
    
    @property
    def bounds(self) -> list[tuple[float, float]]:
        return [
            (self.center.real - self.radius, self.center.real + self.radius),
            (self.center.imag - self.radius, self.center.imag + self.radius)
        ]

class CompositeContour(ContourBase):
    def __init__(
        self, 
        segments: List[np.ndarray],
        n_points: int = 1000,
        dz: int = 1e-3,
        boundary: Polygon = None
    ):
        """
        Represents a composite contour made of multiple segments.
        Each segment is an np.ndarray of complex numbers.
        The full path is flattened and closed automatically.

        Parameters
        ----------
        segments : List[np.ndarray]
            A list of path segments, each segment being an np.ndarray of complex points.
        """
        if not segments or not all(isinstance(seg, np.ndarray) and len(seg) >= 2 for seg in segments):
            raise ValueError("Each segment must be a numpy array with at least two points.")
            
        
        self.segments = segments
        self.Z = self._build_Z()
        self._compute_bounds()
        if boundary:
            self.boundary = boundary
        else:
            self.boundary = Polygon([(z.real, z.imag) for z in self.Z])
            if not self.boundary.is_valid:
                self.boundary = self.boundary.buffer(0)
        
        self.n_points = n_points
        self.dz = dz

    def _build_Z(self) -> np.ndarray:
        """
        Flatten the list of segments into a single closed path.

        Returns
        -------
        np.ndarray
            The concatenated, closed path.
        """
        flat = np.concatenate(self.segments)
        if not np.isclose(flat[0], flat[-1]):
            flat = np.append(flat, flat[0])  # ensure closed path
        return flat

    def _compute_bounds(self):
        """
        Compute bounding box (xmin, xmax), (ymin, ymax).
        """
        self._xmin = np.min(self.Z.real)
        self._xmax = np.max(self.Z.real)
        self._ymin = np.min(self.Z.imag)
        self._ymax = np.max(self.Z.imag)

    @property
    def bounds(self) -> List[tuple]:
        """
        Returns
        -------
        List[Tuple[float, float]]
            Bounding box as [(xmin, xmax), (ymin, ymax)]
        """
        return [(self._xmin, self._xmax), (self._ymin, self._ymax)]
    
    @property
    def center(self) -> complex:
        """
        Returns the center of the bounding box of the contour.
    
        Returns
        -------
        complex
            Center of the contour.
        """
        x_center = (self._xmin + self._xmax) / 2
        y_center = (self._ymin + self._ymax) / 2
        return complex(x_center, y_center)

    @property
    def radius(self) -> float:
        """
        Returns the radius of the smallest circle that contains the bounding box.
    
        Returns
        -------
        float
            Radius of the contour's bounding circle.
        """
        width = self._xmax - self._xmin
        height = self._ymax - self._ymin
        return 0.5 * np.sqrt(width**2 + height**2)


    def __repr__(self):
        return f"CompositeContour(n_segments={len(self.segments)}, n_points={len(self.Z)})"

    def subdivide(self, n_divide: int = 3, overlap_percent: float = 0.3) -> List["CompositeContour"]:
        """
        Subdivide the composite contour using geometry-aware polygon clipping,
        and resample each resulting polygon's boundary using adaptive spacing.
    
        Parameters
        ----------
        n_divide : int
            Number of subdivisions along each axis.
        overlap_percent : float
            Fraction (e.g. 0.1 for 10%) to enlarge each subrectangle in all directions.
    
        Returns
        -------
        List[CompositeContour]
            Subdomains as uniformly-sampled CompositeContours.
        """
        dz = self.dz
    
        if not self.boundary.is_valid:
            poly = self.boundary.buffer(0)
        else:
            poly = self.boundary
    
        subdomains = []
        (xmin, xmax), (ymin, ymax) = self.bounds
    
        x_vals = np.linspace(xmin, xmax, n_divide + 1)
        y_vals = np.linspace(ymin, ymax, n_divide + 1)
        new_width = (xmax - xmin) / n_divide
        new_height = (ymax - ymin) / n_divide
        overlap_x = new_width * overlap_percent / 2
        overlap_y = new_height * overlap_percent / 2
    
        for i in range(n_divide):
            for j in range(n_divide):
                x0 = max(xmin, x_vals[i] - overlap_x)
                x1 = min(xmax, x_vals[i + 1] + overlap_x)
                y0 = max(ymin, y_vals[j] - overlap_y)
                y1 = min(ymax, y_vals[j + 1] + overlap_y)
                rect = shapely_box(x0, y0, x1, y1)
    
                try:
                    inter = poly.intersection(rect)
                except Exception as e:
                    print(f"[Warning] Skipping cell ({i},{j}): intersection failed: {e}")
                    continue
    
                if inter.is_empty:
                    continue
    
                parts = [inter] if inter.geom_type == "Polygon" else (
                    list(inter.geoms) if inter.geom_type == "MultiPolygon" else []
                )
    
                for part in parts:
                    part = orient(part, sign=1.0)
                    resampled_Z = self.adaptive_resample_polygon_boundary(part, dz=dz)
                    subdomains.append(CompositeContour([resampled_Z], boundary=part))
    
        return subdomains

    # def subdivide(self, n_divide: int = 3) -> List["CompositeContour"]:
    #         """
    #         Subdivide the composite contour using geometry-aware polygon clipping
    #         and resample each resulting polygon's boundary using uniform spacing.
    
    #         Parameters
    #         ----------
    #         n_divide : int
    #             Number of subdivisions along each axis.
    #         dz : float
    #             Desired spacing for resampling subdomain boundaries.
    
    #         Returns
    #         -------
    #         List[CompositeContour]
    #             Subdomains as uniformly-sampled CompositeContours.
    #         """
    #         dz = self.dz        
    
    #         poly = Polygon([(z.real, z.imag) for z in self.Z])
    #         if not poly.is_valid:
    #             poly = poly.buffer(0)
    
    #         subdomains = []
    #         (xmin, xmax), (ymin, ymax) = self.bounds
    #         x_vals = np.linspace(xmin, xmax, n_divide + 1)
    #         y_vals = np.linspace(ymin, ymax, n_divide + 1)
    
    #         for i in range(n_divide):
    #             for j in range(n_divide):
    #                 x0, x1 = x_vals[i], x_vals[i + 1]
    #                 y0, y1 = y_vals[j], y_vals[j + 1]
    #                 rect = shapely_box(x0, y0, x1, y1)
    
    #                 try:
    #                     inter = poly.intersection(rect)
    #                 except Exception as e:
    #                     print(f"[Warning] Skipping cell ({i},{j}): intersection failed: {e}")
    #                     continue
    
    #                 if inter.is_empty:
    #                     continue
    
    #                 parts = [inter] if inter.geom_type == "Polygon" else (
    #                     list(inter.geoms) if inter.geom_type == "MultiPolygon" else []
    #                 )
    
    #                 for part in parts:
    #                     part = orient(part, sign=1.0)
    #                     # resampled_Z = self.resample_polygon_boundary(part, dz=dz)
    #                     resampled_Z = self.adaptive_resample_polygon_boundary(part, dz=dz)
    #                     subdomains.append(CompositeContour([resampled_Z], boundary=part))
    
    #         return subdomains

    def resample_polygon_boundary(self, polygon, dz: float) -> np.ndarray:
        """
        Uniformly sample a counter-clockwise boundary along the polygon exterior.
    
        Parameters
        ----------
        polygon : shapely.geometry.Polygon
            The polygon to resample.
        dz : float
            Desired step size.
    
        Returns
        -------
        np.ndarray
            Complex-valued, uniformly spaced counter-clockwise loop.
        """
        polygon = orient(polygon, sign=1.0)
        line = LineString(polygon.exterior.coords)
        length = line.length
        n_points = max(int(np.ceil(length / dz)), 3)
        distances = np.linspace(0, length, n_points)
        points = [line.interpolate(d) for d in distances]
        return np.array([complex(p.x, p.y) for p in points])

    def adaptive_resample_polygon_boundary(self, polygon, dz: float) -> np.ndarray:
        """
        Reuse original boundary points when spacing is fine enough.
        Interpolate intermediate points where gaps exceed `dz`.
    
        Parameters
        ----------
        polygon : shapely.geometry.Polygon
            The polygon to resample adaptively.
        dz : float
            Desired maximum spacing between consecutive boundary points.
    
        Returns
        -------
        np.ndarray
            Complex-valued array of boundary points with adaptive resolution.
        """
        polygon = orient(polygon, sign=1.0)
        coords = list(polygon.exterior.coords)
        result = []
    
        for i in range(len(coords) - 1):
            p0 = coords[i]
            p1 = coords[i + 1]
            z0 = complex(*p0)
            z1 = complex(*p1)
    
            segment_length = abs(z1 - z0)
            n_segment = max(int(np.ceil(segment_length / dz)), 1)
    
            segment = np.linspace(z0, z1, n_segment + 1)[:-1]
            result.extend(segment)
    
        result.append(complex(*coords[-1]))  # close the loop
    
        return np.array(result)


# %% BranchCut class

class BranchCut:
    def __init__(self, points: np.ndarray, branch_point: complex | None = None):
        """
        Parameters
        ----------
        points : List[complex]
            Piecewise linear path representing the branch cut.
        branch_point : complex, optional
            The singular point at the origin of the cut.
        """
        if len(points) < 2:
            raise ValueError("BranchCut requires at least two points.")
        self.points = points
        self.branch_point = branch_point
        self._compute_bounds()
        self._line = LineString([(p.real, p.imag) for p in points])

    def _compute_bounds(self):
        real_parts = [p.real for p in self.points]
        imag_parts = [p.imag for p in self.points]
        self._xmin = min(real_parts)
        self._xmax = max(real_parts)
        self._ymin = min(imag_parts)
        self._ymax = max(imag_parts)

    def intersects(self, contour) -> bool:
        """
        Efficiently check if this branch cut intersects the provided contour.

        Parameters
        ----------
        contour : Any object with `.Z` and `.bounds`.

        Returns
        -------
        bool
            True if the cut intersects the contour.
        """
        (cxmin, cxmax), (cymin, cymax) = contour.bounds

        # Quick bounding box rejection
        if (
            cxmax < self._xmin or cxmin > self._xmax or
            cymax < self._ymin or cymin > self._ymax
        ):
            return False
        
        # create boundary
        # contour_line = LineString([(z.real, z.imag) for z in contour.Z])
        # return self._line.intersects(contour_line)
        
        # use boundary of contour
        return self._line.intersects(contour.boundary)

    def __repr__(self):
        return f"BranchCut(n_points={len(self.points)}, branch_point={self.branch_point})"

