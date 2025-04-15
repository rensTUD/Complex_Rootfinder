# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 17:27:51 2025

@author: RensvanLeijden
"""

# %% IMPORTS

import numpy as np
from typing import Callable
from matplotlib.path import Path
from shapely.geometry import LineString, Point
from typing import List, Tuple
from .contours import CompositeContour
import matplotlib.pyplot as plt

# %% FUNCTIONS

def plot_contour(Z, ax=None, color='black', label=None, lw=1.5):
    """Utility to plot a complex-valued closed path Z."""
    if ax is None:
        ax = plt.gca()
    xy = np.array([[z.real, z.imag] for z in Z])
    ax.plot(xy[:, 0], xy[:, 1], color=color, lw=lw, label=label)
    return ax

def finite_difference_first_derivative_4th_order(
        f: Callable[[np.ndarray], np.ndarray], 
        Z
    ):
        """
        Compute derivative using 4th order accurate finite difference.
        
        Parameters
        ----------
            f : callable 
                the function to differentiate
            Z : array-like
                points at which to evaluate the derivative
        
        Returns:
        -------
            array-like
                The approximated derivative values
        """
        h = 1e-5
        return (-f(Z + 2*h) + 8*f(Z + h) - 8*f(Z - h) + f(Z - 2*h)) / (12*h)
    
def partition_contour_along_cut_general(
    contour_Z: np.ndarray,
    cut_points: List[complex],
    epsilon: float = 1e-2
) -> List[List[complex]]:
    """
    Partition a closed contour along a branch cut into multiple closed paths.

    Parameters
    ----------
    contour_Z : List[complex]
        The original closed contour path.
    cut_points : List[complex]
        Points defining the branch cut as a polyline.
    epsilon : float
        Offset distance to avoid evaluating on the cut directly.

    Returns
    -------
    List[List[complex]]
        List of new closed domain paths as lists of complex points.
    """
    insertions = find_intersections(contour_Z, cut_points)
    if len(insertions) < 2 or len(insertions) % 2 != 0:
        raise ValueError("Expected even number of intersection points (≥2).")

    Z_aug = insert_points_into_path(contour_Z, insertions)
    inserted_points = [pt for _, pt in sorted(insertions)]
    inserted_indices = [Z_aug.index(pt) for pt in inserted_points]

    arcs = split_contour_along_intersections(Z_aug, inserted_indices)

    subdomains = []
    for i in range(0, len(arcs), 2):
        arc = arcs[i]
        pt1 = arc[0]
        pt2 = arc[-1]
        cut_seg = extract_cut_segment(cut_points, pt1, pt2)

        # Offset the cut segment to avoid the true discontinuity
        cut_seg_left = [z + epsilon * 1j for z in cut_seg]
        loop = arc + cut_seg_left[::-1]  # closing counter-clockwise
        if not np.isclose(loop[0], loop[-1]):
            loop.append(loop[0])
        subdomains.append(loop)

    return subdomains

def is_point_inside_contour(Z: np.ndarray, point: complex) -> bool:
    """
    Determine whether a complex point lies inside a closed complex-valued contour.

    Parameters
    ----------
    Z : List[complex]
        Closed path representing the contour.
    point : complex
        Point to test.

    Returns
    -------
    bool
        True if the point lies inside the contour.
    """
    path = Path([(z.real, z.imag) for z in Z])
    return path.contains_point((point.real, point.imag))

# def wrap_contour_around_branch_point(contour, cut, epsilon=1e-4) -> "CompositeContour":
#     """
#     Create a composite contour that wraps around a branch point and its associated cut.

#     Parameters
#     ----------
#     contour : ContourBase or CompositeContour
#         The original enclosing domain.
#     cut : BranchCut
#         A BranchCut object that includes a defined branch point.
#     epsilon : float
#         Offset distance from the cut to avoid evaluating on the discontinuity.

#     Returns
#     -------
#     CompositeContour
#         A new composite contour that encircles the cut and the branch point.
#     """
#     if not cut.branch_point:
#         raise ValueError("wrap_contour_around_branch_point requires a branch point.")

#     boundary = contour.Z
#     left_offset = offset_branchcut(cut.points, epsilon, direction='left')
#     right_offset = offset_branchcut(cut.points, epsilon, direction='right')

#     segments = [boundary, left_offset[::-1], right_offset]
#     return CompositeContour(segments)

def wrap_contour_around_branch_point(contour, cut, epsilon=1e-4) -> "CompositeContour":
    """
    Insert an offset loop around the branch cut and branch point, direction-aware.

    Parameters
    ----------
    contour : ContourBase or CompositeContour
        The enclosing domain.
    cut : BranchCut
        Branch cut object with a defined branch point and .points array.
    epsilon : float
        Distance of the offset away from the cut path.

    Returns
    -------
    CompositeContour
        A new contour with the loop inserted around the cut.
    """
    if not cut.branch_point:
        raise ValueError("wrap_contour_around_branch_point requires a branch point.")

    Z: np.ndarray = contour.Z
    insertions = find_intersections(Z, cut.points)

    if len(insertions) != 1:
        raise ValueError("Expected exactly one intersection between the contour and the branch cut.")

    intersection_point = insertions[0][1]
    insert_index = insertions[0][0]

    before = Z[:insert_index]
    after = Z[insert_index:]

    # Build offset paths in correct direction
    first_segment = build_offset_segment_along_cut(
        np.array(cut.points), pt_start=before[-1], pt_end=cut.branch_point,
        epsilon=epsilon, side="first"
    )

    second_segment = build_offset_segment_along_cut(
        np.array(cut.points), pt_start=cut.branch_point, pt_end=after[0],
        epsilon=epsilon, side="second"
    )

    # Create arc bridging the two segments around the branch point
    arc = create_arc_between_points(
        start=first_segment[-1],
        end=second_segment[0],
        intersection=intersection_point,
        n_points=20
    )

    # Join into new composite contour
    segments = [
        before,
        first_segment,
        arc,
        second_segment,
        after
    ]

    return CompositeContour(segments)

def find_intersections(contour_Z: np.ndarray, cut_points: np.ndarray) -> List[Tuple[int, complex]]:
    contour_line = LineString([(z.real, z.imag) for z in contour_Z])
    cut_line = LineString([(z.real, z.imag) for z in cut_points])
    intersections = contour_line.intersection(cut_line)

    if intersections.is_empty:
        return []

    if isinstance(intersections, Point):
        intersections = [intersections]
    elif hasattr(intersections, 'geoms'):
        intersections = list(intersections.geoms)

    result = []
    for pt in intersections:
        pt_complex = complex(pt.x, pt.y)
        for i in range(len(contour_Z) - 1):
            seg = LineString([(contour_Z[i].real, contour_Z[i].imag),
                              (contour_Z[i + 1].real, contour_Z[i + 1].imag)])
            if seg.distance(pt) < 1e-8:
                result.append((i + 1, pt_complex))
                break
    return result


def insert_points_into_path(Z: np.ndarray, insertions: List[Tuple[int, complex]]) -> np.ndarray:
    Z_out = Z.copy()
    for idx, pt in sorted(insertions, key=lambda x: x[0]):
        Z_out = np.insert(Z_out, idx, pt)
    return Z_out


def split_contour_along_intersections(Z: np.ndarray, intersection_indices: List[int]) -> List[np.ndarray]:
    arcs = []
    n = len(intersection_indices)
    for i in range(n):
        start = intersection_indices[i]
        end = intersection_indices[(i + 1) % n]
        if start < end:
            arc = Z[start:end + 1]
        else:
            arc = np.concatenate((Z[start:], Z[:end + 1]))
        arcs.append(arc)
    return arcs


def extract_cut_segment(cut: np.ndarray, pt_start: complex, pt_end: complex) -> np.ndarray:
    """
    Extract a directional segment from an offset cut path, ensuring it flows
    from pt_start to pt_end along the cut geometry.

    Parameters
    ----------
    cut : np.ndarray
        Full offset cut path.
    pt_start : complex
        Start point (usually insertion point on domain).
    pt_end : complex
        End point (usually branch point).

    Returns
    -------
    np.ndarray
        Segment along the cut from pt_start to pt_end.
    """
    cut_line = LineString([(z.real, z.imag) for z in cut])
    start_dist = cut_line.project(Point(pt_start.real, pt_start.imag))
    end_dist = cut_line.project(Point(pt_end.real, pt_end.imag))

    n_samples = 100
    sample_distances = (
        np.linspace(start_dist, end_dist, n_samples)
        if start_dist <= end_dist
        else np.linspace(start_dist, end_dist, n_samples)[::-1]
    )

    segment = [complex(*cut_line.interpolate(d).coords[0]) for d in sample_distances]
    return np.array(segment, dtype=complex)[1:-1]


def offset_branchcut(points: np.ndarray, epsilon: float, direction: int = +1) -> np.ndarray:
    """
    Offset a polyline path representing the branch cut.

    Parameters
    ----------
    points : np.ndarray
        Array of complex points defining the cut.
    epsilon : float
        Offset distance.
    direction : int
        +1 for one side, -1 for the opposite.

    Returns
    -------
    np.ndarray
        Offset complex-valued path.
    """
    offset_path = []
    n = len(points)
    for i in range(n - 1):
        p0, p1 = points[i], points[i + 1]
        delta = p1 - p0
        perp = delta * 1j
        perp /= abs(perp)
        perp *= direction
        offset_path.extend([p0 + epsilon * perp, p1 + epsilon * perp])
    return np.array(offset_path, dtype=complex)

def build_offset_segment_along_cut(
    cut: np.ndarray,
    pt_start: complex,
    pt_end: complex,
    epsilon: float,
    side: str  # "first" or "second"
) -> np.ndarray:
    """
    Build an offset segment along a cut path, from pt_start to pt_end, and offset
    to the appropriate side ("first" or "second").

    Parameters
    ----------
    cut : np.ndarray
        The full cut as complex-valued points.
    pt_start : complex
        Start point along the cut (projected).
    pt_end : complex
        End point along the cut (projected).
    epsilon : float
        Offset distance perpendicular to the cut direction.
    side : str
        Either "first" or "second". Determines side of the offset.

    Returns
    -------
    np.ndarray
        Offset path segment from pt_start to pt_end.
    """
    cut_line = LineString([(z.real, z.imag) for z in cut])
    start_dist = cut_line.project(Point(pt_start.real, pt_start.imag))
    end_dist = cut_line.project(Point(pt_end.real, pt_end.imag))

    # Ensure we move from pt_start to pt_end along the cut
    sample_distances = np.linspace(start_dist, end_dist, 100)

    segment = [cut_line.interpolate(d) for d in sample_distances]
    segment_coords = np.array([[p.x, p.y] for p in segment])
    segment_complex = segment_coords[:, 0] + 1j * segment_coords[:, 1]

    # Compute tangent vectors and normals
    tangents = np.diff(segment_complex)
    tangents = np.append(tangents, tangents[-1])  # repeat last for shape
    normals = tangents * 1j
    normals /= np.abs(normals)

    # direction = +1 if side == "first" else -1
    offset_segment = segment_complex + epsilon * normals
    return offset_segment.astype(complex)

def create_arc_between_points(
    start: complex,
    end: complex,
    intersection: complex,
    n_points: int = 20
) -> np.ndarray:
    """
    Create an arc between two points, going around the branch point in a direction
    determined by the location of the intersection relative to the arc center.

    Parameters
    ----------
    start : complex
        Start point of the arc.
    end : complex
        End point of the arc.
    intersection : complex
        Intersection point of the branch cut with the domain boundary.
    n_points : int
        Number of points to sample along the arc.

    Returns
    -------
    np.ndarray
        Complex-valued points forming the arc.
    """
    center = 0.5 * (start + end)
    radius = 0.5 * abs(end - start)

    v_start = start - center
    v_end = end - center

    angle_start = np.angle(v_start)
    angle_end = np.angle(v_end)

    # Determine desired rotation direction based on geometry
    if start.imag < end.imag:
        # going bottom to top
        clockwise = intersection.real > center.real
    else:
        # going top to bottom
        clockwise = intersection.real < center.real

    if clockwise:
        if angle_start < angle_end:
            angle_start += 2 * np.pi
        angles = np.linspace(angle_start, angle_end, n_points)
    else:
        if angle_end < angle_start:
            angle_end += 2 * np.pi
        angles = np.linspace(angle_start, angle_end, n_points)

    arc = center + radius * np.exp(1j * angles)
    return arc.astype(complex)



