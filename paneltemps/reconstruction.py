"""Reconstructed CAD-surface quadrature, explicitly distinct from a CFD mesh."""

from collections import defaultdict
from dataclasses import dataclass
import math

import numpy as np
from scipy.spatial import Delaunay, QhullError
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepTools import BRepTools
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCP.gp import gp_Pnt

from .cad import quadrature


def unique_samples(nodes, tolerance=1e-6):
    groups = defaultdict(list)
    for node in nodes.values():
        groups[node.xyz].append(node)
    points = sorted(groups)
    temperatures, ids, duplicates = [], [], []
    for point in points:
        samples = groups[point]
        values = [n.temperature for n in samples]
        if max(values)-min(values) > tolerance:
            raise ValueError(f'conflicting duplicate-coordinate temperatures at {point}: {min(values)} to {max(values)} K')
        temperatures.append(math.fsum(values)/len(values))
        ids.append([n.node_id for n in samples])
        if len(samples) > 1:
            duplicates.append({'xyz_m': list(point), 'node_ids': ids[-1], 'temperature_range_K': max(values)-min(values)})
    return np.array(points), np.array(temperatures), ids, duplicates


def _distance(vertex, shape):
    solver = BRepExtrema_DistShapeShape(vertex, shape)
    if not solver.IsDone() or solver.NbSolution() == 0:
        raise ValueError('exact CAD distance computation failed')
    return solver


def _uv(face, solver):
    closest = solver.PointOnShape2(1)
    projector = GeomAPI_ProjectPointOnSurf(closest, face.surface)
    if projector.NbPoints() == 0:
        raise ValueError(f'cannot project onto face {face.index}')
    uv = list(projector.LowerDistanceParameters())
    umin, umax, vmin, vmax = BRepTools.UVBounds_s(face.shape)
    for axis, periodic, period, lo, hi in (
        (0, face.surface.IsUPeriodic(), face.surface.UPeriod() if face.surface.IsUPeriodic() else 0, umin, umax),
        (1, face.surface.IsVPeriodic(), face.surface.VPeriod() if face.surface.IsVPeriodic() else 0, vmin, vmax),
    ):
        if periodic:
            while uv[axis] < lo-1e-10:
                uv[axis] += period
            while uv[axis] > hi+1e-10:
                uv[axis] -= period
    return uv


def map_samples(points, faces, tolerance, progress=None):
    """Use exact finite-face distances; share samples only across a common CAD edge."""
    bounds = np.array([f.bbox for f in faces])
    membership = [[] for _ in points]
    samples = {f.index: [] for f in faces}
    distances = np.full(len(points), np.inf)
    primary = np.zeros(len(points), dtype=int)
    common_edges = {}
    for a, fa in enumerate(faces):
        for b in range(a+1, len(faces)):
            shared = [ea for ea in fa.edges if any(ea.IsSame(eb) for eb in faces[b].edges)]
            if shared:
                common_edges[a, b] = shared
    for i, point in enumerate(points):
        if progress and i % 5000 == 0:
            progress(f'Mapping source points {i:,}/{len(points):,}')
        lower = np.maximum(np.maximum(bounds[:, 0]-point, point-bounds[:, 1]), 0)
        candidates = np.where(np.linalg.norm(lower, axis=1) <= tolerance)[0]
        vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*point)).Vertex()
        matches = []
        for index in candidates:
            solver = _distance(vertex, faces[index].shape)
            if solver.Value() <= tolerance:
                matches.append((solver.Value(), int(index), solver))
        if not matches:
            continue
        matches.sort(key=lambda x: (x[0], x[1]))
        distance, nearest, _ = matches[0]
        def shares_local_edge(index):
            shared = common_edges.get(tuple(sorted((index, nearest))), [])
            return any(_distance(vertex, edge).Value() <= tolerance for edge in shared)
        # Equally plausible opposite/unrelated surfaces must not be resolved by
        # arbitrary face ordering. A true shared-edge sample can support both.
        if any(abs(other_distance-distance) <= 1e-9 and not shares_local_edge(index)
               for other_distance, index, _ in matches[1:]):
            distances[i] = distance
            continue
        distances[i], primary[i] = distance, faces[nearest].index
        for _, index, solver in matches:
            if index != nearest:
                if not shares_local_edge(index):
                    continue
            uv = _uv(faces[index], solver)
            samples[faces[index].index].append((i, uv))
            membership[i].append(faces[index].index)
    return samples, membership, primary, distances


@dataclass
class Stencil:
    indices: np.ndarray
    weights: np.ndarray
    valid: np.ndarray
    extrapolated: np.ndarray
    nearest_distance: np.ndarray


class FaceInterpolator:
    def __init__(self, face, samples, points, temperatures, duplicate_tolerance=1e-6):
        self.face, self.points = face, points
        # Several distinct source points may project to a single face location.
        # Resolve projection aliases by geometric distance, not temperature.
        grouped = defaultdict(list)
        for index, uv in samples:
            grouped[tuple(np.round(uv, 12))].append(index)
        self.indices, uv = [], []
        self.projected_collisions = []
        for position, ids in grouped.items():
            # Distinct 3D samples can collapse onto the same trimmed edge after
            # projection. They are not duplicate measurements of one source point.
            # Keep the closest sample; reject an equally-close conflicting group.
            distances = np.array([np.linalg.norm(points[i]-face.evaluate(position)) for i in ids])
            best = int(distances.argmin())
            tied = [ids[i] for i in np.where(distances <= distances[best]+1e-9)[0]]
            conflict = np.any(np.ptp(temperatures[:, tied], axis=1) > duplicate_tolerance)
            if len(ids) > 1:
                self.projected_collisions.append({'sample_indices': ids, 'rejected_as_ambiguous': bool(conflict),
                                                 'selected_sample_index': None if conflict else ids[best]})
            if conflict:
                continue
            self.indices.append(ids[best])
            uv.append(position)
        self.indices = np.array(self.indices, dtype=int)
        self.uv = np.array(uv)
        if len(uv) < 3:
            raise ValueError('fewer than three distinct source locations on this face')
        # Use the CAD tangent metric at the sample center, rather than scaling
        # each UV extent to one. Extent normalization distorts long thin panels
        # and can create unnecessarily long interpolation stencils.
        self.origin = self.uv.min(axis=0)
        if np.any(np.ptp(self.uv, axis=0) < 1e-12):
            raise ValueError('source UV coordinates are collinear')
        try:
            self.metric = face.parameter_metric(self.uv.mean(axis=0)) if hasattr(face, 'parameter_metric') else np.eye(2)
        except np.linalg.LinAlgError as exc:
            raise ValueError('singular CAD surface parameter metric') from exc
        try:
            self.triangulation = Delaunay((self.uv-self.origin) @ self.metric)
        except QhullError as exc:
            raise ValueError('source UV coordinates cannot form an interpolation surface') from exc
        self.simplex_ok = np.ones(len(self.triangulation.simplices), dtype=bool)
        # Reject interpolation triangles crossing a trimmed-out region. Multiple
        # edge samples and the centroid guard holes/concave outlines; distance and
        # mesh refinement checks remain necessary for narrow trim features.
        for i, simplex in enumerate(self.triangulation.simplices):
            vertices = self.uv[simplex]
            checks = [vertices.mean(axis=0)]
            for a, b in ((0, 1), (1, 2), (2, 0)):
                checks.extend((1-w)*vertices[a]+w*vertices[b] for w in (.25, .5, .75))
            self.simplex_ok[i] = all(face.inside(p) for p in checks)

    def stencil(self, uv, xyz, max_support, extrapolation_distance):
        from scipy.spatial import cKDTree
        normalized = (uv-self.origin) @ self.metric
        simplex = self.triangulation.find_simplex(normalized, tol=1e-9)
        inside_hull = simplex >= 0
        selected = np.maximum(simplex, 0)
        matrices = self.triangulation.transform[selected]
        first = np.einsum('nij,nj->ni', matrices[:, :2], normalized-matrices[:, 2])
        weights = np.column_stack((first, 1-first.sum(axis=1)))
        indices = self.indices[self.triangulation.simplices[selected]]
        vertex_distance = np.linalg.norm(self.points[indices]-xyz[:, None, :], axis=2).max(axis=1)
        tree = cKDTree(self.points[self.indices])
        nearest_distance, nearest = tree.query(xyz)
        valid = inside_hull & self.simplex_ok[selected] & (vertex_distance <= max_support)
        # Explicit, measured boundary extension, never bridging a rejected hole
        # stencil or borrowing from another CAD face. Its area is separately limited.
        extrapolated = (~inside_hull) & (nearest_distance <= extrapolation_distance)
        for i in np.flatnonzero(extrapolated):
            end = self.uv[nearest[i]]
            path = [uv[i]+w*(end-uv[i]) for w in (0, .25, .5, .75, 1)]
            xyz_path = np.array([self.face.evaluate(p) for p in path])
            path_length = np.linalg.norm(np.diff(xyz_path, axis=0), axis=1).sum()
            if path_length > extrapolation_distance or not all(self.face.inside(p) for p in path):
                extrapolated[i] = False
        indices[extrapolated] = self.indices[nearest[extrapolated], None]
        weights[extrapolated] = (1, 0, 0)
        valid |= extrapolated
        weights = np.maximum(weights, 0)
        weight_sum = weights.sum(axis=1)
        weights[valid] /= weight_sum[valid, None]
        weights[~valid] = 0
        return Stencil(indices, weights, valid, extrapolated, nearest_distance)


def integrate(face, interpolator, temperatures, edge, max_support, extrapolation_distance):
    triangles, uv, xyz, area_weights = quadrature(face, edge)
    if interpolator is None:
        stencil = None
        valid = np.zeros(len(uv), dtype=bool)
        extrapolated = np.zeros(len(uv), dtype=bool)
    else:
        stencil = interpolator.stencil(uv, xyz, max_support, extrapolation_distance)
        valid, extrapolated = stencil.valid, stencil.extrapolated
    area = float(area_weights.sum())
    supported = float(area_weights[valid].sum())
    extrapolated_area = float(area_weights[extrapolated].sum())
    cases = []
    for temperatures_case in temperatures:
        if supported <= 0:
            cases.append(None)
            continue
        values = np.sum(temperatures_case[stencil.indices]*stencil.weights, axis=1)
        mean = float(np.dot(area_weights[valid], values[valid])/supported)
        variance = float(np.dot(area_weights[valid], (values[valid]-mean)**2)/supported)
        cases.append({'mean': mean, 'std': math.sqrt(max(0, variance)),
                      'min': float(values[valid].min()), 'max': float(values[valid].max())})
    return {'triangles': triangles, 'uv': uv, 'xyz': xyz, 'weights': area_weights,
            'stencil': stencil, 'area': area, 'supported_area': supported,
            'coverage': supported/area if area else 0,
            'extrapolated_area': extrapolated_area, 'cases': cases}
