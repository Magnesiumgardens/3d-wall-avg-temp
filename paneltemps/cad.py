"""OpenCascade surface geometry in meters; no SOLIDWORKS installation required."""

from dataclasses import dataclass
import hashlib
import json
import struct

import numpy as np
from OCP.BRep import BRep_Tool
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.BRepTools import BRepTools
from OCP.Bnd import Bnd_Box
from OCP.GProp import GProp_GProps
from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPControl import STEPControl_Reader
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_IN, TopAbs_ON
from OCP.TopExp import TopExp_Explorer
from OCP.TopLoc import TopLoc_Location
from OCP.TopoDS import TopoDS
from OCP.BRepTopAdaptor import BRepTopAdaptor_FClass2d
from OCP.gp import gp_Pnt, gp_Pnt2d, gp_Vec


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coordinate_hash(points):
    """Order-independent hash after exact-coordinate deduplication, in meters."""
    digest = hashlib.sha256()
    for point in sorted({tuple(float(x) if x else 0.0 for x in p) for p in points}):
        digest.update(struct.pack('<ddd', *point))
    return digest.hexdigest()


def load_transform(path, cad_path, points):
    config = json.loads(path.read_text(encoding='utf-8'))
    rotation = np.asarray(config['rotation_BC_to_CAD'], dtype=float)
    translation = np.asarray(config['translation_m'], dtype=float)
    if rotation.shape != (3, 3) or translation.shape != (3,) or not np.isfinite(rotation).all() or not np.isfinite(translation).all():
        raise ValueError('transform requires finite 3x3 rotation and XYZ translation_m')
    if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-8) or not np.isclose(np.linalg.det(rotation), 1, atol=1e-8):
        raise ValueError('transform must be a rigid rotation (no reflection, scale or shear)')
    if config.get('cad_sha256') and config['cad_sha256'] != sha256(cad_path):
        raise ValueError('transform CAD checksum differs: use the exact STEP revision or fit a new transform')
    if config.get('source_coordinates_sha256') and config['source_coordinates_sha256'] != coordinate_hash(points):
        raise ValueError('transform source-coordinate checksum differs: confirm units and the selected BC geometry')
    return points @ rotation.T + translation, config


@dataclass
class Face:
    index: int
    shape: object
    surface: object
    area: float
    center: tuple
    bbox: np.ndarray
    uv_triangles: np.ndarray
    edges: list
    classifier: object

    def inside(self, uv):
        return self.classifier.Perform(gp_Pnt2d(*uv)) in (TopAbs_IN, TopAbs_ON)

    def evaluate(self, uv, derivatives=False):
        if not derivatives:
            return np.array(self.surface.Value(*uv).Coord())
        p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
        self.surface.D1(float(uv[0]), float(uv[1]), p, du, dv)
        return np.array(p.Coord()), np.linalg.norm(np.cross(du.Coord(), dv.Coord()))

    def parameter_metric(self, uv):
        p, du, dv = gp_Pnt(), gp_Vec(), gp_Vec()
        self.surface.D1(float(uv[0]), float(uv[1]), p, du, dv)
        tangent = np.array([du.Coord(), dv.Coord()])
        return np.linalg.cholesky(tangent @ tangent.T)


def load_cad(path, deflection=0.0001):
    reader = STEPControl_Reader()
    if reader.ReadFile(str(path)) != IFSelect_RetDone:
        raise ValueError(f'cannot read STEP: {path}')
    # OpenCascade's setting is measured in mm: 1000 mm per target length unit = m.
    # Set after ReadFile, before TransferRoots; STEP's declared source units remain honored.
    reader.SetSystemLengthUnit(1000.0)
    if reader.TransferRoots() == 0:
        raise ValueError('STEP contains no transferable geometry')
    shape = reader.OneShape()
    if shape.IsNull() or not BRepCheck_Analyzer(shape).IsValid():
        raise ValueError('STEP geometry is empty or invalid')
    mesh = BRepMesh_IncrementalMesh(shape, deflection, False, 0.12, False)
    if not mesh.IsDone():
        raise ValueError('CAD tessellation failed')
    faces = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = TopoDS.Face_s(explorer.Current())
        location = TopLoc_Location()
        triangulation = BRep_Tool.Triangulation_s(face, location)
        if triangulation is None or not triangulation.HasUVNodes():
            raise ValueError('CAD face has no parametric tessellation')
        uv = np.array([triangulation.UVNode(i).Coord() for i in range(1, triangulation.NbNodes()+1)])
        triangles = np.array([triangulation.Triangle(i).Get() for i in range(1, triangulation.NbTriangles()+1)]) - 1
        props = GProp_GProps()
        BRepGProp.SurfaceProperties_s(face, props)
        if not np.isfinite(props.Mass()) or props.Mass() <= 0:
            raise ValueError('CAD face has zero or nonfinite area')
        bounds = Bnd_Box()
        BRepBndLib.AddOptimal_s(face, bounds, False, False)
        box = np.array(bounds.Get()).reshape(2, 3)
        edges = []
        edge_explorer = TopExp_Explorer(face, TopAbs_EDGE)
        while edge_explorer.More():
            edge = edge_explorer.Current()
            if not any(edge.IsSame(other) for other in edges):
                edges.append(edge)
            edge_explorer.Next()
        faces.append(Face(len(faces)+1, face, BRep_Tool.Surface_s(face), props.Mass(),
                          props.CentreOfMass().Coord(), box, uv[triangles], edges,
                          BRepTopAdaptor_FClass2d(face, 1e-9)))
        explorer.Next()
    if not faces:
        raise ValueError('no CAD faces found')
    return faces


def refine_triangles(face, triangles, max_edge, max_triangles=750000):
    """Refine in face UV space using exact surface points to measure edge lengths."""
    pending, completed = np.asarray(triangles), []
    for _ in range(28):
        xyz = np.array([[face.evaluate(uv) for uv in t] for t in pending])
        lengths = np.linalg.norm(xyz - np.roll(xyz, -1, axis=1), axis=2)
        split = lengths.max(axis=1) > max_edge
        completed.extend(pending[~split])
        if not split.any():
            return np.asarray(completed)
        # Longest-edge bisection avoids multiplying the very skinny boundary
        # triangles common in STEP tessellations by four at every refinement.
        selected = pending[split]
        longest = lengths[split].argmax(axis=1)
        row = np.arange(len(selected))
        a = selected[row, longest]
        b = selected[row, (longest+1) % 3]
        c = selected[row, (longest+2) % 3]
        midpoint = (a+b)/2
        pending = np.concatenate([np.stack([a, midpoint, c], axis=1),
                                  np.stack([midpoint, b, c], axis=1)])
        if len(pending) + len(completed) > max_triangles:
            raise ValueError('CAD quadrature refinement exceeded the resource limit; increase --edge-m')
    raise ValueError('CAD quadrature refinement did not converge')


def quadrature(face, max_edge):
    triangles = refine_triangles(face, face.uv_triangles, max_edge)
    # Three-point degree-two Gaussian quadrature on each UV triangle. Jacobians
    # evaluate the exact curved surface, rather than the chordal triangle area.
    barycentric = np.array([[2/3, 1/6, 1/6], [1/6, 2/3, 1/6], [1/6, 1/6, 2/3]])
    uv = np.einsum('qj,tjk->tqk', barycentric, triangles).reshape(-1, 2)
    e1, e2 = triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]
    areas = np.abs(e1[:, 0]*e2[:, 1]-e1[:, 1]*e2[:, 0])/2
    evaluated = [face.evaluate(p, derivatives=True) for p in uv]
    xyz = np.array([p for p, _ in evaluated])
    weights = np.repeat(areas/3, 3)*np.array([j for _, j in evaluated])
    inside = np.array([face.inside(p) for p in uv])
    weights[~inside] = 0
    return triangles, uv, xyz, weights
