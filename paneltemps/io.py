"""Strict readers for six-column nodal data and explicit triangle connectivity."""

import csv
import math
import re
from pathlib import Path

from .model import Node, Triangle

COORDINATE_SCALES = {"m": 1.0, "mm": 0.001, "cm": 0.01, "in": 0.0254, "ft": 0.3048}


def _rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = next(csv.reader([line])) if "," in line else line.split()
            yield line_number, [field.strip() for field in fields]


def _header(fields):
    return [re.sub(r"[^a-z0-9]", "", x.lower()) for x in fields]


def _identifier(value):
    if not re.fullmatch(r"[0-9]+", value) or int(value) <= 0:
        raise ValueError(f"expected a positive integer ID, got {value!r}")
    return int(value)


def _number(value):
    number = float(value.replace("D", "E").replace("d", "e"))
    if not math.isfinite(number):
        raise ValueError(f"expected a finite number, got {value!r}")
    return number


def read_nodes(path, coordinate_unit="m", temperature_unit="K", temperature_column=None):
    if coordinate_unit not in COORDINATE_SCALES or temperature_unit not in ("K", "C", "F"):
        raise ValueError("unsupported coordinate or temperature unit")
    scale = COORDINATE_SCALES[coordinate_unit]
    nodes = {}
    first = True
    indices = [0, 1, 2, 3, 4, 5]
    width = 6
    if temperature_column not in (None, "temperature", "wall-temperature"):
        raise ValueError("temperature_column must be temperature or wall-temperature")
    for line_number, fields in _rows(path):
        if first:
            first = False
            header = _header(fields)
            if "nodenumber" in header:
                if len(header) != len(set(header)):
                    raise ValueError(f"{path}: duplicate column names")
                chosen = _header([temperature_column])[0] if temperature_column else (
                    "walltemperature" if "walltemperature" in header else "temperature")
                required = ["nodenumber", "xcoordinate", "ycoordinate", "zcoordinate", "pressure", chosen]
                if not set(required) <= set(header):
                    raise ValueError(f"{path}: missing required named columns {sorted(set(required)-set(header))}")
                indices = [header.index(key) for key in required]
                width = len(header)
                continue
            if temperature_column == "temperature":
                raise ValueError("explicit temperature column requires a named header")
        try:
            if len(fields) != width:
                raise ValueError(f"expected {width} fields, found {len(fields)}")
            selected = [fields[i] for i in indices]
            node_id = _identifier(selected[0])
            if node_id in nodes:
                raise ValueError(f"duplicate node ID {node_id}")
            # Validate auxiliary CFD columns too; originals retain all fields.
            for i, value in enumerate(fields):
                if i != indices[0]:
                    _number(value)
            x, y, z, pressure, temperature = map(_number, selected[1:])
            x, y, z = x * scale, y * scale, z * scale
            if temperature_unit == "C":
                temperature += 273.15
            elif temperature_unit == "F":
                temperature = (temperature - 32.0) * 5.0 / 9.0 + 273.15
            if not all(math.isfinite(v) for v in (x, y, z, temperature)):
                raise ValueError("unit conversion produced a nonfinite value")
            if temperature < 0:
                raise ValueError("temperature is below absolute zero")
            nodes[node_id] = Node(node_id, x, y, z, pressure, temperature)
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: {exc}") from exc
    if not nodes:
        raise ValueError(f"{path}: no nodal data")
    return nodes


def triangle_area(triangle, nodes):
    a, b, c = (nodes[node_id].xyz for node_id in triangle.node_ids)
    u = tuple(b[i] - a[i] for i in range(3))
    v = tuple(c[i] - a[i] for i in range(3))
    cross = (u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0])
    area = math.hypot(*cross) / 2.0
    if not math.isfinite(area) or area <= 0:
        raise ValueError(f"element {triangle.element_id} has zero or nonfinite area")
    return area


def read_triangles(path, nodes):
    triangles, element_ids, faces = [], set(), set()
    first = True
    for line_number, fields in _rows(path):
        if first:
            first = False
            normalized = _header(fields)
            if normalized in (["elementid", "node1", "node2", "node3"],
                              ["elementid", "node1", "node2", "node3", "panelid"]):
                continue
        try:
            if len(fields) not in (4, 5):
                raise ValueError("expected 4 or 5 connectivity fields")
            element_id = _identifier(fields[0])
            node_ids = tuple(_identifier(x) for x in fields[1:4])
            if element_id in element_ids:
                raise ValueError(f"duplicate element ID {element_id}")
            if len(set(node_ids)) != 3:
                raise ValueError("triangle must have three distinct nodes")
            if not set(node_ids) <= nodes.keys():
                raise ValueError(f"unknown node IDs {sorted(set(node_ids) - nodes.keys())}")
            face = tuple(sorted(node_ids))
            if face in faces:
                raise ValueError(f"duplicate triangle {face}")
            triangle = Triangle(element_id, node_ids, fields[4] or None if len(fields) == 5 else None)
            triangle_area(triangle, nodes)
            triangles.append(triangle)
            element_ids.add(element_id)
            faces.add(face)
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: {exc}") from exc
    if not triangles:
        raise ValueError(f"{path}: no surface triangles")
    return triangles
