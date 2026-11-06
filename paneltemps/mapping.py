"""Explicit panel membership. Bounds require all triangle vertices inside."""

import json
import math
from pathlib import Path

from .model import Panel


def load_panels(path):
    config = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if not isinstance(config, dict) or set(config) != {"panels"}:
        raise ValueError("mapping must contain only a 'panels' array; bounds use meters")
    if not isinstance(config["panels"], list) or not config["panels"]:
        raise ValueError("mapping requires at least one panel")
    panels, seen = [], set()
    for item in config["panels"]:
        if not isinstance(item, dict):
            raise ValueError("each panel must be an object")
        allowed = {"id", "name", "solidworks_component", "solidworks_region", "selection_notes", "selector"}
        if set(item) - allowed:
            raise ValueError(f"unknown panel fields: {sorted(set(item) - allowed)}")
        panel_id, name, selector = item.get("id"), item.get("name"), item.get("selector")
        if not isinstance(panel_id, str) or not panel_id.strip() or panel_id in seen:
            raise ValueError("panel IDs must be unique nonempty strings")
        if not isinstance(name, str) or not name.strip():
            raise ValueError(f"panel {panel_id}: name must be a nonempty string")
        if not isinstance(selector, dict):
            raise ValueError(f"panel {panel_id}: selector must be an object")
        kind = selector.get("type")
        if kind in ("element_ids", "node_ids"):
            ids = selector.get("ids")
            if set(selector) != {"type", "ids"} or not isinstance(ids, list) or not ids:
                raise ValueError(f"panel {panel_id}: selector requires a nonempty ids array")
            if any(type(x) is not int or x <= 0 for x in ids) or len(ids) != len(set(ids)):
                raise ValueError(f"panel {panel_id}: IDs must be unique positive integers")
        elif kind == "bounds":
            if set(selector) != {"type", "min", "max"}:
                raise ValueError(f"panel {panel_id}: bounds require min and max arrays")
            for key in ("min", "max"):
                values = selector[key]
                if not isinstance(values, list) or len(values) != 3 or any(
                    type(x) not in (int, float) or not math.isfinite(x) for x in values
                ):
                    raise ValueError(f"panel {panel_id}: bounds must be finite XYZ triples in meters")
            if any(lo > hi for lo, hi in zip(selector["min"], selector["max"])):
                raise ValueError(f"panel {panel_id}: inverted bounds")
        elif kind == "mesh_panel":
            if set(selector) != {"type", "value"} or not isinstance(selector.get("value"), str) or not selector["value"]:
                raise ValueError(f"panel {panel_id}: mesh_panel requires a nonempty string value")
        else:
            raise ValueError(f"panel {panel_id}: unsupported selector type {kind!r}")
        metadata = {key: item.get(key, "") for key in ("solidworks_component", "solidworks_region", "selection_notes")}
        if any(not isinstance(value, str) for value in metadata.values()):
            raise ValueError(f"panel {panel_id}: CAD metadata must be strings")
        panels.append(Panel(panel_id, name, selector=selector, **metadata))
        seen.add(panel_id)
    return panels


def assign_panels(nodes, triangles, panels, method, allow_unmapped=False):
    if method not in ("area", "nodal"):
        raise ValueError("method must be area or nodal")
    if method == "area" and not triangles:
        raise ValueError("area averaging requires triangle connectivity")
    entities = {t.element_id: t for t in triangles} if method == "area" else nodes
    assignments, owners = {}, {}
    for panel in panels:
        selector = panel.selector
        kind = selector["type"]
        if kind in ("node_ids", "element_ids"):
            expected_kind = "element_ids" if method == "area" else "node_ids"
            if kind != expected_kind:
                raise ValueError(f"panel {panel.panel_id}: {method} requires {expected_kind}, bounds, or compatible mesh labels")
            missing = set(selector["ids"]) - entities.keys()
            if missing:
                raise ValueError(f"panel {panel.panel_id}: unknown IDs {sorted(missing)}")
            selected = list(selector["ids"])
        elif kind == "mesh_panel":
            if method != "area":
                raise ValueError("mesh_panel selectors require area mode with triangle panel_id labels")
            selected = [key for key, entity in entities.items() if entity.panel_id == selector["value"]]
        elif kind == "bounds":
            selected = []
            for key, entity in entities.items():
                points = [nodes[n].xyz for n in entity.node_ids] if method == "area" else [entity.xyz]
                if all(all(lo <= x <= hi for lo, x, hi in zip(selector["min"], xyz, selector["max"]))
                       for xyz in points):
                    selected.append(key)
        else:
            raise ValueError(f"unknown selector {kind!r}")
        if not selected:
            raise ValueError(f"panel {panel.panel_id}: no contributing {'elements' if method == 'area' else 'nodes'}")
        for entity_id in selected:
            if entity_id in owners:
                raise ValueError(f"overlap at ID {entity_id}: panels {owners[entity_id]} and {panel.panel_id}")
            owners[entity_id] = panel.panel_id
        assignments[panel.panel_id] = selected
    missing = set(entities) - owners.keys()
    if missing and not allow_unmapped:
        raise ValueError(f"{len(missing)} unmapped {'elements' if method == 'area' else 'nodes'} (first IDs: {sorted(missing)[:10]}). Define panels or explicitly use --allow-unmapped")
    return assignments
