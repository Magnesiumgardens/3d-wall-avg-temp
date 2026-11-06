"""Exact first and second surface moments for linear triangular temperature fields."""

import math

from .io import triangle_area


def summarize(nodes, triangles, panels, assignments, method):
    if method not in ("area", "nodal"):
        raise ValueError("method must be area or nodal")
    element_lookup = {t.element_id: t for t in triangles}
    results = []
    for panel in panels:
        ids = assignments[panel.panel_id]
        if not ids:
            raise ValueError(f"panel {panel.panel_id}: empty assignment")
        if method == "area":
            elements = [element_lookup[key] for key in ids]
            node_ids = set(n for t in elements for n in t.node_ids)
            areas = [triangle_area(t, nodes) for t in elements]
            area = math.fsum(areas)
            if not math.isfinite(area):
                raise ValueError(f"panel {panel.panel_id}: nonfinite total area")
            weights = [a / area for a in areas]
            mean = math.fsum(w * math.fsum(nodes[n].temperature / 3 for n in t.node_ids)
                             for w, t in zip(weights, elements))
            # E[(T-mean)^2] on a linear triangle = (sum(d_i^2)+sum(d_i*d_j))/6.
            # The equivalent expression below avoids cancellation between cross terms.
            moments = []
            for t in elements:
                d = [nodes[n].temperature - mean for n in t.node_ids]
                moments.append((math.fsum(x*x for x in d) + math.fsum(d)**2) / 12)
            variance = math.fsum(w * v for w, v in zip(weights, moments))
            element_count = len(elements)
        else:
            node_ids = set(ids)
            mean = math.fsum(nodes[n].temperature / len(node_ids) for n in node_ids)
            variance = math.fsum((nodes[n].temperature - mean)**2 / len(node_ids) for n in node_ids)
            area, element_count = None, 0
        values = [nodes[n].temperature for n in node_ids]
        if not all(math.isfinite(x) for x in (mean, variance)):
            raise ValueError(f"panel {panel.panel_id}: nonfinite temperature statistics")
        results.append({
            "panel_id": panel.panel_id, "panel_name": panel.name, "area_m2": area,
            "temperature_avg_K": mean, "temperature_min_K": min(values),
            "temperature_max_K": max(values), "temperature_std_K": math.sqrt(variance),
            "node_count": len(node_ids), "element_count": element_count,
            "averaging_method": "area_linear_triangles" if method == "area" else "provisional_nodal_arithmetic",
        })
    return results
