"""Offline, interactive views for checking panel coverage and averages."""

from html import escape
from pathlib import Path


def write_visualization(path, nodes, triangles, panels, assignments, results, method):
    """Write two 3D views with Plotly embedded, requiring no network to open.

    The first view shows the requested panel assignments; the second colors the
    same geometry by its panel average. Unassigned geometry is visible in gray.
    In ``area`` mode these are the supplied triangles, not a triangulation or a
    surface reconstruction from point locations.
    """
    if method not in {"area", "nodal"}:
        raise ValueError("Visualization method must be 'area' or 'nodal'.")
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError as exc:
        raise RuntimeError(
            "Plotly is required for visualization; install requirements.txt."
        ) from exc

    palette = [
        "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
        "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
        "#2E91E5", "#E15F99", "#1CA71C", "#FB0D0D", "#DA16FF",
        "#222A2A", "#B68100", "#750D86", "#EB663B", "#511CFB",
    ]
    summaries = {str(row["panel_id"]): row for row in results}
    triangles_by_id = {triangle.element_id: triangle for triangle in triangles}
    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{"type": "scene"}, {"type": "scene"}]],
        subplot_titles=("Panel assignments", "Panel-average temperature (K)"),
        horizontal_spacing=0.02,
    )

    def points(node_ids):
        selected = [nodes[node_id] for node_id in node_ids]
        return dict(
            x=[node.x for node in selected],
            y=[node.y for node in selected],
            z=[node.z for node in selected],
        )

    def mesh(element_ids):
        # Repeated triangle vertices let hover text identify the element while
        # preserving shared panel boundaries without changing any topology.
        ids = [
            node_id
            for element_id in element_ids
            for node_id in triangles_by_id[element_id].node_ids
        ]
        geometry = points(ids)
        geometry.update(
            i=list(range(0, len(ids), 3)),
            j=list(range(1, len(ids), 3)),
            k=list(range(2, len(ids), 3)),
        )
        return geometry, ids

    def hover_suffix():
        return "<br>x=%{x:.6g} m<br>y=%{y:.6g} m<br>z=%{z:.6g} m<extra></extra>"

    def add_geometry(ids, label, group, color, average=None, hover=None):
        safe_label = escape(str(label), quote=True)
        if method == "area":
            geometry, vertex_ids = mesh(ids)
            element_labels = [element_id for element_id in ids for _ in range(3)]
            details = [
                f"{hover or safe_label}<br>Element {escape(str(element_id))}"
                f"<br>Node {node_id}"
                for element_id, node_id in zip(element_labels, vertex_ids)
            ]
            for column in (1, 2):
                styling = dict(color=color)
                if column == 2 and average is not None:
                    styling = dict(
                        intensity=[average] * len(vertex_ids),
                        intensitymode="vertex", coloraxis="coloraxis",
                    )
                fig.add_trace(
                    go.Mesh3d(
                        **geometry, **styling, name=safe_label,
                        legendgroup=group, showlegend=column == 1,
                        text=details,
                        hovertemplate="%{text}" + hover_suffix(),
                        flatshading=True,
                        lighting=dict(ambient=0.85, diffuse=0.4, specular=0.05),
                    ), row=1, col=column,
                )
        else:
            geometry = points(ids)
            details = [f"{hover or safe_label}<br>Node {node_id}" for node_id in ids]
            for column in (1, 2):
                marker = dict(size=3, color=color, opacity=0.9)
                if column == 2 and average is not None:
                    marker.update(color=[average] * len(ids), coloraxis="coloraxis")
                fig.add_trace(
                    go.Scatter3d(
                        **geometry, mode="markers", marker=marker,
                        name=safe_label, legendgroup=group, showlegend=column == 1,
                        text=details, hovertemplate="%{text}" + hover_suffix(),
                    ), row=1, col=column,
                )

    used_ids = set()
    averages = []
    for index, panel in enumerate(panels):
        panel_id = str(panel.panel_id)
        ids = list(assignments[panel_id])
        summary = summaries[panel_id]
        average = float(summary["temperature_avg_K"])
        averages.append(average)
        label = f"{panel_id}: {panel.name}"
        hover = (
            f"{escape(label, quote=True)}"
            f"<br>Average: {average:.6g} K"
            f"<br>Range: {summary['temperature_min_K']:.6g}"
            f" – {summary['temperature_max_K']:.6g} K"
        )
        add_geometry(
            ids, label, f"panel-{index}", palette[index % len(palette)],
            average=average, hover=hover,
        )
        used_ids.update(ids)

    available_ids = set(triangles_by_id) if method == "area" else set(nodes)
    unmapped = sorted(available_ids - used_ids)
    if unmapped:
        add_geometry(unmapped, "Unassigned", "unassigned", "#999999")

    # In area mode nodes not participating in any triangle also need to remain
    # visible. These cannot establish surface coverage or a surface area.
    orphan_nodes = set()
    if method == "area":
        connected_nodes = {
            node_id for triangle in triangles for node_id in triangle.node_ids
        }
        orphan_nodes = set(nodes) - connected_nodes
        if orphan_nodes:
            ids = sorted(orphan_nodes)
            for column in (1, 2):
                fig.add_trace(
                    go.Scatter3d(
                        **points(ids), mode="markers",
                        marker=dict(size=3, color="#999999"),
                        name="Nodes without elements", legendgroup="orphan-nodes",
                        showlegend=column == 1,
                        text=[f"Node {node_id}: no surface element" for node_id in ids],
                        hovertemplate="%{text}" + hover_suffix(),
                    ), row=1, col=column,
                )

    if averages:
        lower, upper = min(averages), max(averages)
        if lower == upper:
            lower, upper = max(0.0, lower - 0.5), upper + 0.5
    else:
        lower, upper = 0.0, 1.0
    scene = dict(
        xaxis_title="x (m)", yaxis_title="y (m)", zaxis_title="z (m)",
        aspectmode="data",
    )
    fig.update_layout(
        scene=scene, scene2=scene,
        coloraxis=dict(
            colorscale="Turbo", cmin=lower, cmax=upper,
            colorbar=dict(title="Average (K)", x=1.02, len=0.65),
        ),
        margin=dict(l=0, r=90, t=55, b=110), height=720,
        legend=dict(orientation="h", y=-0.06, x=0, groupclick="togglegroup"),
        template="plotly_white",
    )
    note = (
        "Area mode: the supplied triangular elements are shown. The temperature "
        "view displays one uniform average per panel, not the original nodal field."
        if method == "area" else
        "PROVISIONAL NODAL MODE: these arithmetic averages are sensitive to mesh "
        "density. Points do not define surface areas or prove surface coverage."
    )
    entity = "elements" if method == "area" else "nodes"
    coverage = f"Unassigned {entity}: {len(unmapped)}."
    if method == "area":
        coverage += f" Nodes without elements: {len(orphan_nodes)}."
    fragment = fig.to_html(
        full_html=False, include_plotlyjs=True,
        config={"responsive": True, "displaylogo": False},
        div_id="panel-temperature-views",
    )
    document = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<title>Panel assignments and average temperatures</title>'
        '<style>body{font-family:system-ui,sans-serif;margin:1.5rem;}'
        'p{max-width:80rem;line-height:1.5}.note{background:#f3f4f6;'
        'padding:1rem;border-radius:.4rem}</style></head><body>'
        '<h1>Panel assignments and average temperatures</h1>'
        f'<p class="note">{escape(note)} {escape(coverage)}</p>'
        '<p>Drag each view to rotate it; scroll to zoom. Use the legend to toggle '
        'panels in both views. Distances are in meters. Confirm the coordinate '
        'frame and physical regions against the source geometry and CAD model.</p>'
        + fragment + '</body></html>'
    )
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(document, encoding="utf-8")
