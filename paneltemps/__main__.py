"""Command-line entry point: python -m paneltemps --help."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys

from . import __version__
from .io import COORDINATE_SCALES, read_nodes, read_triangles
from .mapping import assign_panels, load_panels
from .stats import summarize
from .visualize import write_visualization


def _csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be a finite positive number")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description="Export per-panel surface temperatures and an offline 3D mapping check.")
    parser.add_argument("input", type=Path, help="named-column CFD export or six-column CSV/whitespace file")
    parser.add_argument("--panels", type=Path, required=True, help="explicit JSON panel mapping")
    parser.add_argument("--elements", type=Path, help="linear triangle connectivity CSV")
    parser.add_argument("--method", choices=("area", "nodal"), default="area",
                        help="default: area, which requires --elements; nodal is explicitly provisional")
    parser.add_argument("--coordinate-unit", choices=tuple(COORDINATE_SCALES), default="m")
    parser.add_argument("--temperature-unit", choices=("K", "C", "F"), default="K")
    parser.add_argument("--temperature-column", choices=("temperature", "wall-temperature"),
                        help="named source field; default prefers wall-temperature, otherwise temperature")
    parser.add_argument("--wall-temperature-kind", choices=("unknown", "prescribed", "adiabatic", "structural"), default="unknown")
    parser.add_argument("--case-label", default="unspecified", help="source case identity, e.g. synthetic_demo")
    parser.add_argument("--allow-unmapped", action="store_true", help="allow partial coverage; report every omitted ID")
    parser.add_argument("--warn-range-K", type=_positive, help="flag panel temperature ranges at or above this threshold")
    parser.add_argument("--output-dir", type=Path, required=True, help="new or empty directory; existing outputs are never overwritten")
    args = parser.parse_args(argv)
    try:
        if args.method == "area" and args.elements is None:
            raise ValueError("area method requires --elements; select --method nodal explicitly for provisional arithmetic averages")
        if args.output_dir.exists() and (not args.output_dir.is_dir() or any(args.output_dir.iterdir())):
            raise ValueError("--output-dir must be a new or empty directory")
        nodes = read_nodes(args.input, args.coordinate_unit, args.temperature_unit, args.temperature_column)
        triangles = read_triangles(args.elements, nodes) if args.elements else []
        panels = load_panels(args.panels)
        assignments = assign_panels(nodes, triangles, panels, args.method, args.allow_unmapped)
        results = summarize(nodes, triangles, panels, assignments, args.method)

        used = {key for ids in assignments.values() for key in ids}
        all_ids = {t.element_id for t in triangles} if args.method == "area" else set(nodes)
        unmapped = sorted(all_ids - used)
        connected = {n for t in triangles for n in t.node_ids}
        orphans = sorted(set(nodes) - connected) if args.method == "area" else []
        cad_unmapped = [p.panel_id for p in panels if not p.solidworks_component.strip() or not p.solidworks_region.strip()]
        warnings = [
            "Source units are user-specified assumptions, not independently verified.",
            "CAD face correspondence and exterior-skin coverage have not been verified by this tool.",
            "Area-averaging temperature does not preserve heat flux or total heat input.",
        ]
        if args.wall_temperature_kind == "unknown":
            warnings.append("Wall-temperature physical definition is unknown; confirm it with the source-data owner.")
        elif args.wall_temperature_kind in ("prescribed", "adiabatic"):
            warnings.append("CFD prescribed or adiabatic wall temperatures are not independently predicted conducting-skin temperatures.")
        if args.method == "nodal":
            warnings.append("PROVISIONAL: arithmetic nodal averages depend on mesh density; area is unavailable.")
        if any(p.selector["type"] == "bounds" for p in panels):
            warnings.append("Bounds select whole triangles only when every vertex is inside; crossing triangles are not clipped.")
        if unmapped:
            warnings.append(f"Partial coverage explicitly allowed: {len(unmapped)} unmapped {'elements' if args.method == 'area' else 'nodes'}.")
        if orphans:
            warnings.append(f"{len(orphans)} nodes have no supplied surface elements; their surface coverage cannot be assessed.")
        if cad_unmapped:
            warnings.append(f"{len(cad_unmapped)} panels lack a complete CAD component/region mapping.")
        ranges = [{"panel_id": r["panel_id"], "range_K": r["temperature_max_K"] - r["temperature_min_K"]} for r in results]
        flagged = [r["panel_id"] for r in ranges if args.warn_range_K is not None and r["range_K"] >= args.warn_range_K]
        if flagged:
            warnings.append("Review temperature gradients and consider justified smaller regions for: " + ", ".join(flagged))

        args.output_dir.mkdir(parents=True, exist_ok=True)
        # Write visualization first, so a missing dependency does not leave seemingly complete tables.
        write_visualization(args.output_dir / "panel_visualization.html", nodes, triangles, panels, assignments, results, args.method)
        _csv(args.output_dir / "panel_temperatures.csv", list(results[0]), results)
        mapping_fields = ["panel_id", "panel_name", "solidworks_component", "solidworks_region", "selection_notes"]
        _csv(args.output_dir / "solidworks_mapping.csv", mapping_fields, [
            dict(panel_id=p.panel_id, panel_name=p.name, solidworks_component=p.solidworks_component,
                 solidworks_region=p.solidworks_region, selection_notes=p.selection_notes) for p in panels
        ])
        entity_key = "element_ids" if args.method == "area" else "node_ids"
        (args.output_dir / "panel_assignments.json").write_text(json.dumps(
            {"entity_type": entity_key, "assignments": assignments, "unmapped_ids": unmapped,
             "nodes_without_elements": orphans}, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        inputs_dir = args.output_dir / "inputs"
        inputs_dir.mkdir()
        preserved = []
        for role, source in (("boundary_conditions", args.input), ("panel_mapping", args.panels), ("connectivity", args.elements)):
            if source is None:
                continue
            destination = inputs_dir / (role + source.suffix)
            shutil.copyfile(source, destination)
            preserved.append({"role": role, "source_name": source.name,
                              "preserved_path": str(destination.relative_to(args.output_dir)),
                              "sha256": hashlib.sha256(destination.read_bytes()).hexdigest()})
        report = {
            "tool_version": __version__, "case_label": args.case_label, "method": args.method,
            "input_coordinate_unit": args.coordinate_unit, "input_temperature_unit": args.temperature_unit,
            "temperature_column_selection": args.temperature_column or "auto: wall-temperature preferred, otherwise temperature",
            "output_coordinate_unit": "m", "output_temperature_unit": "K",
            "wall_temperature_kind": args.wall_temperature_kind,
            "input_node_count": len(nodes), "input_element_count": len(triangles), "panel_count": len(panels),
            "unmapped_entity_ids": unmapped, "nodes_without_elements": orphans,
            "panels_without_complete_cad_mapping": cad_unmapped,
            "cad_mapping_verified": False, "source_units_verified": False,
            "panel_temperature_ranges": ranges, "range_warning_threshold_K": args.warn_range_K,
            "panels_flagged_for_range_review": flagged, "warnings": warnings, "preserved_inputs": preserved,
            "standard_deviation_definition": "population surface standard deviation of a piecewise-linear triangular field" if args.method == "area" else "population nodal standard deviation (provisional)",
        }
        (args.output_dir / "run_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError, RuntimeError, OverflowError, csv.Error) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(f"Exported {len(results)} panels to {args.output_dir}")
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
