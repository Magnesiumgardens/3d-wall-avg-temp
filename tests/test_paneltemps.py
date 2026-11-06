"""Behavioral checks for parsing, panel coverage, and exact surface integration."""

import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from paneltemps.io import read_nodes, read_triangles
from paneltemps.mapping import assign_panels, load_panels
from paneltemps.model import Node, Triangle
from paneltemps.stats import summarize


HEADER = "nodenumber,x-coordinate,y-coordinate,z-coordinate,pressure,wall-temperature\n"
ELEMENT_HEADER = "element_id,node_1,node_2,node_3,panel_id\n"


class TemporaryInputs(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding="utf-8")
        return path

    def panels(self, *definitions):
        return load_panels(self.write("panels.json", json.dumps({"panels": definitions})))

    @staticmethod
    def triangle_nodes():
        return {
            1: Node(1, 0.0, 0.0, 0.0, 101.0, 300.0),
            2: Node(2, 2.0, 0.0, 0.0, 102.0, 300.0),
            3: Node(3, 0.0, 1.0, 0.0, 103.0, 303.0),
        }


class NodeParsingTests(TemporaryInputs):
    def test_named_isothermal_export_keeps_temperature_separate_from_flux(self):
        path = self.write("iso", "nodenumber,x-coordinate,y-coordinate,z-coordinate,pressure,temperature,heat-flux\n1,0,0,0,100,300,-10000\n")
        self.assertEqual(read_nodes(path)[1].temperature, 300)
        with self.assertRaises(ValueError):
            read_nodes(path, temperature_column="wall-temperature")

    def test_named_adiabatic_export_prefers_wall_temperature(self):
        path = self.write("adiab", "nodenumber,x-coordinate,y-coordinate,z-coordinate,pressure,dynamic-pressure,temperature,wall-temperature\n1,0,0,0,100,4,2200,2300\n")
        self.assertEqual(read_nodes(path)[1].temperature, 2300)
        self.assertEqual(read_nodes(path, temperature_column="temperature")[1].temperature, 2200)

    def test_named_columns_allow_reordering_but_reject_duplicate_headers(self):
        path = self.write("reordered", "temperature,pressure,z-coordinate,y-coordinate,x-coordinate,nodenumber\n300,100,3,2,1,5\n")
        self.assertEqual(read_nodes(path)[5].xyz, (1, 2, 3))
        bad = self.write("duplicate", "nodenumber,x-coordinate,y-coordinate,z-coordinate,pressure,temperature,temperature\n1,0,0,0,100,300,301\n")
        with self.assertRaises(ValueError):
            read_nodes(bad)

    def test_comma_scientific_and_fortran_notation(self):
        path = self.write("nodes.csv", HEADER + "12,1.642018079E+00,-6.75E-02,-2.448D-01,1.0175E+04,1.843398156E+03\n")
        nodes = read_nodes(path)
        self.assertEqual(list(nodes), [12])
        self.assertAlmostEqual(nodes[12].x, 1.642018079)
        self.assertAlmostEqual(nodes[12].y, -0.0675)
        self.assertAlmostEqual(nodes[12].z, -0.2448)
        self.assertAlmostEqual(nodes[12].pressure, 10175.0)
        self.assertAlmostEqual(nodes[12].temperature, 1843.398156)

    def test_whitespace_file_and_unit_conversion(self):
        path = self.write("nodes.dat", HEADER.replace(",", " ") + "1\t1000  -20  500   10  26.85\n")
        node = read_nodes(path, coordinate_unit="mm", temperature_unit="C")[1]
        self.assertEqual((node.x, node.y, node.z), (1.0, -0.02, 0.5))
        self.assertAlmostEqual(node.temperature, 300.0)

    def test_invalid_values_and_duplicate_ids_fail(self):
        cases = {
            "duplicate": "1,0,0,0,10,300\n1,1,0,0,10,301\n",
            "fractional ID": "1.5,0,0,0,10,300\n",
            "zero ID": "0,0,0,0,10,300\n",
            "nan coordinate": "1,nan,0,0,10,300\n",
            "infinite temperature": "1,0,0,0,10,inf\n",
            "negative Kelvin": "1,0,0,0,10,-1\n",
            "wrong column count": "1,0,0,10,300\n",
            "empty data": "",
        }
        for label, rows in cases.items():
            with self.subTest(label=label), self.assertRaises(ValueError):
                read_nodes(self.write("bad.csv", HEADER + rows))

    def test_celsius_below_absolute_zero_fails(self):
        with self.assertRaises(ValueError):
            read_nodes(self.write("cold.csv", HEADER + "1,0,0,0,0,-274\n"), temperature_unit="C")


class ElementParsingTests(TemporaryInputs):
    def test_triangle_connectivity_and_supplied_panel_label(self):
        triangles = read_triangles(self.write("surface.csv", ELEMENT_HEADER + "20,1,2,3,skin_a\n"), self.triangle_nodes())
        self.assertEqual(len(triangles), 1)
        self.assertEqual(triangles[0].element_id, 20)
        self.assertEqual(tuple(triangles[0].node_ids), (1, 2, 3))
        self.assertEqual(triangles[0].panel_id, "skin_a")

    def test_invalid_connectivity_fails(self):
        cases = {
            "missing node": "20,1,2,99,a\n",
            "repeated node": "20,1,2,2,a\n",
            "duplicate ID": "20,1,2,3,a\n20,1,3,2,b\n",
            "duplicate geometry": "20,1,2,3,a\n21,3,2,1,b\n",
        }
        for label, rows in cases.items():
            with self.subTest(label=label), self.assertRaises(ValueError):
                read_triangles(self.write("bad_elements.csv", ELEMENT_HEADER + rows), self.triangle_nodes())

    def test_collinear_triangle_fails(self):
        nodes = self.triangle_nodes()
        nodes[3] = Node(3, 1.0, 0.0, 0.0, 0.0, 303.0)
        with self.assertRaises(ValueError):
            read_triangles(self.write("flat.csv", ELEMENT_HEADER + "20,1,2,3,a\n"), nodes)


class MappingTests(TemporaryInputs):
    def test_explicit_ids_and_mesh_labels(self):
        nodes = self.triangle_nodes()
        elements = [Triangle(20, (1, 2, 3), "skin_a")]
        panels = self.panels({"id": "a", "name": "Skin A", "selector": {"type": "mesh_panel", "value": "skin_a"}})
        self.assertEqual(assign_panels(nodes, elements, panels, "area"), {"a": [20]})
        panels = self.panels({"id": "a", "name": "Skin A", "selector": {"type": "element_ids", "ids": [20]}})
        self.assertEqual(assign_panels(nodes, elements, panels, "area"), {"a": [20]})

    def test_node_mapping_does_not_allow_silent_omission(self):
        panels = self.panels({"id": "a", "name": "Partial", "selector": {"type": "node_ids", "ids": [1, 2]}})
        with self.assertRaises(ValueError):
            assign_panels(self.triangle_nodes(), [], panels, "nodal")
        self.assertEqual(assign_panels(self.triangle_nodes(), [], panels, "nodal", allow_unmapped=True), {"a": [1, 2]})

    def test_overlapping_selections_fail(self):
        panels = self.panels(
            {"id": "a", "name": "A", "selector": {"type": "node_ids", "ids": [1, 2]}},
            {"id": "b", "name": "B", "selector": {"type": "node_ids", "ids": [2, 3]}},
        )
        with self.assertRaises(ValueError):
            assign_panels(self.triangle_nodes(), [], panels, "nodal")

    def test_unknown_ids_and_empty_regions_fail(self):
        for selector in [
            {"type": "node_ids", "ids": [1, 2, 3, 99]},
            {"type": "bounds", "min": [9, 9, 9], "max": [10, 10, 10]},
        ]:
            with self.subTest(selector=selector), self.assertRaises(ValueError):
                panels = self.panels({"id": "a", "name": "A", "selector": selector})
                assign_panels(self.triangle_nodes(), [], panels, "nodal", allow_unmapped=True)

    def test_bounds_are_in_meters(self):
        panels = self.panels({"id": "a", "name": "Whole", "selector": {"type": "bounds", "min": [-0.1, -0.1, -0.1], "max": [2.1, 1.1, 0.1]}})
        self.assertEqual(set(assign_panels(self.triangle_nodes(), [], panels, "nodal")["a"]), {1, 2, 3})

    def test_area_bounds_leave_crossing_triangles_unmapped(self):
        nodes = {
            1: Node(1, 0, 0, 0, 0, 300),
            2: Node(2, 1, 0, 0, 0, 300),
            3: Node(3, 1, 1, 0, 0, 300),
            4: Node(4, 1.2, 0, 0, 0, 300),
        }
        elements = [Triangle(1, (1, 2, 3)), Triangle(2, (1, 4, 3))]
        panels = self.panels({"id": "left", "name": "Left", "selector": {"type": "bounds", "min": [0, 0, 0], "max": [1, 1, 0]}})
        # Element 2's centroid is inside, but one vertex crosses the boundary.
        with self.assertRaises(ValueError):
            assign_panels(nodes, elements, panels, "area")
        self.assertEqual(assign_panels(nodes, elements, panels, "area", allow_unmapped=True), {"left": [1]})

    def test_area_panels_may_share_boundary_nodes_without_element_overlap(self):
        nodes = {
            1: Node(1, 0, 0, 0, 0, 300),
            2: Node(2, 1, 0, 0, 0, 301),
            3: Node(3, 1, 1, 0, 0, 302),
            4: Node(4, 2, 0, 0, 0, 303),
        }
        elements = [Triangle(1, (1, 2, 3)), Triangle(2, (2, 4, 3))]
        panels = self.panels(
            {"id": "left", "name": "Left", "selector": {"type": "bounds", "min": [0, 0, 0], "max": [1, 1, 0]}},
            {"id": "right", "name": "Right", "selector": {"type": "bounds", "min": [1, 0, 0], "max": [2, 1, 0]}},
        )
        self.assertEqual(assign_panels(nodes, elements, panels, "area"), {"left": [1], "right": [2]})


class IntegrationTests(TemporaryInputs):
    def area_summary(self, nodes, elements):
        panels = self.panels({"id": "a", "name": "Skin A", "selector": {"type": "element_ids", "ids": [t.element_id for t in elements]}})
        assignments = assign_panels(nodes, elements, panels, "area")
        return summarize(nodes, elements, panels, assignments, "area")[0]

    def test_exact_linear_triangle_moments(self):
        result = self.area_summary(self.triangle_nodes(), [Triangle(20, (1, 2, 3))])
        self.assertAlmostEqual(result["area_m2"], 1.0)
        self.assertAlmostEqual(result["temperature_avg_K"], 301.0)
        self.assertAlmostEqual(result["temperature_std_K"], math.sqrt(0.5))
        self.assertEqual(result["temperature_min_K"], 300.0)
        self.assertEqual(result["temperature_max_K"], 303.0)
        self.assertEqual(result["node_count"], 3)
        self.assertEqual(result["element_count"], 1)

    def test_nonuniform_refinement_preserves_mean_area_and_std(self):
        nodes = self.triangle_nodes()
        original = self.area_summary(nodes, [Triangle(20, (1, 2, 3))])
        nodes[4] = Node(4, 0.2, 0.1, 0.0, 0.0, 300.3)
        refined = self.area_summary(nodes, [Triangle(21, (1, 2, 4)), Triangle(22, (2, 3, 4)), Triangle(23, (3, 1, 4))])
        for field in ["area_m2", "temperature_avg_K", "temperature_std_K"]:
            self.assertAlmostEqual(refined[field], original[field], places=8, msg=field)
        self.assertEqual(refined["node_count"], 4)
        self.assertEqual(refined["element_count"], 3)
        self.assertNotAlmostEqual(sum(n.temperature for n in nodes.values()) / 4, refined["temperature_avg_K"])

    def test_nodal_mode_reports_provisional_statistics_without_area(self):
        nodes = self.triangle_nodes()
        panels = self.panels({"id": "a", "name": "Skin A", "selector": {"type": "node_ids", "ids": [1, 2, 3]}})
        result = summarize(nodes, [], panels, assign_panels(nodes, [], panels, "nodal"), "nodal")[0]
        self.assertAlmostEqual(result["temperature_avg_K"], 301.0)
        self.assertAlmostEqual(result["temperature_std_K"], math.sqrt(2.0))
        self.assertIn(result["area_m2"], (None, ""))
        self.assertEqual(result["node_count"], 3)
        self.assertEqual(result["element_count"], 0)
        self.assertIn("nodal", result["averaging_method"])

    def test_uniform_field_has_zero_surface_standard_deviation(self):
        nodes = {n: Node(n, node.x, node.y, node.z, node.pressure, 1843.398156) for n, node in self.triangle_nodes().items()}
        result = self.area_summary(nodes, [Triangle(20, (1, 2, 3))])
        self.assertAlmostEqual(result["temperature_avg_K"], 1843.398156)
        self.assertAlmostEqual(result["temperature_std_K"], 0.0, places=10)


class CommandLineTests(TemporaryInputs):
    def inputs(self):
        nodes = self.write("nodes.csv", HEADER + "1,0,0,0,10,300\n2,2,0,0,20,300\n3,0,1,0,30,303\n")
        panels = self.write("map.json", json.dumps({"panels": [{"id": "skin", "name": "Skin", "solidworks_component": "skin_part", "solidworks_region": "outside", "selector": {"type": "node_ids", "ids": [1, 2, 3]}}]}))
        return nodes, panels

    def run_cli(self, *args):
        return subprocess.run([sys.executable, "-m", "paneltemps", *map(str, args)], cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True)

    def test_missing_elements_cannot_silently_fall_back_to_nodal(self):
        nodes, panels = self.inputs()
        result = self.run_cli(nodes, "--panels", panels, "--output-dir", self.root / "output")
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_explicit_nodal_cli_produces_csv_and_preserves_input(self):
        nodes, panels = self.inputs()
        before = nodes.read_bytes()
        output = self.root / "output"
        result = self.run_cli(nodes, "--panels", panels, "--method", "nodal", "--output-dir", output)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(nodes.read_bytes(), before)
        tables = []
        for path in output.glob("*.csv"):
            with path.open(newline="", encoding="utf-8") as stream:
                reader = csv.DictReader(stream)
                if "temperature_avg_K" in (reader.fieldnames or []):
                    tables.extend(reader)
        self.assertEqual(len(tables), 1)
        self.assertEqual(tables[0]["panel_id"], "skin")
        self.assertAlmostEqual(float(tables[0]["temperature_avg_K"]), 301.0)
        self.assertTrue(list(output.glob("*.html")), "Expected an assignment/temperature visualization")

    def test_output_directory_is_never_overwritten(self):
        nodes, panels = self.inputs()
        output = self.root / "existing"
        output.mkdir()
        sentinel = output / "my-results.txt"
        sentinel.write_text("keep this", encoding="utf-8")
        result = self.run_cli(nodes, "--panels", panels, "--method", "nodal", "--output-dir", output)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep this")
        self.assertEqual(list(output.iterdir()), [sentinel])

    def test_cli_overlap_fails_without_creating_output(self):
        nodes, _ = self.inputs()
        panels = self.write("overlap.json", json.dumps({"panels": [
            {"id": "a", "name": "A", "selector": {"type": "node_ids", "ids": [1, 2]}},
            {"id": "b", "name": "B", "selector": {"type": "node_ids", "ids": [2, 3]}},
        ]}))
        output = self.root / "invalid"
        result = self.run_cli(nodes, "--panels", panels, "--method", "nodal", "--output-dir", output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("overlap", result.stderr.lower())
        self.assertFalse(output.exists())

    def test_area_cli_records_mapping_units_and_exact_preserved_sources(self):
        nodes, _ = self.inputs()
        elements = self.write("elements.csv", ELEMENT_HEADER + "20,1,2,3,outer\n")
        panels = self.write("area_map.json", json.dumps({"panels": [
            {"id": "a", "name": "Exterior", "solidworks_component": "skin_part", "solidworks_region": "outer_set", "selector": {"type": "mesh_panel", "value": "outer"}},
        ]}))
        output = self.root / "area_output"
        result = self.run_cli(nodes, "--panels", panels, "--elements", elements, "--output-dir", output,
                              "--case-label", "synthetic_test", "--wall-temperature-kind", "structural", "--warn-range-K", "2")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with (output / "panel_temperatures.csv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(float(rows[0]["area_m2"]), 1.0)
        self.assertAlmostEqual(float(rows[0]["temperature_avg_K"]), 301.0)
        self.assertAlmostEqual(float(rows[0]["temperature_std_K"]), math.sqrt(0.5))
        with (output / "solidworks_mapping.csv").open(newline="", encoding="utf-8") as stream:
            mapping = list(csv.DictReader(stream))
        self.assertEqual(mapping[0]["solidworks_component"], "skin_part")
        self.assertEqual(mapping[0]["solidworks_region"], "outer_set")
        report = json.loads((output / "run_report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["case_label"], "synthetic_test")
        self.assertEqual(report["input_coordinate_unit"], "m")
        self.assertEqual(report["input_temperature_unit"], "K")
        self.assertEqual(report["wall_temperature_kind"], "structural")
        self.assertFalse(report["cad_mapping_verified"])
        self.assertFalse(report["source_units_verified"])
        self.assertEqual(report["panels_flagged_for_range_review"], ["a"])
        originals = {"boundary_conditions": nodes, "panel_mapping": panels, "connectivity": elements}
        for item in report["preserved_inputs"]:
            original = originals[item["role"]].read_bytes()
            self.assertEqual((output / item["preserved_path"]).read_bytes(), original)
            self.assertEqual(item["sha256"], hashlib.sha256(original).hexdigest())

    def test_partial_coverage_is_explicit_and_reports_omitted_ids(self):
        nodes, _ = self.inputs()
        panels = self.write("partial.json", json.dumps({"panels": [
            {"id": "a", "name": "Partial", "selector": {"type": "node_ids", "ids": [1, 2]}},
        ]}))
        output = self.root / "partial"
        result = self.run_cli(nodes, "--panels", panels, "--method", "nodal", "--allow-unmapped", "--output-dir", output)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads((output / "run_report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["unmapped_entity_ids"], [3])
        assignments = json.loads((output / "panel_assignments.json").read_text(encoding="utf-8"))
        self.assertEqual(assignments["unmapped_ids"], [3])
        self.assertIn("Partial coverage", result.stderr)


if __name__ == "__main__":
    unittest.main()
