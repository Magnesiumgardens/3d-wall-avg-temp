# Surface panel temperature averaging

A standalone Python post-processor for six-column nodal boundary-condition data.
It exports panel temperatures, an auditable panel/CAD association table, coverage
reports, and a self-contained interactive 3D visualization. It does not modify or
run aerodynamic or thermal solvers.

**Current data status:** only the written handoff was supplied. There is no actual
vehicle dataset, mesh connectivity, or CAD model in this checkout. Everything in
`examples/` and the generated `outputs/synthetic_*` directories is a synthetic
software demonstration. No real panel definitions, units, CAD matches, exterior
skin, component coverage, or wall-temperature interpretation have been verified.

## Run

Python 3.10+ is required (validated here on Python 3.12).
Use the existing checkout; cloud tasks already have an isolated workspace and do
not need a separate Git worktree.

```bash
cd /workspace/3d-wall-avg-temp
python3 -m venv /workspace/.venvs/paneltemps
/workspace/.venvs/paneltemps/bin/python -m pip install -r requirements.txt
/workspace/.venvs/paneltemps/bin/python -m unittest discover -s tests -v
```

Run the surface-area demonstration:

```bash
/workspace/.venvs/paneltemps/bin/python -m paneltemps examples/synthetic_nodes.csv \
  --elements examples/synthetic_triangles.csv \
  --panels examples/panels_area.json \
  --coordinate-unit m --temperature-unit K \
  --method area --case-label SYNTHETIC_DEMO \
  --output-dir outputs/synthetic_area
```

For a point-only dataset, choose the provisional method explicitly:

```bash
/workspace/.venvs/paneltemps/bin/python -m paneltemps examples/synthetic_nodes.csv \
  --panels examples/panels_nodal.json --method nodal \
  --case-label SYNTHETIC_NODAL_DEMO --output-dir outputs/synthetic_nodal
```

Use a **new or empty output directory** for each run; previous outputs are never
overwritten. `python -m paneltemps --help` lists all options. No background
services, credentials, or SOLIDWORKS installation are needed for post-processing.
Only the plotting dependency needs installation; exported HTML works offline.

## Inputs and units

Nodal data may be comma-separated or whitespace-separated, with this optional
header, blank lines, and whole-line `#` comments:

```csv
nodenumber,x-coordinate,y-coordinate,z-coordinate,pressure,wall-temperature
1,1.642018079E+00,-6.750000000E-02,-2.448393844E-01,1.017518533E+04,1.843398156E+03
```

Exactly six fields are required in this order. Scientific notation including
Fortran `D` exponents is accepted. IDs must be unique positive integers. Invalid
numbers, nonfinite values, and temperatures below absolute zero are rejected.
Pressure is preserved but not used in the temperature calculation; its unit need
not be inferred for this workflow.

Coordinates default to **meters**, and temperatures to **Kelvin**. These are
assumptions until checked against source metadata. Override with
`--coordinate-unit m|mm|cm|in|ft` and `--temperature-unit K|C|F`. Calculations,
visualizations, areas, and exports always use meters, square meters, and Kelvin.
Use `--case-label` to record the source case and `--wall-temperature-kind
unknown|prescribed|adiabatic|structural` to record its declared interpretation.
This option records metadata; it cannot establish the physical meaning of a CFD
field. Default `unknown` remains an explicit unresolved warning.

Area mode requires **linear surface triangles**, supplied as:

```csv
element_id,node_1,node_2,node_3,panel_id
10,1,2,3,skin_A
```

`panel_id` is optional; the first four fields are required. Connectivity may also
be whitespace-separated and headerless. Missing nodes, duplicate element IDs,
repeated-node/zero-area triangles, and repeated triangles with the same node IDs
in any order are rejected. This checks the supplied connectivity, not physical
mesh validity: coincident faces with different node IDs, self-intersections,
holes, interior faces, and nonmanifold surfaces require source/CAD inspection.
No surface is reconstructed from an unconnected point cloud.

Other source formats, quadrilaterals, high-order elements, volumes, CAD files,
and proprietary solver files are not parsed directly. Obtain/export linear
surface connectivity with preserved node IDs and panel labels. Any triangulation
of other element types must use a justified interpolation convention before
using this tool; splitting a quad arbitrarily can change the represented field.

## Panel definitions

The JSON mapping explicitly defines the intended regions. Preserve trustworthy
surface labels from the source when available; do not guess parts from point
positions or claim that an unlabeled point cloud identifies CAD faces.

```json
{
  "panels": [
    {
      "id": "skin_A",
      "name": "Reviewed exterior region A",
      "solidworks_component": "assembly/component-instance",
      "solidworks_region": "reviewed-exterior-face-set",
      "selection_notes": "Record CAD revision and coordinate alignment here",
      "selector": {"type": "mesh_panel", "value": "skin_A"}
    }
  ]
}
```

Each panel has one selector:

| Selector | Fields | Use |
| --- | --- | --- |
| `mesh_panel` | `value`: label string from connectivity | Area mode; use existing source labels. |
| `element_ids` | `ids`: positive integer array | Area mode; explicit reviewed element selection. |
| `node_ids` | `ids`: positive integer array | Provisional nodal mode. |
| `bounds` | `min`: XYZ array, `max`: XYZ array | Either mode; axis-aligned limits **always in meters**, inclusive. |

Bounds select a triangle only when **all three vertices** are inside the box.
They do not clip triangles or split them across panels. Crossing elements remain
unmapped and fail the coverage check by default. Bounds are a coarse configurable
selection aid, not a CAD face-matching algorithm. Adjacent nodal bounds sharing a
boundary can overlap and will fail; use explicit IDs for ambiguous cases.

Panel IDs and names must be nonempty; IDs are unique. Empty panels, unknown IDs,
and overlapping entity selections are errors. Area mode permits shared edge
nodes between different panels, since the triangles rather than the nodes must
be disjoint. All supplied triangles (or all nodes in nodal mode) must be assigned
unless `--allow-unmapped` is explicitly given. Omitted entities are always listed
and shown gray. Nodes without connectivity are reported separately in area mode;
their presence cannot establish a surface area. Coverage of a supplied file does
not establish coverage of the entire vehicle or of all exterior CAD faces.

CAD fields are optional until reviewed. Empty fields are reported as unresolved;
even filled fields remain unverified descriptions until checked in SOLIDWORKS.
No geometric registration or automatic nearest-face matching is performed.

## Calculation

For triangle `e`, compute its 3D area `A_e` from the cross product and interpolate
temperature linearly from its three nodes. The exact area integral of that
interpolated field gives:

```text
T_mean,e = (T1 + T2 + T3) / 3
T_mean,p = sum(A_e * T_mean,e) / sum(A_e)
```

The reported standard deviation is the **population surface standard deviation**,
using the exact second moment of the same linear field. With
`d_i = T_i - T_mean,p`, each triangle contributes:

```text
mean_e((T - T_mean,p)^2) = (sum(d_i^2) + (sum(d_i))^2) / 12
std_p = sqrt(sum(A_e * mean_e((T - T_mean,p)^2)) / sum(A_e))
```

Minimum and maximum are the extrema among contributing vertices, which are also
the extrema of the piecewise-linear field. Node counts are unique within each
panel; shared nodes may contribute to more than one area panel. Element count
counts contributing triangles. Refining the same linear field nonuniformly leaves
the area, mean, and standard deviation unchanged; the tests verify this.

Nodal mode instead reports an arithmetic mean and population **nodal** standard
deviation. Its method is labeled `provisional_nodal_arithmetic`, area is blank,
and element count is zero. Nonuniform point density can bias these statistics;
they are not accurate area averages. Area mode never silently falls back to this
method.

Review temperature ranges and gradients before applying one average per panel.
Optional `--warn-range-K 100` flags ranges of 100 K or more for review (100 is an
example threshold, not a universal acceptance criterion). Choose the threshold
for the intended analysis. No automatic subdivision is performed.

## Outputs

| File | Contents |
| --- | --- |
| `panel_temperatures.csv` | One row per panel with the requested fields below, plus `element_count` and `averaging_method`. |
| `solidworks_mapping.csv` | Panel names, component-instance paths, face-set descriptions, selection notes. |
| `panel_visualization.html` | Offline interactive 3D assignment and average-temperature views; actual triangles in area mode, points in nodal mode. |
| `panel_assignments.json` | Exact contributing entity IDs, omitted IDs, nodes without elements. |
| `run_report.json` | Units, case interpretation, coverage, CAD mapping gaps, temperature ranges, warnings, input hashes. |
| `inputs/` | Byte-for-byte copies of nodal data, connectivity when supplied, and panel configuration. Originals are read-only inputs. |

```csv
panel_id,panel_name,area_m2,temperature_avg_K,temperature_min_K,temperature_max_K,temperature_std_K,node_count,element_count,averaging_method
```

Open the HTML file locally in a WebGL-capable browser. Both views show the same
geometry; the temperature view deliberately displays the exported averages,
not the original detailed field. Check the original inputs as well when
assessing localized hot spots. Very large meshes can produce large HTML files.

For the synthetic two-panel example, each rectangle has area 2 m². Its means are
325 K (base) and 315 K (side), with ranges 300–350 K and 300–330 K respectively.
These are known demonstration values, not vehicle results.

## SOLIDWORKS and next inputs

Follow [the manual thermal-study and CAD-mapping guide](docs/SOLIDWORKS.md).
It explains prescribed exterior temperatures, material and contact inputs,
missing skin resistance/capacity, steady state versus constant-boundary transient
conduction, and why Simulation API automation remains unverified here.

For actual results, supply the complete boundary-condition file, any surface
connectivity and named regions, confirmed units and wall-temperature definition,
and the CAD model/revision or a reviewed face-association table. Without these,
no actual panel-temperature table, CAD assignment, or internal thermal result
can be validated. The tool supplies boundary values; SOLIDWORKS performs the
separate thermal solve. No convection coefficients, new aerodynamic heat fluxes,
temperature history, or trajectory simulation are generated.
