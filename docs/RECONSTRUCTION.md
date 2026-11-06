# STEP-based temperature reconstruction

This optional workflow reconstructs a temperature field on the supplied CAD
surfaces when original CFD surface connectivity is unavailable. Its numerical
checks do not establish the physical meaning of source temperatures or approve
SOLIDWORKS face selections. Outputs are labeled as reconstructed CAD averages.
If the original meshing/solver software can export the matching surface elements,
coordinates and region labels, prefer those original elements. A Pointwise `.pw`
project is not itself an ANSYS solver mesh. CGNS, solver-native or surface-only
exports can be inspected for conversion to the original triangle workflow;
arbitrary quad/polygon splitting requires its own interpolation convention.

## Windows: selected MK3 cases

Use Python **3.12**, which has wheels for the pinned CAD dependency. The original
point-only workflow can keep using its existing environment. Put these three
original files in `data/` under your local repository:

```text
data/MK3_watertight.STEP
data/M08_ALT30_A00_B0_base_Iso_5spe_2_redo
data/M08_ALT30_A00_B0_base_Adiab_5spe_9
```

From the repository folder in PowerShell:

```powershell
git pull origin main
powershell -NoProfile -File .\scripts\run_mk3.ps1
```

The helper creates `.venv-cad` using Python 3.12, installs `requirements-cad.txt`,
and runs both cases into a new timestamped `outputs/mk3_*` folder. No activation
script is needed. If your PowerShell policy blocks the helper, use the manual
commands below. It expects the exact
filenames above. CAD packages download on the first installation and can be large.

For custom paths/filenames or explicit manual steps:

```powershell
py -3.12 -m venv .venv-cad
.\.venv-cad\Scripts\python.exe -m pip install -r requirements-cad.txt
$runOutput = "outputs/mk3_$(Get-Date -Format 'yyyyMMdd_HHmmss')"
.\.venv-cad\Scripts\python.exe -m paneltemps reconstruct --step ".\data\MK3_watertight.STEP" --bc ".\data\M08_ALT30_A00_B0_base_Iso_5spe_2_redo" ".\data\M08_ALT30_A00_B0_base_Adiab_5spe_9" --transform examples/mk3_alignment.json --coordinate-unit m --temperature-unit K --temperature-kinds prescribed adiabatic --output-dir $runOutput
Invoke-Item $runOutput
```

The two case folders each contain `panel_temperatures.csv`,
`solidworks_mapping.csv`, `panel_visualization.html`, the original BC file, and
the original node-ID associations. Open the HTML in a WebGL-capable browser. It
works offline. Gray temperature faces have unresolved results or were explicitly
excluded. The common output directory contains:

- `reconstruction_report.json`: thresholds, warnings, individual face coverage,
  checks for both cases, and any unresolved panels.
- `source_mapping.csv`: every unique BC location, transformed CAD position,
  candidate face memberships, distance, and unmatched status.
- `face_inventory.json`: exact face areas, centers, bounding boxes, and CAD hash.
- `face_groups_template.json`: editable physical panel grouping and CAD association
  fields tied to the exact STEP revision.
- `inputs/`: copies of geometry, transform and grouping configuration.

Exit code **0** means the selected numerical checks passed; **3** means diagnostic
outputs were written but some panel results or source points remain unresolved;
**2** means an input, dependency, or execution error. With exit 3, valid panel rows
remain usable for review, while failed averages stay blank. `--allow-partial`
changes exit 3 to 0 without changing any status, thresholds, or blank values.

## Alignment and panel grouping

`examples/mk3_alignment.json` holds the fitted BC-to-CAD rigid transform, in meters.
It is guarded by hashes of both the STEP file and the normalized unique source
coordinate set. Row reordering and consistent coordinate duplicates do not change
that coordinate hash. Different units, CAD revisions, or BC point sets fail the
guard rather than silently reusing this fit. For other geometries supply your own
validated transform with `rotation_BC_to_CAD` and `translation_m`; identity is an
explicit 3x3 identity rotation and zero translation. Hash fields are optional for
general transforms but are included for the supplied MK3 case.

By default every CAD face is one **candidate** panel. This avoids inventing
physical panel boundaries. It is not a claim that the CAD face decomposition
matches manufacturing panels or reviewed SOLIDWORKS named regions. To group
several CAD faces into one physical panel:

1. Inspect face IDs and positions in the visualization and inventory.
2. Copy `face_groups_template.json` outside the output directory to edit it.
3. Combine the desired `face_indices` under a panel ID/name and fill in the CAD
   component-instance, reviewed region, and selection notes.
4. Account for every face exactly once: group it or list it in
   `excluded_face_indices`. Unlisted, duplicate, or overlapping faces are errors.
5. Rerun with `--face-groups path/to/reviewed_groups.json` and a new output folder.

The file's `cad_sha256` must match the STEP. Face indices are inspection indices,
not persistent SOLIDWORKS selection IDs; rerun/review after CAD or CAD-reader
changes. Group results combine surface integrals, not arithmetic means of faces.
No SOLIDWORKS API loads are created.

## Reconstruction method

1. Parse named source columns in the confirmed input units. Prefer
   `wall-temperature` when present, otherwise `temperature`. Pressure and heat
   flux are never substituted for temperature.
2. Merge exact-coordinate duplicates only when temperatures agree within
   `--duplicate-tolerance-K` (default 0.000001 K). Cases are aligned by sorted
   coordinate triples, independently of their node IDs. Batch cases must have the
   same unique coordinates; process other point sets separately.
3. Read the STEP's declared length units and transfer all geometry into meters.
   Apply the explicit rigid transform; scale, shear and reflections are rejected.
4. Check every source position against finite CAD faces using bounding-box lower
   bounds and exact surface distances. Default mapping tolerance is 0.0005 m.
   Use the nearest finite face. Share a boundary sample with another face only
   when both faces share a topological edge and the sample is within tolerance
   of that edge. Nearby unrelated surfaces do not exchange samples. Shared-edge
   assignment tolerance itself can blur a small boundary strip and needs review.
   Equidistant competing surfaces without a local common edge are flagged as
   ambiguous, rather than chosen by face order.
5. Project each face's samples into its own surface UV coordinates and build a
   Delaunay piecewise-linear interpolant there. The tangent metric at the sample
   center scales UV distances into physical lengths for triangulation; this is
   exact for a planar surface and a local approximation on curved surfaces.
   Distinct 3D points can collapse to one trimmed-edge location. The closest
   source point supplies that projected location; equally close conflicting
   values are excluded as ambiguous, with the collision groups reported. Exact
   source-coordinate duplicates still require agreeing temperatures. Fewer than
   three distinct noncollinear retained face samples leave the face unsupported.
6. Reject interpolation triangles whose edge/centroid trim checks cross outside
   the finite face. Target quadrature points use only their own face's data.
   Every source vertex of an interpolation stencil must lie within
   `--max-support-m` (default 0.04 m) of the target point.
7. A small extension outside the sampled UV convex hull can use the nearest
   same-face source value within `--extrapolation-m` (default 0.003 m), provided
   the checked UV path stays inside the trim and its surface length also meets
   that limit. This is constant boundary extrapolation, explicitly measured and
   limited; it does not fill holes rejected by an interior stencil.
8. Integrate on CAD UV triangles using three-point Gaussian quadrature and the
   exact surface Jacobian. Longest-edge subdivision limits the physical triangle
   edges. Repeat at two spatial resolutions: default maximum edges 0.02 m and
   0.01 m. Compare integrated area with the exact BRep surface area, coverage,
   mean, and surface standard deviation. A face failing mean/std or supported-area
   convergence receives up to two further edge halvings by default, comparing
   each new result with its predecessor. `--adaptive-levels 0` disables additional
   levels; values 1–3 set the limit. Per-face final resolutions are reported.

The default acceptance checks are at least 98% supported quadrature area, at most
2% extrapolated quadrature area, at most 1% error relative to exact CAD area,
at most 1 K change in both mean and standard deviation, and at most 0.5% of exact
area changing support between resolutions. These are numerical review thresholds,
not universal engineering acceptance criteria. Defaults can be set explicitly via
`--min-coverage`, `--max-extrapolated-fraction`, `--max-area-error`,
`--convergence-K`, and `--coverage-change`. Do not relax them simply to make rows
pass; inspect the reason and obtain better coverage or refine the reconstruction.

## Reading the CSV and limitations

`area_m2` is the selected exact CAD area. `integrated_area_m2` is quadrature area;
`supported_area_m2` is the portion with accepted interpolation support. The
reported `temperature_avg_K` and `temperature_std_K` integrate **supported area**
and remain blank when checks fail. Even passing defaults can permit up to 2%
unsupported area: this omission must be assessed, not implicitly assigned the
same mean. The `supported_temperature_*` fields retain diagnostic conditional
statistics on unresolved rows and must not be mistaken for accepted full-panel
results. Minimum and maximum temperatures are extrema among unique source samples
assigned to the selected faces, not proof of extrema in the unsampled physical
field. Node count counts those assigned unique coordinates, including samples
that may be discarded by projected-collision or interpolation-support checks;
it is not a count of nonzero quadrature weights. Interpolation sample counts and
projection-collision exclusions are reported separately for each face. Original
IDs and source-coordinate duplicates are preserved.

Refinement checks numerical integration of the chosen interpolant, not the error
between that interpolant and the unknown physical field. Trim tests sample
stencil paths and cannot guarantee detection of arbitrarily narrow trim features.
Periodic seams, holes, thin regions, incomplete data, nearby surfaces, and sharp
temperature gradients deserve closer inspection. No temperature is borrowed from
another face to fill missing data. Nonuniform source density no longer directly
sets quadrature weights, but still affects interpolation quality. For higher
confidence compare more integration resolutions, support limits and source
interpolation choices or obtain original connectivity.

Keep the isothermal and adiabatic cases separate. A constant 300 K isothermal wall
does not become a predicted hot skin; adiabatic wall temperatures also do not
independently predict conducting-skin temperatures. Neither area averaging nor
this reconstruction conserves the original heat flux. Skin thickness, material
properties, contacts, initial conditions/duration, and internal heat-transfer
paths remain separate modeling requirements described in
[the SOLIDWORKS guide](SOLIDWORKS.md).
