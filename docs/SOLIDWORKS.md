# Apply panel temperatures in SOLIDWORKS Simulation

This is a manual first-pass workflow. The supplied handoff contains no CAD model,
surface mesh, material definitions, or temperature-source metadata. Consequently,
neither exterior-skin coverage nor correspondence with actual CAD faces has been
verified. Demonstration outputs are synthetic and must not be applied to a real
assembly.

## Establish what the temperatures mean

Confirm the coordinate frame, coordinate units, temperature units, case identifier,
and temperature definition with the source-data owner before applying any loads:

| Source quantity | Interpretation for this workflow |
| --- | --- |
| Prescribed CFD wall temperature | A solver input; it does not independently predict the skin temperature. Imposing it in Simulation is a stipulated temperature scenario. |
| Adiabatic-wall temperature | A reference temperature for aerodynamic heat transfer under an adiabatic-wall condition. It is not generally the temperature of a conducting exterior skin. Imposing it directly is an additional assumption, not a coupled heating prediction. |
| Calculated structural surface temperature | A candidate prescribed surface temperature for the matching location and case, subject to the assumptions and units of the calculation that produced it. |
| Unknown wall-temperature definition | Keep the interpretation unresolved in the analysis record; do not call the resulting calculation a validated heating prediction. |

Prescribing an area-average temperature does not preserve the original heat flux
or total heat input. Heat flux depends on temperature gradients, materials, contacts,
and other boundary conditions. A temperature boundary can supply or remove whatever
heat the mathematical model requires to maintain its specified temperature.

## Map panels to exterior faces

Use the exported mapping CSV as a reviewed association table. Its identifying
columns are:

```csv
panel_id,panel_name,solidworks_component,solidworks_region,selection_notes
```

`solidworks_component` should identify the assembly component instance, not just
the part filename. `solidworks_region` should identify a saved selection set or an
explicitly documented set of exterior faces. Use `selection_notes` to record the
CAD revision, coordinate transform, selection procedure, and any unresolved match.
Blank regions remain unmapped; an exported row alone does not establish a match.
Several CAD faces may belong to one panel, but a physical face must not receive
conflicting temperatures from overlapping panels.

1. Open the correct assembly revision and resolve the intended component instances.
   Compare overall dimensions and reference landmarks with the mesh visualization.
   Check units, handedness, axis orientation, origin, and assembly placement. A
   temperature conversion does not fix a coordinate-frame mismatch.
2. Compare the panel-assignment visualization with the CAD. Start with a small,
   easily identified region and check both its location and extent. Coincident
   shared-edge nodes can belong to adjacent panels; overlapping surface areas
   cannot be accepted merely because their points look close.
3. Select the physical outer face or faces for each panel and save a named selection
   where the installed version supports it. Record the full component-instance
   path and selection name in the mapping CSV. A part face in a repeated component
   is not sufficient to identify an assembly face.
4. Record every unresolved panel and every exterior CAD region absent from the
   temperature dataset. Do not infer full-vehicle coverage from node counts.
   Compare areas and boundaries when both mesh and CAD areas are available.
5. Recheck the mapping after CAD edits. Simple face indices, visual ordering, or
   nearest-point matches are not reliable persistent identifiers across revisions.

If coordinates are the only available geometry, obtain surface connectivity or
CAD geometry for an area average. Arithmetic nodal averages are provisional and
can be biased by nonuniform node density. A scatter plot alone does not establish
surface coverage, face identity, or panel area.

## Check the physical model

Inspect whether the assembly includes the exterior skin, its thickness, and its
thermal connection to the structure. If the skin is missing, placing the exterior
temperature directly on internal members bypasses the skin's thermal resistance
and, in a transient study, its heat capacity. The effect on a particular internal
temperature depends on the whole thermal network; it is not a generally valid
replacement. Add an appropriate skin representation or a justified equivalent
thermal model before drawing conclusions about internal heating. Verify how any
shell representation handles through-thickness conduction and face temperatures.

Assign appropriate material properties to skin, structural members, and internal
components. Thermal conductivity is needed for conduction. Density and specific
heat are also needed for transient heat storage. Use temperature-dependent values
where the study requires them and suitable source data are available.

Review contact definitions rather than accepting mechanical assembly contact as
proof of a thermal path. Model thermal contact resistance or conductance as
appropriate. Gaps and cavities do not conduct heat like bonded solids; radiation
or fluid convection can be important between separated components. A conduction-only
model cannot predict those missing paths. No new external convection coefficients,
aerodynamic heat fluxes, or trajectory heating are calculated by this tool.

## Create and run the study manually

Menu wording and available features vary with the installed SOLIDWORKS release
and Simulation license. Use the installed release's thermal-study help for its
exact dialogs.

1. Enable the SOLIDWORKS Simulation add-in and create a thermal study of the
   intended assembly configuration. Choose steady state or transient deliberately.
2. Define materials and thermal contacts, and include any other justified heat
   sources or thermal boundaries. Do not invent values to complete the model.
3. For each reviewed mapping row, select its physical exterior face set. Add a
   prescribed **Temperature** thermal boundary using `temperature_avg_K` from the
   panel-results CSV. Name the load with the panel ID so it can be audited.
4. Select Kelvin in the temperature dialog when available. If the UI is set to
   Celsius, use `T_C = T_K - 273.15`. Verify the displayed unit and the value actually
   stored; do not paste a Kelvin number into a Celsius field. Do not assign these
   exterior temperatures directly to the payload.
5. Inspect every load's selected faces. Check that all intended exterior regions
   are represented and no face receives conflicting prescribed values. Keep the
   temperature CSV, reviewed mapping, source-case identity, and CAD revision
   together in the analysis record.
6. For transient studies, supply justified initial temperatures, an exposure
   duration, and time-stepping controls. Apply the exterior values as constant over
   that duration. The input file contains no trajectory history; do not create a
   time-varying temperature curve from it.
7. Mesh the thermal model. Examine convergence at relevant internal locations and
   thermal contacts, refining the spatial mesh and transient time steps as needed.
   Run the study and inspect solver warnings, energy balance, and heat flow at the
   prescribed-temperature boundaries as well as temperatures.
8. Probe the payload and other critical components. For transient studies, review
   their temperature histories over the supplied duration, not just the final
   frame. Record the assumptions and unresolved mappings with the results.

| Study | What it answers | Additional inputs |
| --- | --- | --- |
| Fixed-temperature steady state | The equilibrium distribution under the specified constant boundaries; it does not describe time to heat up. | Conductivity, contacts, and a well-posed set of boundaries and heat sources. |
| Transient conduction with constant exterior temperatures | How the modeled structure evolves from its specified initial state during a chosen exposure to these fixed boundaries. | Initial temperatures, density, specific heat, duration, and adequate time resolution. |
| Transient aerodynamic heating | A changing external heating environment coupled to thermal response. This is outside this implementation. | Time-resolved aerodynamic and thermal inputs or a separate coupled model. |

One panel average suppresses local hot regions and gradients. Review each panel's
minimum, maximum, standard deviation, area, and assignment plot before accepting
that simplification. Subdivide only where the source field, physical geometry,
and intended decision justify it; do not automatically generate many small panels.
Adjacent constant-temperature regions can introduce abrupt boundary jumps and
locally large gradients near their shared edge.

## Automation feasibility and limits

SOLIDWORKS exposes Windows automation interfaces, and the separately documented
Simulation API is the relevant route for thermal-study operations. General model
selection through the main SOLIDWORKS API alone is not a verified thermal-load
importer. Consult the vendor's [SOLIDWORKS API Help](https://help.solidworks.com/)
for the installed release, including its SOLIDWORKS Simulation API reference and
thermal examples. Check the installed product's API libraries, licensing, and
supported interfaces before implementing a macro.

An automated importer is plausible once the following are available:

- A Windows machine with the intended SOLIDWORKS release and the required
  Simulation license/add-in.
- A known assembly revision and a reviewed association between each panel ID and
  its component instance and face set.
- A tested selection strategy using saved selections or persistent references
  where supported, with explicit checks when references no longer resolve.
- A version-tested way to create or update thermal loads, set temperature units,
  and read back both selected entities and values.

That environment and CAD mapping are not available in this Linux cloud workspace.
No SOLIDWORKS API calls or generated macro have been executed or verified here.
The manual method above is the supported deliverable for this first implementation.
A later importer should reject unmapped or duplicate selections, offer a dry-run
report, and validate one known face and temperature by reading them back before
applying the rest. It must not silently choose the nearest CAD face or rely on
unstable face numbers.
