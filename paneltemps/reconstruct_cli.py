"""CLI for auditable STEP-based reconstruction of one or more BC cases."""

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from importlib.metadata import version
import platform
import shutil
import sys

import numpy as np

from .cad import coordinate_hash, load_cad, load_transform, sha256
from .io import COORDINATE_SCALES, read_nodes, _rows, _header
from .reconstruction import FaceInterpolator, integrate, map_samples, unique_samples
from . import __version__


def _positive(value):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError('must be finite and positive')
    return value


def _fraction(value):
    value = float(value)
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise argparse.ArgumentTypeError('must be between 0 and 1')
    return value


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def panels_for_faces(path, faces, cad_hash):
    if path is None:
        return [{'id': f'face_{f.index:03}', 'name': f'CAD face {f.index} (unreviewed)',
                 'face_indices': [f.index], 'solidworks_component': '', 'solidworks_region': '',
                 'selection_notes': 'Candidate CAD face; inspection index is not a persistent SOLIDWORKS ID.'}
                for f in faces], []
    config = json.loads(path.read_text(encoding='utf-8'))
    if config.get('cad_sha256') != cad_hash:
        raise ValueError('face-group file must contain the exact cad_sha256 from face_inventory.json')
    panels, seen_ids, seen_faces = config.get('panels', []), set(), set()
    if not panels:
        raise ValueError('face-group file contains no panels')
    for panel in panels:
        if not isinstance(panel.get('id'), str) or not panel['id'].strip() or panel['id'] in seen_ids:
            raise ValueError('face-group panel IDs must be unique nonempty strings')
        if not isinstance(panel.get('name'), str) or not panel['name'].strip():
            raise ValueError('each face-group panel requires a name')
        ids = panel.get('face_indices')
        if not isinstance(ids, list) or not ids or any(type(i) is not int or i < 1 or i > len(faces) for i in ids):
            raise ValueError('face_indices must be a nonempty list of valid CAD face indices')
        if len(set(ids)) != len(ids) or set(ids) & seen_faces:
            raise ValueError('overlapping or duplicate CAD faces in panel groups')
        for key in ('solidworks_component', 'solidworks_region', 'selection_notes'):
            if not isinstance(panel.setdefault(key, ''), str):
                raise ValueError('CAD association fields must be strings')
        seen_ids.add(panel['id'])
        seen_faces.update(ids)
    excluded = config.get('excluded_face_indices', [])
    if not isinstance(excluded, list) or any(type(i) is not int or i < 1 or i > len(faces) for i in excluded):
        raise ValueError('invalid excluded_face_indices')
    if len(excluded) != len(set(excluded)) or set(excluded) & seen_faces:
        raise ValueError('duplicate or simultaneously selected/excluded CAD faces')
    if seen_faces | set(excluded) != {f.index for f in faces}:
        raise ValueError('every CAD face must be grouped or explicitly excluded')
    return panels, excluded


def _panel_results(panels, faces, integrations, temperatures, mappings, thresholds):
    results = []
    for panel in panels:
        face_ids = panel['face_indices']
        exact_area = sum(faces[i-1].area for i in face_ids)
        fine_area = sum(integrations[i][1]['area'] for i in face_ids)
        supported_area = sum(integrations[i][1]['supported_area'] for i in face_ids)
        coverage = supported_area/fine_area if fine_area else 0
        extrapolated = sum(integrations[i][1]['extrapolated_area'] for i in face_ids)
        node_indices = sorted({index for i in face_ids for index, _ in mappings[i]})
        for case_index in range(len(temperatures)):
            def moments(level):
                parts = [integrations[i][level] for i in face_ids]
                parts = [p for p in parts if p['cases'][case_index] is not None]
                total = sum(p['supported_area'] for p in parts)
                if total == 0:
                    return None
                mean = sum(p['supported_area']*p['cases'][case_index]['mean'] for p in parts)/total
                variance = sum(p['supported_area']*(p['cases'][case_index]['std']**2+
                               (p['cases'][case_index]['mean']-mean)**2) for p in parts)/total
                return mean, float(np.sqrt(max(0, variance)))
            coarse, fine = moments(0), moments(1)
            mean_change = abs(coarse[0]-fine[0]) if coarse and fine else None
            std_change = abs(coarse[1]-fine[1]) if coarse and fine else None
            area_error = abs(fine_area-exact_area)/exact_area
            reasons = []
            if fine is None:
                reasons.append('no supported surface area')
            if coverage < thresholds.min_coverage:
                reasons.append('insufficient area coverage')
            if fine_area and extrapolated/fine_area > thresholds.max_extrapolated_fraction:
                reasons.append('too much boundary extrapolation')
            if area_error > thresholds.max_area_error:
                reasons.append('CAD integration area mismatch')
            if mean_change is None or mean_change > thresholds.convergence_K:
                reasons.append('mean refinement check failed')
            if std_change is None or std_change > thresholds.convergence_K:
                reasons.append('standard deviation refinement check failed')
            # Coverage changes can accidentally hide a mean error; check its area too.
            coarse_supported = sum(integrations[i][0]['supported_area'] for i in face_ids)
            if abs(supported_area-coarse_supported)/exact_area > thresholds.coverage_change:
                reasons.append('supported area changed under refinement')
            source_values = temperatures[case_index, node_indices]
            results.append({'case_index': case_index, 'panel_id': panel['id'], 'panel_name': panel['name'],
                'area_m2': exact_area, 'integrated_area_m2': fine_area, 'supported_area_m2': supported_area,
                'coverage_fraction': coverage, 'extrapolated_area_m2': extrapolated,
                'temperature_avg_K': fine[0] if fine and not reasons else None,
                'temperature_std_K': fine[1] if fine and not reasons else None,
                'temperature_min_K': float(source_values.min()) if len(source_values) else None,
                'temperature_max_K': float(source_values.max()) if len(source_values) else None,
                'supported_temperature_avg_K': fine[0] if fine else None,
                'supported_temperature_std_K': fine[1] if fine else None,
                'node_count': len(node_indices),
                'quadrature_point_count': sum(len(integrations[i][1]['uv']) for i in face_ids),
                'area_relative_error': area_error,
                'refinement_mean_change_K': mean_change, 'refinement_std_change_K': std_change,
                'averaging_method': 'reconstructed_CAD_UV_linear_quadrature',
                'status': 'numerical_checks_passed_CAD_selection_unverified' if not reasons else 'UNRESOLVED',
                'issues': '; '.join(reasons)})
    return results


def visualization(path, faces, integrations, panels, rows, case_label=''):
    from html import escape
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'scene'}, {'type': 'scene'}]],
                        subplot_titles=('Candidate CAD panel assignments', 'Reconstructed averages (K)'))
    palette = ['#636efa', '#ef553b', '#00cc96', '#ab63fa', '#ffa15a', '#19d3f3', '#ff6692']
    by_id = {r['panel_id']: r for r in rows}
    valid_means = [r['temperature_avg_K'] for r in rows if r['temperature_avg_K'] is not None]
    selected_faces = {i: panel for panel in panels for i in panel['face_indices']}
    panel_positions = {p['id']: i for i,p in enumerate(panels)}
    for face in faces:
        uv_triangles = integrations[face.index][1]['triangles']
        # Bound the HTML size: a coarser CAD-only preview is sufficient to show
        # face assignments. Integration itself always uses the fine quadrature.
        if len(uv_triangles) > 2000:
            uv_triangles = face.uv_triangles
        xyz = np.array([[face.evaluate(p) for p in t] for t in uv_triangles]).reshape(-1, 3)
        panel = selected_faces.get(face.index)
        row = by_id.get(panel['id']) if panel else None
        label = escape(panel['name'] if panel else 'Explicitly excluded')
        hover = f'{label}<br>CAD inspection face {face.index}<br>'
        hover += (f"Coverage {row['coverage_fraction']:.2%}<br>{escape(row['status'])}" if row else 'Excluded')
        mean = row['temperature_avg_K'] if row else None
        if mean is not None:
            hover += f'<br>Reconstructed mean {mean:.6g} K'
        geometry = dict(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], i=np.arange(0,len(xyz),3),
                        j=np.arange(1,len(xyz),3), k=np.arange(2,len(xyz),3))
        for column in (1, 2):
            style = {'color': palette[panel_positions[panel['id']] % len(palette)] if panel else '#999999'}
            if column == 2:
                style = {'intensity': [mean]*len(xyz), 'coloraxis': 'coloraxis'} if mean is not None else {'color': '#999999'}
            fig.add_trace(go.Mesh3d(**geometry, **style, name=label, showlegend=False,
                hovertemplate=hover+'<extra></extra>', flatshading=True), row=1, col=column)
    low, high = (min(valid_means), max(valid_means)) if valid_means else (0, 1)
    if low == high:
        low, high = low-.5, high+.5
    fig.update_scenes(aspectmode='data', camera_eye=dict(x=3, y=3, z=2),
                      xaxis_title='CAD x (m)', yaxis_title='CAD y (m)', zaxis_title='CAD z (m)')
    fig.update_layout(height=850, coloraxis=dict(colorscale='Inferno', cmin=low, cmax=high,
                                              colorbar=dict(title='K')),
                      title=escape(case_label)+' — CAD reconstruction; candidate selections require review')
    html = fig.to_html(include_plotlyjs=True)
    note = '<p style="font:16px system-ui;padding:1em">Reconstructed CAD-surface averages, not original CFD element averages. Gray temperature faces are unresolved or explicitly excluded. Passing numerical checks does not verify CAD selections or temperature interpretation. Hover to inspect face IDs and panel coverage.</p>'
    path.write_text(html.replace('<body>', '<body>'+note), encoding='utf-8')


def run(args):
    output = args.output_dir
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError('--output-dir must be new or empty')
    if args.extrapolation_m > args.max_support_m:
        raise ValueError('--extrapolation-m must not exceed --max-support-m')
    if len({p.name for p in args.bc}) != len(args.bc):
        raise ValueError('BC filenames must be distinct for separate case output directories')
    if args.temperature_kinds and len(args.temperature_kinds) != len(args.bc):
        raise ValueError('--temperature-kinds must have one entry per BC file')
    print('Reading CAD and boundary-condition files...', flush=True)
    cases = []
    for path in args.bc:
        nodes = read_nodes(path, args.coordinate_unit, args.temperature_unit)
        points, temperatures, node_ids, duplicates = unique_samples(nodes, args.duplicate_tolerance_K)
        if cases and not np.array_equal(cases[0]['points'], points):
            raise ValueError('batch reconstruction requires identical unique coordinate sets; run different geometries separately')
        header = _header(next(_rows(path))[1])
        column = ('wall-temperature' if 'walltemperature' in header else 'temperature') if 'nodenumber' in header else 'headerless column 6'
        cases.append(dict(path=path, points=points, temperatures=temperatures, ids=node_ids, duplicates=duplicates, column=column,
                          row_count=len(nodes)))
    points, transform = load_transform(args.transform, args.step, cases[0]['points'])
    temperatures = np.array([c['temperatures'] for c in cases])
    faces = load_cad(args.step, args.deflection_m)
    cad_hash = sha256(args.step)
    panels, excluded = panels_for_faces(args.face_groups, faces, cad_hash)
    output.mkdir(parents=True, exist_ok=True)
    inventory = [{'face_index': f.index, 'area_m2': f.area, 'center_m': list(f.center),
                  'bbox_m': f.bbox.tolist()} for f in faces]
    write_json(output/'face_inventory.json', {'cad_sha256': cad_hash, 'faces': inventory,
        'warning': 'Indices refer to this STEP checksum and pinned CAD reader, not persistent SOLIDWORKS face IDs.'})
    write_json(output/'face_groups_template.json', {'cad_sha256': cad_hash, 'panels': panels,
                                                 'excluded_face_indices': excluded})
    mappings, membership, primary, distances = map_samples(points, faces, args.surface_tolerance_m,
                                                          lambda s: print(s, flush=True))
    write_csv(output/'source_mapping.csv', [dict(sample_index=i, source_x_m=source[0], source_y_m=source[1],
        source_z_m=source[2], cad_x_m=p[0], cad_y_m=p[1], cad_z_m=p[2],
        primary_face_index=int(primary[i]) or '', candidate_face_indices=';'.join(map(str,membership[i])),
        distance_m=float(distances[i]) if np.isfinite(distances[i]) else '',
        status='mapped' if membership[i] else ('AMBIGUOUS' if np.isfinite(distances[i]) else 'UNMAPPED'))
        for i, (source, p) in enumerate(zip(cases[0]['points'], points))])
    integrations, face_diagnostics = {}, []
    for face in faces:
        print(f'Integrating face {face.index}/{len(faces)} ({len(mappings[face.index])} source samples)...', flush=True)
        reason = ''
        try:
            interpolator = FaceInterpolator(face, mappings[face.index], points, temperatures, args.duplicate_tolerance_K)
        except ValueError as exc:
            interpolator, reason = None, str(exc)
        coarse = integrate(face, interpolator, temperatures, 2*args.edge_m, args.max_support_m, args.extrapolation_m)
        fine = integrate(face, interpolator, temperatures, args.edge_m, args.max_support_m, args.extrapolation_m)
        fine_edge = args.edge_m
        refinement_count = 0
        for _ in range(args.adaptive_levels):
            changes = [abs(a[key]-b[key]) for a,b in zip(coarse['cases'], fine['cases'])
                       if a is not None and b is not None for key in ('mean','std')]
            support_change = abs(fine['supported_area']-coarse['supported_area'])/face.area
            if (not changes or max(changes) <= args.convergence_K) and support_change <= args.coverage_change:
                break
            fine_edge /= 2
            print(f'  Refining face {face.index} to {fine_edge:g} m quadrature edges...', flush=True)
            coarse, fine = fine, integrate(face, interpolator, temperatures, fine_edge, args.max_support_m, args.extrapolation_m)
            refinement_count += 1
        integrations[face.index] = (coarse, fine)
        face_diagnostics.append({'face_index': face.index, 'source_sample_count': len(mappings[face.index]),
            'projected_collision_groups': interpolator.projected_collisions if interpolator else [],
            'interpolation_sample_count': len(interpolator.indices) if interpolator else 0,
            'fine_edge_m': fine_edge, 'coarse_edge_m': 2*fine_edge,
            'adaptive_refinement_count': refinement_count,
            'interpolator_issue': reason, 'exact_area_m2': face.area,
            'fine_integrated_area_m2': fine['area'], 'fine_supported_area_m2': fine['supported_area'],
            'fine_coverage_fraction': fine['coverage'], 'fine_extrapolated_area_m2': fine['extrapolated_area']})
    results = _panel_results(panels, faces, integrations, temperatures, mappings, args)
    mapped = distances[primary > 0]
    total_area = sum(f.area for f in faces)
    warnings = [
        'This is a reconstructed CAD-surface field, not an original CFD mesh area average.',
        'Face indices and candidate panel definitions require CAD review; no SOLIDWORKS selection was verified.',
        'Prescribed and adiabatic CFD temperatures are not independently predicted conducting-skin temperatures.',
        'No exterior skin thickness, internal conduction model, heat flux, or convection coefficient is established.',
        'Trim checks sample interpolation-triangle edges and centroids; narrow trim features may need further validation.',
        'No interpolation error bound follows from the quadrature refinement check.',
    ]
    if args.face_groups is None:
        warnings.append('Each CAD face is a candidate panel. Physical panel grouping has not been reviewed.')
    report = {'method': 'reconstructed_CAD_UV_linear_quadrature', 'tool_version': __version__,
        'dependency_versions': {'python': platform.python_version(), 'numpy': np.__version__,
                                'scipy':version('scipy'), 'cadquery-ocp':version('cadquery-ocp')},
        'created_UTC': datetime.now(timezone.utc).isoformat(),
        'coordinate_unit': 'm', 'temperature_unit': 'K', 'input_coordinate_unit': args.coordinate_unit,
        'input_temperature_unit': args.temperature_unit, 'cad_sha256': cad_hash,
        'coordinate_sha256': coordinate_hash(cases[0]['points']), 'transform': transform,
        'cad_selections_verified': False, 'source_point_count': len(points),
        'unmapped_source_count': int((primary == 0).sum()), 'shared_boundary_sample_count': sum(len(x)>1 for x in membership),
        'ambiguous_source_count': int(((primary == 0) & np.isfinite(distances)).sum()),
        'mapped_source_distance_median_p95_max_m': np.quantile(mapped, [.5,.95,1]).tolist() if len(mapped) else [],
        'all_source_points_checked_against_surface_tolerance': True,
        'mapped_points_use_exact_face_distances': True, 'surface_tolerance_m': args.surface_tolerance_m,
        'cad_area_m2': total_area, 'excluded_face_indices': excluded,
        'fine_max_triangle_edge_m': args.edge_m, 'coarse_max_triangle_edge_m': 2*args.edge_m,
        'max_adaptive_refinements': args.adaptive_levels,
        'max_interpolation_vertex_distance_m': args.max_support_m,
        'max_boundary_extrapolation_distance_m': args.extrapolation_m,
        'min_panel_coverage': args.min_coverage, 'max_extrapolated_fraction': args.max_extrapolated_fraction,
        'max_area_relative_error': args.max_area_error, 'max_mean_std_refinement_change_K': args.convergence_K,
        'max_supported_area_refinement_change_fraction': args.coverage_change,
        'face_diagnostics': face_diagnostics, 'warnings': warnings, 'cases': []}
    common_inputs = output/'inputs'
    common_inputs.mkdir()
    shutil.copyfile(args.step, common_inputs/'geometry.step')
    shutil.copyfile(args.transform, common_inputs/'transform.json')
    if args.face_groups:
        shutil.copyfile(args.face_groups, common_inputs/'face_groups.json')
    unresolved_total = 0
    for case_index, case in enumerate(cases):
        directory = output/case['path'].name
        directory.mkdir()
        rows = [{k:v for k,v in row.items() if k != 'case_index'} for row in results if row['case_index'] == case_index]
        unresolved = [r['panel_id'] for r in rows if r['status'] == 'UNRESOLVED']
        unresolved_total += len(unresolved)
        write_csv(directory/'panel_temperatures.csv', rows)
        write_csv(directory/'solidworks_mapping.csv', [{'panel_id': p['id'], 'panel_name': p['name'],
            'cad_face_indices': ';'.join(map(str,p['face_indices'])), 'solidworks_component': p['solidworks_component'],
            'solidworks_region': p['solidworks_region'], 'selection_notes': p['selection_notes'], 'cad_sha256':cad_hash} for p in panels])
        visualization(directory/'panel_visualization.html', faces, integrations, panels, rows, case['path'].name)
        shutil.copyfile(case['path'], directory/'original_boundary_conditions.txt')
        write_json(directory/'node_ids.json', {'sample_index_to_original_node_ids': case['ids'], 'coordinate_duplicates': case['duplicates']})
        report['cases'].append({'filename': case['path'].name, 'sha256': sha256(case['path']),
            'input_row_count': case['row_count'], 'duplicate_coordinate_groups': case['duplicates'],
            'temperature_column':case['column'],
            'temperature_kind': args.temperature_kinds[case_index] if args.temperature_kinds else 'unknown',
            'unresolved_panels': unresolved, 'panel_count':len(rows)})
    report['all_numerical_checks_passed'] = unresolved_total == 0 and not report['unmapped_source_count']
    write_json(output/'reconstruction_report.json', report)
    print(f'Wrote {len(cases)} cases to {output}. Unresolved panel/case rows: {unresolved_total}; unmapped source points: {report["unmapped_source_count"]}.', flush=True)
    print('CAD face selections and physical interpretation still require review.', flush=True)
    return 0 if report['all_numerical_checks_passed'] or args.allow_partial else 3


def main(argv=None):
    parser = argparse.ArgumentParser(description='Reconstruct temperature on STEP faces, with coverage and two-level integration checks.')
    parser.add_argument('--step', type=Path, required=True)
    parser.add_argument('--bc', type=Path, nargs='+', required=True, help='one or more BC files with identical unique coordinate sets')
    parser.add_argument('--transform', type=Path, required=True, help='rigid BC-to-CAD transform JSON, translation in meters')
    parser.add_argument('--face-groups', type=Path, help='reviewed groups from face_groups_template.json; default each CAD face is a candidate panel')
    parser.add_argument('--coordinate-unit', choices=tuple(COORDINATE_SCALES), default='m')
    parser.add_argument('--temperature-unit', choices=('K','C','F'), default='K')
    parser.add_argument('--temperature-kinds', choices=('unknown','prescribed','adiabatic','structural'), nargs='+')
    parser.add_argument('--surface-tolerance-m', type=_positive, default=.0005)
    parser.add_argument('--edge-m', type=_positive, default=.01, help='fine maximum CAD quadrature triangle edge (coarse is twice this)')
    parser.add_argument('--deflection-m', type=_positive, default=.0001)
    parser.add_argument('--adaptive-levels', type=int, choices=range(4), default=2,
                        help='additional face integration refinements when convergence checks fail (0-3, default 2)')
    parser.add_argument('--max-support-m', type=_positive, default=.04)
    parser.add_argument('--extrapolation-m', type=_positive, default=.003)
    parser.add_argument('--min-coverage', type=_fraction, default=.98)
    parser.add_argument('--max-extrapolated-fraction', type=_fraction, default=.02)
    parser.add_argument('--max-area-error', type=_fraction, default=.01)
    parser.add_argument('--convergence-K', type=_positive, default=1.0)
    parser.add_argument('--coverage-change', type=_fraction, default=.005)
    parser.add_argument('--duplicate-tolerance-K', type=_positive, default=1e-6)
    parser.add_argument('--allow-partial', action='store_true', help='exit 0 with unresolved rows (they stay blank and explicitly flagged)')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        parser.exit(2, f'error: {exc}\n')


if __name__ == '__main__':
    raise SystemExit(main())
