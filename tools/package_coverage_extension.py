#!/usr/bin/env python3
"""Merge compatible native coverage grids and package a validated offline report."""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import zipfile

from render_coverage_explorer import find_plotly, render
from validate_coverage_grid import AXES, FIELDS, validate_document, validate_grid


POLICY = 'log_geocentric_perigee_radius'
SCIENCE_METADATA = (*AXES[1:], 'samples', 'seed', 'coverage_definition', 'model_scope',
                    'uncertainty_convention', 'altitude_definition', 'earth_radius_km',
                    'raan_deg', 'argument_of_perigee_deg', 'mean_anomaly_deg')


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def merge_studies(documents, manifests):
    """Reuse measurements only when force settings, executable and axes agree."""
    if not documents or len(documents) != len(manifests):
        raise ValueError('Provide one manifest for each nonempty source study')
    first = documents[0]['metadata']
    for document, manifest in zip(documents, manifests):
        validate_document(document, 'source')
        for key in ('executable_sha256', 'config_sha256', 'config_text'):
            if not manifest.get(key) or manifest[key] != manifests[0].get(key):
                raise ValueError('Source provenance differs for ' + key)
        for key in SCIENCE_METADATA:
            if document['metadata'].get(key) != first.get(key):
                raise ValueError('Source metadata differs for ' + key)
        for key in AXES:
            if document['metadata'][key] != manifest['metadata'][key]:
                raise ValueError('Manifest axes differ from study: ' + key)
        for key in SCIENCE_METADATA:
            if document['metadata'].get(key) != manifest['metadata'].get(key):
                raise ValueError('Manifest metadata differs from study: ' + key)
    meta = copy.deepcopy(first)
    meta['perigee_altitudes_km'] = sorted({h for doc in documents for h in doc['metadata'][AXES[0]]})
    meta['altitude_interpolation'] = POLICY
    meta['interpolation'] = ('Interpolate log coverage time in log geocentric perigee radius, '
                             'inclination and eccentricity; unavailable corners preserve holes.')
    rows, seen = [], set()
    for source_index, document in enumerate(documents):
        for original in document['runs']:
            key = tuple(original[name] for name in FIELDS)
            if key in seen:
                raise ValueError('Source studies overlap at ' + repr(key))
            seen.add(key)
            row = dict(original, source_study=source_index, source_run_id=original['run_id'])
            rows.append(row)
    rows.sort(key=lambda row: tuple(row[name] for name in FIELDS))
    for index, row in enumerate(rows):
        row['run_id'] = f'extended_{index:06d}'
    meta.update(expected_cases=math.prod(len(meta[key]) for key in AXES),
                completed_cases=len(rows), complete=True)
    combined = {'metadata': meta, 'runs': rows}
    validate_document(combined, 'combined')
    return combined


def original_reference(path):
    """Keep the PDF's actual measured values and audit its initial tail census."""
    path = Path(path)
    original = read_json(path)
    source = original['metadata']
    expected = {'samples': 20000, 'seed': 20260919, 'force_model': 'j2_j2sq',
                'initial_type': 'osculating', 'output_type': 'mean',
                'orbit_keplerian_deg': [26560000, 0.02, 55, 20, 30, 10]}
    for key, value in expected.items():
        if source.get(key) != value:
            raise ValueError('Original MEO reference metadata differs for ' + key)
    positions, velocities = source['position_sigmas_km'], source['velocity_sigmas_m_s']
    for axis, endpoints in ((positions, (1, 100)), (velocities, (0.01, 1))):
        if (len(axis) != 11 or axis[0] != endpoints[0] or axis[-1] != endpoints[1]
                or any(not isinstance(value, (int, float)) or isinstance(value, bool)
                       or not math.isfinite(value) for value in axis)
                or any(left >= right for left, right in zip(axis, axis[1:]))):
            raise ValueError('Original MEO reference requires its 11-by-11 uncertainty axes')
    expected_points = {(position, velocity) for position in positions for velocity in velocities}
    points = [(row['position_sigma_km'], row['velocity_sigma_m_s']) for row in original['runs']]
    identifiers = [row['run_id'] for row in original['runs']]
    if (len(points) != 121 or set(points) != expected_points or
            len(set(identifiers)) != len(identifiers)):
        raise ValueError('Original MEO reference requires 121 unique measured cases')
    a, e, inc, raan, arg, anomaly = source['orbit_keplerian_deg']
    radius = 6378.137
    h = a / 1000 * (1 - e) - radius
    metadata = dict(source, perigee_altitudes_km=[h], inclinations_deg=[inc],
                    eccentricities=[e], earth_radius_km=radius, raan_deg=raan,
                    argument_of_perigee_deg=arg, mean_anomaly_deg=anomaly,
                    original_study_sha256=sha256(path))
    rows = []
    for original_row in original['runs']:
        if any(original_row.get(key) != source[key] for key in ('samples', 'seed')):
            raise ValueError('Reference particle count or seed differs from metadata')
        time, lower, upper, confirmation = [original_row.get(key) for key in (
            'coverage_time_days', 'coverage_lower_days', 'coverage_upper_days',
            'coverage_confirmation_days')]
        if (any(not isinstance(value, (int, float)) or isinstance(value, bool)
                or not math.isfinite(value) for value in (time, lower, upper, confirmation))
                or not 0 <= lower <= time == upper <= confirmation or time == 0):
            raise ValueError('Reference contains invalid measured coverage bounds')
        row = dict(original_row, status='ok', perigee_altitude_km=h,
                   inclination_deg=inc, eccentricity=e)
        samples_path = path.parent / original_row['output_directory'] / 'initial_samples.csv'
        counts, minimum, samples = [0, 0, 0], math.inf, 0
        with samples_path.open(encoding='utf-8', newline='') as stream:
            for particle in csv.DictReader(stream):
                perigee = float(particle['a_m']) / 1000 * (
                    1 - math.hypot(float(particle['ex']), float(particle['ey']))) - radius
                minimum = min(minimum, perigee)
                counts = [count + (perigee < limit) for count, limit in zip(counts, (300, 500, 1000))]
                samples += 1
        if samples != source['samples']:
            raise ValueError('Reference sample census differs for ' + row['run_id'])
        row.update(min_initial_perigee_km=minimum, particles_below_300km=counts[0],
                   particles_below_500km=counts[1], particles_below_1000km=counts[2])
        rows.append(row)
    return {'label': 'Original MEO / PDF',
            'description': ('The original 121 measured PDF cases: a = 26,560 km, '
                            'perigee = 19,650.663 km, inclination = 55 degrees, eccentricity = 0.02. '
                            'These retain the original output cadence and are separate from slider interpolation.'),
            'metadata': metadata, 'runs': rows}


def audit_regression_holdouts(data, directory):
    """Audit prior direct measurements, retaining their distinct executable provenance."""
    document = read_json(directory / 'study.json')
    manifest = read_json(directory / 'manifest.json')
    for key in (*AXES, *SCIENCE_METADATA):
        if document['metadata'].get(key) != manifest['metadata'].get(key):
            raise ValueError('Regression manifest metadata differs from study: ' + key)
    # These holdouts may predate a validated implementation optimization. They
    # must match scientific settings, but their executable hash is not rewritten.
    report = validate_grid(data, document)
    report['source'] = {'directory': str(directory.resolve()),
                        'study_sha256': sha256(directory / 'study.json'),
                        'manifest_sha256': sha256(directory / 'manifest.json'),
                        'executable_sha256': manifest.get('executable_sha256'),
                        'config_sha256': manifest.get('config_sha256')}
    report['provenance_note'] = (
        'Previously simulated direct holdouts, compared under the current interpolation policy. '
        'Scientific metadata and seed must match; the historical executable hash is retained '
        'and is not required to equal the current optimized implementation.')
    return report


def package(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    documents = [read_json(path / 'study.json') for path in args.study_dirs]
    manifests = [read_json(path / 'manifest.json') for path in args.study_dirs]
    data = merge_studies(documents, manifests)
    holdouts = read_json(args.holdouts / 'study.json')
    holdout_manifest = read_json(args.holdouts / 'manifest.json')
    for key in ('executable_sha256', 'config_sha256', 'config_text'):
        if holdout_manifest.get(key) != manifests[0].get(key):
            raise ValueError('Holdout provenance differs for ' + key)
    validation = validate_grid(data, holdouts)
    write_json(output / 'validation.json', validation)
    if not validation['passed'] or validation['provisional']:
        raise ValueError('Complete, passing withheld interpolation validation is required')
    regression = None
    if args.regression_holdouts is not None:
        regression = audit_regression_holdouts(data, args.regression_holdouts)
        write_json(output / 'regression_validation.json', regression)
        if not regression['passed'] or regression['provisional']:
            raise ValueError('Complete, passing legacy holdout regression validation is required')
    reference = original_reference(args.reference)
    data['reference_surface'] = reference
    data['metadata']['initial_view'] = 'reference'
    stats = validation['statistics']
    sources = [{'path': str(path.resolve()), 'study_sha256': sha256(path / 'study.json'),
                'manifest_sha256': sha256(path / 'manifest.json'), 'cases': len(doc['runs'])}
               for path, doc in zip(args.study_dirs, documents)]
    data['metadata']['validation_summary'] = stats
    if regression is not None:
        data['metadata']['regression_validation_summary'] = regression['statistics']
    data['metadata']['source_studies'] = sources
    data['metadata']['provenance'] = (
        f"{len(data['runs']):,} C++ ensembles at 20,000 particles each; "
        f"{len(holdouts['runs'])} withheld ensembles, {stats['n']} supported comparisons. "
        f"Maximum withheld interpolation error {100 * stats['max_relative_error']:.3f}%. "
        'The original MEO/PDF measurements are available as a separate reference surface.')
    valid = [row for row in data['runs'] if row['status'] == 'ok']
    summary = {'main_ensembles': len(data['runs']), 'particles_per_ensemble': 20000,
               'main_particle_draws': len(data['runs']) * 20000,
               'reference_ensembles': len(reference['runs']),
               'holdout_ensembles': len(holdouts['runs']),
               'status_counts': dict(Counter(row['status'] for row in data['runs'])),
               'minimum_time_case': min(valid, key=lambda row: row['coverage_time_days']),
               'maximum_time_case': max(valid, key=lambda row: row['coverage_time_days']),
               'reference_maximum_days': max(row['coverage_time_days'] for row in reference['runs']),
               'interpolation_validation': stats, 'source_studies': sources,
               'executable_sha256': manifests[0]['executable_sha256'],
               'holdouts_directory': str(args.holdouts.resolve()),
               'native_scope': data['metadata']['model_scope']}
    if regression is not None:
        summary['regression_validation'] = {'statistics': regression['statistics'],
                                            'source': regression['source']}
    write_json(output / 'study.json', data)
    write_json(output / 'summary.json', summary)
    (output / 'explorer.html').write_text(render(data, find_plotly(args.plotly_js)), encoding='utf-8')
    with (output / 'coverage_grid.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data['runs'][0]))
        writer.writeheader()
        writer.writerows(data['runs'])
    regression_note = ''
    if regression is not None:
        regression_stats = regression['statistics']
        regression_note = (
            f"\nThe earlier {regression['holdouts']['validated_cases']} direct holdouts were also "
            f"rechecked using the new interpolation coordinate. The {regression_stats['n']} "
            f"supported comparisons passed the unchanged 5% gate, with maximum error "
            f"{100 * regression_stats['max_relative_error']:.3f}%. Their original executable "
            "and source-file hashes remain in regression_validation.json; these are reused "
            "measurements, not new simulations.\n")
    readme = f'''# Orbital coverage explorer: 1,000-30,000 km perigee altitude

Open **explorer.html** in a browser. This file embeds all data and plotting code;
it works offline and can be shared on its own.

## Controls and comparison with the original PDF

The altitude control is **perigee altitude**, the lowest nominal altitude above
Earth's equatorial radius of 6,378.137 km. Inclination spans 0-90 degrees and
eccentricity spans 0-0.1. Position and velocity uncertainties span 1-100 km and
0.01-1 m/s, with 11 logarithmic nodes on each uncertainty axis.

Use **Original MEO / PDF** to load the original 121 measured cases exactly:
semimajor axis 26,560 km, perigee altitude 19,650.663 km, inclination 55 degrees,
eccentricity 0.02. Its maximum is {summary['reference_maximum_days']:.6f} days.
Moving an orbital slider leaves that reference and selects the expanded grid.
The reference retains the PDF's output epochs; it is not silently mixed into
the new grid. Matching the orbit manually instead displays a labeled interpolation.

The plot can fit its height/color scale to the selected surface or keep a common
scale for comparisons. Snap to simulated orbital nodes for measured surfaces.
Hover for coverage times and onset brackets; export the current surface as CSV.

## Measurements and interpolation

Measured perigee nodes: {', '.join(f'{h:g}' for h in data['metadata']['perigee_altitudes_km'])} km.
Inclination nodes: 0, 45, 90 degrees. Eccentricity nodes: 0, 0.05, 0.1.
The full grid contains {len(data['runs']):,} cases, {len(data['runs']) * 20000:,}
particle draws, and status counts {json.dumps(summary['status_counts'], sort_keys=True)}.
Measured coverage times range from {summary['minimum_time_case']['coverage_time_days']:.6f}
to {summary['maximum_time_case']['coverage_time_days']:.6f} days.

Intermediate orbits interpolate log time using log(Earth radius + perigee altitude),
linear inclination and linear eccentricity. This respects the leading power-law
dependence on orbit size. It does not replace the native propagation at measured
nodes. A missing/masked supporting corner remains a gap.

The {len(holdouts['runs'])} independently simulated holdouts check original and
refined altitude intervals, inclination 22.5/67.5 degrees, eccentricity .025/.075,
and position sigma 1/10/100 km by velocity sigma .01/.1/1 m/s. Altitude checks
use geometric midpoints in geocentric perigee radius. An original altitude
midpoint can become a measured altitude after refinement; its independent
inclination/eccentricity holdouts remain in the audit alongside new half-interval
checks. Refinement does not discard earlier failing comparisons or widen the gate.
The {stats['n']} supported comparisons have maximum relative error
{100 * stats['max_relative_error']:.3f}%, 95th percentile
{100 * stats['p95_relative_error']:.3f}%, and RMS
{100 * stats['rms_relative_error']:.3f}%. Validation applies to these sampled
holdouts, not every continuously selectable point. See validation.json.
{regression_note}

## Scientific conventions

Every case uses 20,000 independent Gaussian draws with per-axis 1-sigma RTN
uncertainties, no initial correlations, seed 20260919, osculating initialization,
RAAN 20 degrees, argument of perigee 30 degrees, and mean anomaly 10 degrees.
Coverage requires all 72 relative mean-longitude bins occupied and maximum gap
at most 5 degrees for three consecutive saved epochs. The reported time is the
first qualifying epoch; temporal brackets are not Monte Carlo confidence intervals.
Coverage does not imply uniform phase density or a particular tube thickness.

The native C++ DSST model includes J2 and J2-squared only. Its mean-phase flow is
evaluated exactly after native osculating-to-mean initialization. Drag, higher
harmonics, Sun/Moon gravity, radiation pressure, maneuvers, and process noise
are omitted. Long-duration high-altitude cases remain conditional experiments
under this force model. Earth-crossing ensembles are masked whole; no particles
are discarded or redrawn. Above-Earth low-perigee tails are explicitly counted.

## Provenance

Existing measurements are reused byte-for-byte in their scientific fields;
source_study/source_run_id retain their origin. Executable and config hashes
must match before merging. The executable SHA-256 is
`{summary['executable_sha256']}`. Source paths and file hashes are in summary.json.
All particle propagation uses C++; Python orchestrates, audits and packages it.
The repository tool tools/package_coverage_extension.py reproduces this package
from source study directories, the holdout directory, and the original MEO study.

Files: explorer.html (standalone), study.json (embedded dataset), summary.json,
coverage_grid.csv (measured grid), validation.json, checksums.sha256, explorer.zip.
'''
    (output / 'README.md').write_text(readme, encoding='utf-8')
    names = ('explorer.html', 'study.json', 'summary.json', 'coverage_grid.csv', 'validation.json', 'README.md')
    if regression is not None:
        names += ('regression_validation.json',)
    (output / 'checksums.sha256').write_text(''.join(sha256(output / name) + '  ' + name + '\n' for name in names), encoding='utf-8')
    with zipfile.ZipFile(output / 'explorer.zip', 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for name in ('explorer.html', 'README.md'):
            archive.write(output / name, name)
    print(json.dumps({key: value for key, value in summary.items() if key != 'source_studies'}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study-dirs', nargs='+', type=Path, required=True)
    parser.add_argument('--holdouts', type=Path, required=True)
    parser.add_argument('--regression-holdouts', type=Path,
                        help='Optional historical direct holdouts, also required to pass the 5%% gate')
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plotly-js', type=Path)
    package(parser.parse_args())


if __name__ == '__main__':
    main()
