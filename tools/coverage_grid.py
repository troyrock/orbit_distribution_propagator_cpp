#!/usr/bin/env python3
"""Resumable parameter-grid orchestration; all physics executes in C++."""
import argparse
import concurrent.futures
from contextlib import contextmanager
import csv
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import time
import uuid

FIELDS = ('run_id', 'perigee_altitude_km', 'inclination_deg', 'eccentricity',
          'position_sigma_km', 'velocity_sigma_m_s')
AXES = ('perigee_altitudes_km', 'inclinations_deg', 'eccentricities',
        'position_sigmas_km', 'velocity_sigmas_m_s')


@contextmanager
def output_lock(directory):
    """Prevent two writers; OS locks are released even after process failure."""
    with (directory / '.coverage-grid.lock').open('a+b') as stream:
        if stream.seek(0, 2) == 0:
            stream.write(b'0')
            stream.flush()
        stream.seek(0)
        if os.name == 'nt':
            import msvcrt
            acquire = lambda: msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            release = lambda: msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            acquire = lambda: fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            release = lambda: fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
        try:
            acquire()
        except OSError as error:
            raise RuntimeError('Study directory is locked by another coverage-grid process') from error
        try:
            yield
        finally:
            stream.seek(0)
            release()


def numbers(text):
    values = [float(x) for x in text.split(',')]
    if not values or any(not math.isfinite(x) for x in values):
        raise argparse.ArgumentTypeError('Axes require finite comma-separated numbers')
    if values != sorted(set(values)):
        raise argparse.ArgumentTypeError('Axes must be strictly increasing')
    return values


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def make_cases(axes):
    return [dict(zip(FIELDS, (f'case_{index:06d}', *values)))
            for index, values in enumerate(itertools.product(*axes))]


def load_completed(directory, cases):
    expected = {case['run_id']: case for case in cases}
    completed = {}
    for path in sorted(directory.glob('batch_*.jsonl')):
        for line in path.read_text(encoding='utf-8').splitlines():
            if not line.strip():
                continue
            row = json.loads(line)  # Fail loudly on an incomplete/corrupt checkpoint.
            name = row['run_id']
            if name not in expected or name in completed:
                raise ValueError(f'Unexpected or duplicate completed case: {name}')
            if any(row[field] != expected[name][field] for field in FIELDS[1:]):
                raise ValueError(f'Checkpoint parameters do not match: {name}')
            completed[name] = row
    return completed


def prepare_reuse(source, destination, manifest, cases):
    """Validate a complete source study before adopting matching parameter tuples."""
    source = source.resolve(strict=True)
    if source == destination.resolve():
        raise ValueError('Reuse source must differ from the output directory')
    manifest_bytes = (source / 'manifest.json').read_bytes()
    study_bytes = (source / 'study.json').read_bytes()
    source_manifest = json.loads(manifest_bytes)
    document = json.loads(study_bytes)
    for key in ('schema_version', 'executable_sha256', 'config_sha256', 'config_text'):
        if source_manifest.get(key) != manifest[key]:
            raise ValueError(f'Reuse rejected: source {key} differs from the new study')
    source_metadata = source_manifest.get('metadata', {})
    # Grid axes and presentation interpolation may change; the physics contract may not.
    presentation = set(AXES) | {'interpolation', 'altitude_interpolation'}
    for key, value in manifest['metadata'].items():
        if key not in presentation and source_metadata.get(key) != value:
            raise ValueError(f'Reuse rejected: source metadata {key} differs')
    for key, value in source_metadata.items():
        if document.get('metadata', {}).get(key) != value:
            raise ValueError(f'Reuse rejected: source study metadata {key} differs from its manifest')
    # Loading by path keeps this standalone script importable from test tools too.
    spec = importlib.util.spec_from_file_location(
        'coverage_reuse_validation', Path(__file__).with_name('validate_coverage_grid.py'))
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    validator.validate_document(document, str(source))
    by_parameters = {tuple(row[field] for field in FIELDS[1:]): row
                     for row in document['runs']}
    adopted = []
    for case in cases:
        row = by_parameters.get(tuple(case[field] for field in FIELDS[1:]))
        if row is not None:
            adopted.append(dict(row, run_id=case['run_id']))
    provenance = {'source_directory': str(source),
                  'manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(),
                  'study_sha256': hashlib.sha256(study_bytes).hexdigest(),
                  'adopted_cases': len(adopted)}
    return provenance, adopted


def initialize_manifest(directory, raw, manifest, cases, reuse_study):
    """Install an atomic reused checkpoint, retaining ordinary offline resume."""
    path = directory / 'manifest.json'
    checkpoint = raw / 'batch_reused.jsonl'
    adopted = None
    if path.exists():
        saved = json.loads(path.read_text(encoding='utf-8'))
        if 'reuse_study' in saved:
            manifest['reuse_study'] = saved['reuse_study']
        if saved != manifest:
            raise ValueError('Resume rejected: executable, configuration or axes differ from manifest')
        provenance = saved.get('reuse_study')
        if reuse_study is not None:
            if provenance is None or reuse_study.resolve() != Path(provenance['source_directory']):
                raise ValueError('Resume rejected: reuse source differs from manifest')
        if provenance is not None and (reuse_study is not None or not checkpoint.exists()):
            verified, adopted = prepare_reuse(Path(provenance['source_directory']),
                                              directory, manifest, cases)
            if verified != provenance:
                raise ValueError('Resume rejected: reuse source content differs from manifest')
    else:
        if list(raw.iterdir()):
            raise ValueError('Cannot adopt existing raw data without a matching manifest')
        if reuse_study is not None:
            manifest['reuse_study'], adopted = prepare_reuse(reuse_study, directory, manifest, cases)
        atomic_json(path, manifest)
    if adopted is not None:
        if checkpoint.exists():
            existing = [json.loads(line) for line in checkpoint.read_text(encoding='utf-8').splitlines()
                        if line.strip()]
            if existing != adopted:
                raise ValueError('Resume rejected: adopted checkpoint differs from source study')
        else:
            temporary = checkpoint.with_suffix('.jsonl.tmp')
            temporary.write_text(''.join(json.dumps(row, allow_nan=False) + '\n' for row in adopted),
                                 encoding='utf-8')
            temporary.replace(checkpoint)
    return manifest


def export(directory, manifest, completed, cases):
    ordered = [completed[c['run_id']] for c in cases if c['run_id'] in completed]
    metadata = dict(manifest['metadata'], completed_cases=len(ordered), expected_cases=len(cases),
                    complete=len(ordered) == len(cases))
    atomic_json(directory / 'study.json', {'metadata': metadata, 'runs': ordered})
    if ordered:
        with (directory / 'coverage_grid.csv').open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(ordered[0]))
            writer.writeheader()
            writer.writerows(ordered)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', required=True, type=Path)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--reuse-study', type=Path,
                        help='Reuse compatible records from a completed grid directory by parameter tuple')
    parser.add_argument('--altitude-interpolation', default='linear',
                        choices=('linear', 'log_geocentric_perigee_radius'),
                        help='Altitude coordinate for visualization interpolation (default: legacy linear)')
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--batch-size', type=int, default=16)
    parser.add_argument('--altitudes', type=numbers, default=numbers('1000,2000,3000'))
    parser.add_argument('--inclinations', type=numbers, default=numbers('0,45,90'))
    parser.add_argument('--eccentricities', type=numbers, default=numbers('0,0.05,0.1'))
    parser.add_argument('--positions', type=numbers, default=[10 ** (k / 5) for k in range(11)])
    parser.add_argument('--velocities', type=numbers, default=[10 ** (-2 + k / 5) for k in range(11)])
    args = parser.parse_args()
    if min(args.jobs, args.threads, args.batch_size) < 1:
        parser.error('jobs, threads and batch-size must be positive')
    axes = [args.altitudes, args.inclinations, args.eccentricities, args.positions, args.velocities]
    if (min(axes[0]) < 0 or min(axes[1]) < 0 or max(axes[1]) > 90 or
            min(axes[2]) < 0 or max(axes[2]) >= 1 or min(axes[3]) <= 0 or min(axes[4]) <= 0):
        parser.error('Axes outside the supported prograde, elliptic, positive-sigma domain')
    args.exe = args.exe.resolve(strict=True)
    args.config = args.config.resolve(strict=True)
    args.output.mkdir(parents=True, exist_ok=True)
    with output_lock(args.output):
        return run_study(args, axes)


def run_study(args, axes):
    raw = args.output / 'raw'
    raw.mkdir(exist_ok=True)
    config = {}
    for line in args.config.read_text(encoding='utf-8').splitlines():
        line = line.split('#', 1)[0].strip()
        if '=' in line:
            key, value = line.split('=', 1)
            config[key.strip()] = value.strip()
    # The explorer contract is deliberately explicit; do not mislabel arbitrary configs.
    for key, value in {'samples': '20000', 'phase_bins': '72', 'coverage_max_gap_deg': '5',
                       'coverage_occupied_fraction': '1', 'persistence': '3',
                       'force_model': 'j2_j2sq', 'initial_type': 'osculating'}.items():
        if config.get(key) != value:
            raise ValueError(f'Explorer requires {key} = {value} explicitly in base config')
    manifest = {
        'schema_version': 1, 'executable_sha256': digest(args.exe),
        'config_sha256': digest(args.config), 'config_text': args.config.read_text(encoding='utf-8'),
        'metadata': {
            **dict(zip(AXES, axes)), 'samples': 20000, 'seed': int(config['seed']),
            'coverage_definition': 'All 72 relative mean-longitude bins occupied and maximum circular gap <= 5 degrees for 3 consecutive epochs; reported time is first qualifying epoch.',
            'model_scope': 'Native C++ DSST degree-2 J2 plus Zeis J2-squared, gravity only; exact constant mean-longitude rates after native osculating-to-mean conversion. No drag, higher harmonics, third bodies, or radiation pressure.',
            'uncertainty_convention': 'Independent Gaussian per-axis 1-sigma in RTN position and velocity; no initial correlations; common random seed.',
            'altitude_definition': 'Nominal perigee altitude above the equatorial Earth radius.',
            'earth_radius_km': float(config.get('earth_radius_m', '6378137')) / 1000,
            'raan_deg': 20, 'argument_of_perigee_deg': 30, 'mean_anomaly_deg': 10,
            'interpolation': 'Trilinear interpolation of log coverage time between orbital nodes; unavailable corners preserve holes.',
            'domain_policy': 'Mask whole ensembles containing initial or mean Earth-crossing ellipses. No particles clipped, removed, or redrawn. Above-Earth low-perigee tails remain idealized gravity-only cases.',
        },
    }
    # Omitting the default key preserves manifests and resume for legacy studies.
    if args.altitude_interpolation != 'linear':
        manifest['metadata']['altitude_interpolation'] = args.altitude_interpolation
    cases = make_cases(axes)
    manifest = initialize_manifest(args.output, raw, manifest, cases, args.reuse_study)
    completed = load_completed(raw, cases)
    pending = [case for case in cases if case['run_id'] not in completed]
    # Complete uncertainty corners and center across all orbital settings first.
    # These nodes permit an early independent midpoint-interpolation audit.
    position_anchors = {args.positions[k] for k in (0, len(args.positions) // 2, -1)}
    velocity_anchors = {args.velocities[k] for k in (0, len(args.velocities) // 2, -1)}
    pending.sort(key=lambda case: not (case['position_sigma_km'] in position_anchors and
                                      case['velocity_sigma_m_s'] in velocity_anchors))
    chunks = [pending[index:index + args.batch_size] for index in range(0, len(pending), args.batch_size)]
    export(args.output, manifest, completed, cases)
    start = time.perf_counter()
    initially_complete = len(completed)

    def execute(chunk):
        stem = raw / ('batch_' + uuid.uuid4().hex)
        with stem.with_suffix('.csv').open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(chunk)
        command = [str(args.exe), '--config', str(args.config), '--cases', str(stem.with_suffix('.csv')),
                   '--output', str(stem.with_suffix('.jsonl')), '--threads', str(args.threads)]
        with stem.with_suffix('.log').open('w', encoding='utf-8') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        return [json.loads(line) for line in stem.with_suffix('.jsonl').read_text(encoding='utf-8').splitlines()]

    print(f'{len(cases)} cases; {len(completed)} already complete; {args.jobs} jobs x {args.threads} threads', flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(execute, chunk) for chunk in chunks}
        try:
            while futures:
                done, futures = concurrent.futures.wait(futures, timeout=30,
                                                        return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    for row in future.result():
                        completed[row['run_id']] = row
                if done:
                    export(args.output, manifest, completed, cases)
                elapsed = time.perf_counter() - start
                finished = len(completed) - initially_complete
                remaining = (len(cases) - len(completed)) * elapsed / finished if finished else None
                print(json.dumps({'completed': len(completed), 'total': len(cases),
                                  'percent': round(100 * len(completed) / len(cases), 1),
                                  'elapsed_seconds': round(elapsed, 1),
                                  'eta_seconds': round(remaining, 1) if remaining is not None else None}), flush=True)
        except BaseException:
            for future in futures:
                future.cancel()
            raise
    export(args.output, manifest, completed, cases)


if __name__ == '__main__':
    main()
