#!/usr/bin/env python3
"""Acceptance check for campaign era 0831_a1 (fixed-algorithm rerun).

0831_a1 pins every published Q-less row to RANDLAPACK_CHOL_MAX_RETRIES=0, so a
Cholesky breakdown must surface as a FAIL row (qr_status != 0) instead of
silently measuring the shift-rescued variant. This script states the prediction
made from 0829_a1 BEFORE the rerun and checks the new data against it, so the
outcome is a pass/fail verdict rather than a post-hoc reading of the figures.

Prediction, per row, from the 0829_a1 baseline:
  chol_retries > 0  in 0829  ->  qr_status != 0 (FAIL) in 0831
  chol_retries == 0 in 0829  ->  qr_status == 0 and numbers match 0829

The one row exempt from the prediction is the FEM2 native_ill CholQR family:
its POTRF verdict sits on the rounding boundary and is decided by BLAS
summation order (see reports/2026-08-31-iterations-vs-walltime-canonical.md),
so either outcome is admissible there and is reported, not failed.

Usage:  python3 verify_0831_fixed_algorithm.py [results_dir]
"""
import csv
import glob
import os
import sys

RESULTS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', 'results')

EXPECTED_COMMIT = '92061e2'
EXPECTED_KNOB = 'RANDLAPACK_CHOL_MAX_RETRIES=0'
BOUNDARY_EXEMPT = {('native_ill_dd', 'CholQR'), ('native_ill_dd', 'CholQR2')}

FAMILIES = [
    ('toeplitz_ls_{era}_pcg_ne', ['small', 'fixedm', 'middle', 'large'],
     '*_toeplitz_ls_results.csv'),
    ('irlsq_reg_{era}', ['small_kc1e10_dd', 'native_ill_dd', 'large_kc1e10_dd'],
     '*_irlsq_reg_results.csv'),
]


def load(pattern_dir, cell, glob_pat):
    """Return (rows, header_lines) for a cell, preferring the trsm_left arm."""
    for arm in ('trsm_left', 'merged', ''):
        hits = glob.glob(os.path.join(RESULTS, pattern_dir, cell, arm, glob_pat))
        if hits:
            path = sorted(hits)[-1]
            head = [l.strip() for l in open(path) if l.startswith('#')]
            rows = [r for r in csv.DictReader(
                l for l in open(path) if not l.startswith('#'))]
            return rows, head, path
    return None, None, None


def main():
    problems, notes, checked = [], [], 0

    for tmpl, cells, pat in FAMILIES:
        for cell in cells:
            new, head, path = load(tmpl.format(era='0831_a1'), cell, pat)
            if new is None:
                notes.append(f'{cell}: no 0831_a1 data yet (cell not pulled or still queued)')
                continue
            old, _, _ = load(tmpl.format(era='0829_a1'), cell, pat)
            if old is None:
                notes.append(f'{cell}: no 0829_a1 baseline to predict from; skipped')
                continue

            hdr = ' '.join(head)
            if EXPECTED_KNOB not in hdr:
                problems.append(f'{cell}: CSV header missing {EXPECTED_KNOB} '
                                f'(knob did not reach the run) -> {path}')
            if EXPECTED_COMMIT not in hdr:
                problems.append(f'{cell}: header commit is not {EXPECTED_COMMIT} -> {path}')

            base = {r['algorithm']: r for r in old if r['run'] == '0'}
            for r in (x for x in new if x['run'] == '0'):
                alg = r['algorithm']
                if alg not in base or alg.startswith('Blendenpik') or alg == 'unpreconditioned':
                    continue
                retries = int(base[alg].get('chol_retries', 0) or 0)
                failed = int(r['qr_status']) != 0
                checked += 1
                if (cell, alg) in BOUNDARY_EXEMPT:
                    notes.append(f'{cell}/{alg}: boundary-exempt, observed '
                                 f'{"FAIL" if failed else "success"} '
                                 f'(0829 retries={retries})')
                elif retries > 0 and not failed:
                    problems.append(f'{cell}/{alg}: predicted FAIL (0829 retries={retries}) '
                                    f'but qr_status=0; knob may not be in effect')
                elif retries == 0 and failed:
                    problems.append(f'{cell}/{alg}: predicted success (0829 retries=0) '
                                    f'but qr_status={r["qr_status"]}; regression')

    print(f'checked {checked} predicted rows')
    if notes:
        print('\nNOTES')
        for n in notes:
            print('  -', n)
    if problems:
        print('\nFAILURES')
        for p in problems:
            print('  !', p)
        return 1
    print('\nAll predictions held.' if checked else '\nNo rows checked yet.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
