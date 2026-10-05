#!/usr/bin/env python3
"""Compute the fixed-layout uniform-line-permutation reference.

This script does NOT claim that the WTL-DH affine line schedule is a uniform
random permutation of all L lines.  The reference (L+1)/(K+1) is exact only
under the model stated in Proposition 2 of the manuscript: conditional on a
fixed layout with K nonfull lines, query lines are visited in an independent
uniformly random permutation and each visited line is exhausted before moving
on.  We report the WTL-DH observation next to this reference only as an
empirical comparison.
"""
import csv
import sys

LAYOUT = sys.argv[1] if len(sys.argv) > 1 else 'results_layout_wtl.csv'
MAIN = sys.argv[2] if len(sys.argv) > 2 else 'results_main.csv'
OUT = sys.argv[3] if len(sys.argv) > 3 else 'results_line_reference.csv'

with open(MAIN, newline='', encoding='utf-8') as f:
    wtl = {round(float(r['delta']), 9): r for r in csv.DictReader(f) if r['scheme'] == 'wtl'}

rows = []
with open(LAYOUT, newline='', encoding='utf-8') as f:
    for r in csv.DictReader(f):
        d = round(float(r['delta']), 9)
        if d not in wtl:
            continue
        n = 2 ** int(r['lg2n'])
        b = 8
        L = (n + b - 1) // b
        K = int(r['nonfull_lines'])
        e = int(r['empty_slots'])

        # Expected number of nonfull lines if the e empty slots were instead a
        # uniformly random e-subset of the n slots.  This is a layout reference.
        p_full = 1.0
        for i in range(b):
            p_full *= (n - e - i) / (n - i)
        uniform_nonfull = L * (1.0 - p_full)

        # Exact Proposition-2 reference under an independent uniformly random
        # permutation of the L lines.  WTL-DH uses a structured affine family,
        # so this value is NOT a theorem-derived WTL-DH prediction.
        reference = (L + 1) / (K + 1)
        observed = float(wtl[d]['neg_lines'])
        rows.append({
            'delta': f'{d:.6f}',
            'nlines': L,
            'nonfull': K,
            'uniform_expected_nonfull': f'{uniform_nonfull:.6f}',
            'uniform_permutation_reference_lines': f'{reference:.9f}',
            'observed_wtl_lines': f'{observed:.4f}',
            'wtl_minus_reference_pct': f'{100 * (observed - reference) / reference:.6f}',
        })

if not rows:
    raise SystemExit('no matching WTL rows')

with open(OUT, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
