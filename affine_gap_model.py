#!/usr/bin/env python3
"""Sharp K-only bounds for a fixed full-cycle line order.

For L payload lines with K nonfull lines and F=L-K full lines, let f_1,...,f_K
be the numbers of consecutive full lines immediately preceding the K nonfull
lines in the cyclic order induced by a fixed full-cycle step.  Under a uniform
starting line, the exact expected number of visited lines is

    (1/L) * sum_i (f_i+1)(f_i+2)/2.

Given only L and K, convexity gives a sharp minimum when the F full lines are
as evenly distributed as possible among the K gaps, and a sharp maximum when
all F full lines lie in one gap.  This script evaluates those bounds for the
realized WTL-DH layouts reported by experiment.c.
"""
import csv
import sys

LAYOUT = sys.argv[1] if len(sys.argv) > 1 else 'results_layout_wtl.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'results_affine_gap_bounds.csv'

rows=[]
with open(LAYOUT,newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f):
        L=(2**int(r['lg2n']))//8
        K=int(r['nonfull_lines'])
        if K <= 0:
            continue
        F=L-K
        q,rem=divmod(F,K)
        # gap f=q contributes (q+1)(q+2)/2; gap f=q+1 contributes
        # (q+2)(q+3)/2.
        min_num=(K-rem)*(q+1)*(q+2)/2 + rem*(q+2)*(q+3)/2
        max_num=(F+1)*(F+2)/2 + (K-1)
        rows.append({
            'delta': r['delta'],
            'nlines': L,
            'nonfull': K,
            'full': F,
            'balanced_gap_q': q,
            'balanced_gap_remainder': rem,
            'sharp_min_uniform_start_lines': f'{min_num/L:.9f}',
            'sharp_max_uniform_start_lines': f'{max_num/L:.9f}',
            'uniform_permutation_reference_lines': f'{(L+1)/(K+1):.9f}',
        })

with open(OUT,'w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
