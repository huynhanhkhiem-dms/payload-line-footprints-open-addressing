#!/usr/bin/env python3
"""Exact line-width sensitivity for the permutation-order envelope in Theorem 4.

The calculation uses the same stable binomial-ratio recurrence as
ordering_envelope.py.  It changes only the slots-per-line parameter b.
"""
import csv
import sys
from ordering_envelope import envelope

N = 2 ** 24
DELTA = 0.001
WIDTHS = (4, 8, 16)


def main():
    output = sys.argv[1] if len(sys.argv) > 1 else "results_line_width_envelope.csv"
    e = int(DELTA * N)
    rows = []
    for b in WIDTHS:
        lower, upper, probes = envelope(N, e, b)
        rows.append({
            "n": N,
            "delta": f"{DELTA:.3f}",
            "empty_slots": e,
            "slots_per_line": b,
            "grouped_line_lower": f"{lower:.9f}",
            "scatter_first_upper": f"{upper:.9f}",
            "same_expected_probes": f"{probes:.9f}",
            "upper_lower_ratio": f"{upper/lower:.9f}",
        })
    with open(output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
