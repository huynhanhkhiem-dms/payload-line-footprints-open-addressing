#!/usr/bin/env python3
"""Exact probe-order footprint envelope from Theorem 4.

Uses the recurrence
  h(q+1) = h(q) * (n-q-e)/(n-q),
where h(q) = C(n-q,e)/C(n,e), avoiding large integer binomials.
Standard library only.
"""
import csv
import sys

N = 2 ** 24
B = 8
DELTAS = [0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]


def envelope(n: int, e: int, b: int):
    if n % b:
        raise ValueError("n must be divisible by b")
    if not (1 <= e <= n):
        raise ValueError("e must lie in [1,n]")
    lines = n // b
    h = 1.0
    grouped = 0.0
    scatter = 0.0
    last_q = n - e
    for q in range(last_q + 1):
        if q < lines:
            scatter += h
        if q % b == 0 and q // b < lines:
            grouped += h
        if q < last_q:
            h *= (n - q - e) / (n - q)
    probes = (n + 1) / (e + 1)
    return grouped, scatter, probes


def main():
    output = sys.argv[1] if len(sys.argv) > 1 else "results_ordering_envelope.csv"
    rows = []
    for delta in DELTAS:
        e = int(delta * N)
        grouped, scatter, probes = envelope(N, e, B)
        rows.append({
            "n": N,
            "b": B,
            "delta": f"{delta:.3f}",
            "empty_slots": e,
            "grouped_line_lower": f"{grouped:.9f}",
            "scatter_first_upper": f"{scatter:.9f}",
            "same_expected_probes": f"{probes:.9f}",
            "upper_lower_ratio": f"{scatter/grouped:.9f}",
        })
    with open(output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
