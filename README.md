# Payload-Line Footprints in Open Addressing

Reproducibility artifact for the manuscript **“Payload-Line Footprints in Open Addressing: Probe-Order Bounds and Cyclic-Gap Identities.”**

**Author:** Huynh Anh Khiem  
**Affiliation:** Faculty of Information Technology, Ton Duc Thang University, Ho Chi Minh City, Vietnam  
**Email:** huynhanhkhiem@tdtu.edu.vn  
**ORCID:** 0009-0007-7210-174X

## Purpose

This repository reproduces the finite table-array payload-line results reported in the manuscript. A payload line contains eight 64-bit table slots (64 bytes). Payload-line counts are address-footprint measurements; they are **not** hardware cache-miss measurements.

## Requirements

- POSIX shell
- C11 compiler (the reproduction script defaults to `gcc`)
- Python 3
- Matplotlib only if PDF figures are regenerated

No external dataset or network access is required.

## Quick integrity check

```bash
FAST=1 sh run.sh
```

`FAST=1` uses deliberately small tables and query counts to exercise all code paths and checks. It is **not** used for manuscript numbers.

## Full reproduction

```bash
sh run.sh
```

The script compiles `experiment.c` with `-O3 -std=c11 -Wall -Wextra -DINSTRUMENT`. Any insertion failure, missing present key, false positive, Elastic spill, Funnel geometry-bound violation, or Funnel state-aware miss-bound violation is fatal.

## Implemented schemes

| Code | Scheme |
|---|---|
| `up` | Random-slot probing (RP), with replacement |
| `dh` | Full-cycle double hashing (DH) |
| `wtl` | Walking-the-line double hashing (WTL-DH): scan all eight slots of one payload line, then advance by an odd double-hash step over payload lines; insertion and lookup use the same schedule |
| `eh` | Finite Elastic Hashing realization |
| `ehb` | Experimental blocked Elastic variant |
| `fh` | Finite Funnel Hashing realization |

## Model distinction

The sharp permutation envelope assumes a uniformly random empty-slot set independent of a fixed probe permutation. It is a conditional theorem. In the large-table limit at free fraction `delta`, the upper/lower ratio is `[1-(1-delta)^b]/delta` and tends to the line width `b` as `delta` tends to zero. A self-consistent WTL-DH build reshapes the empty-line layout and need not preserve either the conditional probe count or the conditional lower-envelope value.

For a fixed table with `L` payload lines and `K` nonfull lines, Proposition 2 gives `(L+1)/(K+1)` only when query lines are visited in an independent uniformly random permutation and each visited line is exhausted. WTL-DH does not sample uniformly from all `L!` line permutations: it uses a structured full-cycle affine family `a + j*d (mod L)` with an odd step `d`. Proposition 3 treats that structure directly for a fixed step: if `f_1,...,f_K` are the counts of full lines in the cyclic gaps immediately preceding the nonfull lines, a uniform starting line has exact expected footprint `(1/L) sum_i (f_i+1)(f_i+2)/2`. This shows why `K` alone is insufficient.

`line_profile_model.py` computes the uniform-permutation reference. `affine_gap_model.py` computes the sharp lower and upper bounds implied by `L` and `K` alone for any fixed full-cycle cyclic order.

## Funnel state certificate

With `FH_PROFILE=1`, `experiment.c` emits the main-level fill profile. At `n=2^24` and `delta=0.001`, levels 0–25 have full buckets, level 26 still has full buckets, and level 27 has none. Each 20-slot main bucket spans at most three eight-slot payload lines in the realized alignment. Therefore every miss terminates by level 27 and the realized state-aware bound is `28 x 3 = 84` lines. The executable computes and enforces this bound in addition to the geometry-only bound 211.

## Supplied result files

| File | Content |
|---|---|
| `results_main.csv` | Primary RP/DH/WTL-DH/EH/EHB/FH sweep |
| `results_theory.csv` | Exact RP model check |
| `results_ordering_envelope.csv` | Conditional probe-order envelope |
| `results_layout_wtl.csv` | Realized WTL-DH empty-line layout |
| `results_line_reference.csv` | Uniform-line-permutation reference and WTL-DH observation |
| `results_affine_gap_bounds.csv` | Sharp `L,K` bounds for fixed full-cycle line orders |
| `results_fh_profile.csv` | Headline FH fill profile |
| `results_funnel_bound.csv` | Headline correctness and both FH bounds |
| `results_wtl_robustness.csv` | Five independent WTL-DH headline builds |
| `results_large_robustness.csv` | Five independent RP/EHB/FH headline builds |
| `results_sensitivity.csv` | Elastic finite-constant sensitivity |
| `results_scaling.csv` | Table-size sweep |
| `results_theorem_aligned.csv` | `delta=1/1024` check |
| `results_line_width_envelope.csv` | Exact Theorem 4 line-width sensitivity for `b = 4, 8, 16` |

## Derived-analysis scripts

```text
ordering_envelope.py [output.csv]
line_width_envelope.py [output.csv]
line_profile_model.py [layout.csv] [main.csv] [output.csv]
affine_gap_model.py [layout.csv] [output.csv]
make_figures.py [main.csv]
```

## Provenance

The C executable emits primary rows, correctness counters, WTL-DH layout summaries, FH level profiles, and deterministic FH bounds. The numerical Python scripts contain the documented derived calculations. No hidden numerical post-processing is needed for the manuscript claims.

The common CSV format retains elapsed-time fields emitted by the executable. The manuscript does not use them as scientific evidence because the validation host is shared and the artifact does not provide a controlled multi-processor hardware-counter study.

## Headline checks

For `n=2^24`, `delta=0.001`, `seed=0`:

| Quantity | Value |
|---|---:|
| RP miss lines | 998.7196 |
| WTL-DH miss probes | 1514.9214 |
| WTL-DH miss lines | 189.8706 |
| WTL-DH nonfull payload lines | 10943 |
| Uniform-line-permutation reference from `K` | 191.625822368 lines |
| Sharp fixed-cycle lower bound from `L,K` | 96.322219849 lines |
| FH miss probes | 547.7310 |
| FH miss lines | 82.6556 |
| FH observed miss maximum | 84 |
| FH geometry-only bound | 211 |
| FH state-aware miss bound | 84 |
| Correctness counters | all zero |

Conditional permutation-envelope values at the same headline empty count are approximately `125.432323037` and `999.953331744` lines, with expected first-empty position `999.953331744`. They are model quantities, not a self-consistent WTL-DH performance ratio.

## Line-width sensitivity

`line_width_envelope.py` generates `results_line_width_envelope.csv` for `b = 4, 8, 16` at `n = 2^24` and `delta = 0.001` using the exact Theorem 4 recurrence. The reproduction script generates the matching `reproduced_line_width_envelope.csv`.

## Reproducibility note

A local integrity run of

```bash
FAST=1 sh run.sh
```

completed successfully before this repository was populated, with all requested studies generated and all correctness/bound checks passing.
