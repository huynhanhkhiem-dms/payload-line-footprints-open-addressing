#!/bin/sh
set -eu

CC=${CC:-gcc}
CFLAGS=${CFLAGS:--O3 -std=c11 -Wall -Wextra}
FAST=${FAST:-0}

$CC $CFLAGS -DINSTRUMENT experiment.c -lm -o experiment_instrumented

HEADER='scheme,lg2n,delta,bsize,ehC,ins_probes,ins_lines,ins_pages,ins_p99,ins_maxp,ins_maxl,pos_probes,pos_lines,pos_p99,pos_maxp,neg_probes,neg_lines,neg_p99,neg_maxp,t_build,t_pos,t_neg,m,nq'
THEORY_HEADER='lg2n,delta_target,empty_slots,ideal_probe_expectation,ideal_line_expectation_b8,observed_neg_probes,observed_neg_lines,nq'
CHECK_HEADER='scheme,lg2n,delta,fails,missing,falsepos,spills,ins_maxl,pos_maxl,neg_maxl,fh_line_bound,fh_state_miss_bound'
LARGE_HEADER='seed,scheme,lg2n,delta,pos_probes,pos_lines,pos_maxp,neg_probes,neg_lines,neg_maxp,nq,fails,missing,falsepos,spills,pos_maxl,neg_maxl,fh_line_bound'

if [ "$FAST" = 1 ]; then
  # Deliberately small integrity configuration: compiles every path, checks
  # correctness/bounds, exercises the line-reference calculation, and finishes
  # quickly on ordinary hardware.  It is not used for manuscript numbers.
  LG=17
  MAIN_DELTAS='0.01 0.001'
  MAIN_NQ=40
  ROB_SEEDS='1'
  ROB_NQ=30
  LARGE_ROB_LG=17
  LARGE_ROB_SEEDS='1'
  LARGE_ROB_NQ=30
  SCALE_LGS='17 18'
  SENS='0.01:2 0.001:2'
  SENS_NQ=30
  ALIGN_NQ=30
else
  LG=24
  MAIN_DELTAS='0.10 0.05 0.02 0.01 0.005 0.002 0.001'
  MAIN_NQ=5000
  ROB_SEEDS='1 2 3 4 5'
  ROB_NQ=2000
  LARGE_ROB_LG=24
  LARGE_ROB_SEEDS='1 2 3 4 5'
  LARGE_ROB_NQ=2000
  SCALE_LGS='20 22 24'
  SENS='0.01:1 0.01:2 0.01:4 0.01:8 0.001:2 0.001:4 0.001:8'
  SENS_NQ=1000
  ALIGN_NQ=5000
fi

rm -f reproduced_*.csv
printf '%s\n' "$CHECK_HEADER" > reproduced_checks.csv

emit_check() {
  printf '%s\n' "$1" | awk '/^CHECKCSV,/{sub(/^CHECKCSV,/,""); print}' >> reproduced_checks.csv
}

row_for() {
  scheme=$1; nq=$2; lg=$3; delta=$4; c=${5:-2}; seed=${6:-0}
  out=$(NQ="$nq" SEED="$seed" ./experiment_instrumented "$scheme" "$lg" "$delta" 8 "$c")
  emit_check "$out"
  printf '%s\n' "$out" | awk '/^CSV,/{sub(/^CSV,/,""); print}'
}

# 1) Primary sweep, RP exact-model check, and WTL-DH realized layouts.
printf '%s\n' "$HEADER" > reproduced_main.csv
printf '%s\n' "$THEORY_HEADER" > reproduced_theory.csv
printf '%s\n' 'scheme,lg2n,delta,empty_slots,nonfull_lines,empty_lines,max_empty_per_line,sum_sq_empty' > reproduced_layout_wtl.csv
for d in $MAIN_DELTAS; do
  for s in up dh wtl eh ehb fh; do
    out=$(NQ="$MAIN_NQ" SEED=0 ./experiment_instrumented "$s" "$LG" "$d" 8 2)
    emit_check "$out"
    row=$(printf '%s\n' "$out" | awk '/^CSV,/{sub(/^CSV,/,""); print}')
    printf '%s\n' "$row" >> reproduced_main.csv
    if [ "$s" = up ]; then
      th=$(printf '%s\n' "$out" | awk -F'[=,]' '/^UPTHEORY/{print $5 "," $9 "," $11}')
      empty=$(printf '%s' "$th" | cut -d, -f1)
      pexp=$(printf '%s' "$th" | cut -d, -f2)
      lexp=$(printf '%s' "$th" | cut -d, -f3)
      negp=$(printf '%s' "$row" | cut -d, -f16)
      negl=$(printf '%s' "$row" | cut -d, -f17)
      nq=$(printf '%s' "$row" | cut -d, -f24)
      printf '%s,%s,%s,%s,%s,%s,%s,%s\n' "$LG" "$d" "$empty" "$pexp" "$lexp" "$negp" "$negl" "$nq" >> reproduced_theory.csv
    fi
    if [ "$s" = wtl ]; then
      printf '%s\n' "$out" | awk '/^LAYOUTCSV,/{sub(/^LAYOUTCSV,/,""); print}' >> reproduced_layout_wtl.csv
    fi
  done
done
python3 line_profile_model.py reproduced_layout_wtl.csv reproduced_main.csv reproduced_line_reference.csv
python3 affine_gap_model.py reproduced_layout_wtl.csv reproduced_affine_gap_bounds.csv
python3 ordering_envelope.py reproduced_ordering_envelope.csv
python3 line_width_envelope.py reproduced_line_width_envelope.csv

# 2) FH headline fill profile and both deterministic bounds.
printf '%s\n' 'level,size_slots,occupied_slots,nonempty_buckets,full_buckets' > reproduced_fh_profile.csv
fhout=$(FH_PROFILE=1 NQ="$MAIN_NQ" SEED=0 ./experiment_instrumented fh "$LG" 0.001 8 2)
emit_check "$fhout"
printf '%s\n' "$fhout" | awk '/^FHPROFILE,/{sub(/^FHPROFILE,/,""); print}' >> reproduced_fh_profile.csv
printf '%s\n' 'scheme,lg2n,delta,fails,missing,falsepos,spills,ins_maxl,pos_maxl,neg_maxl,general_bound,state_miss_bound,nq' > reproduced_funnel_bound.csv
chk=$(printf '%s\n' "$fhout" | awk '/^CHECKCSV,/{sub(/^CHECKCSV,/,""); print}')
nq=$(printf '%s\n' "$fhout" | awk -F, '/^CSV,/{print $25}')
printf '%s,%s\n' "$chk" "$nq" >> reproduced_funnel_bound.csv

# 3) Matched WTL-DH robustness at the headline point.
printf 'seed,%s\n' "$HEADER" > reproduced_wtl_robustness.csv
for seed in $ROB_SEEDS; do
  out=$(NQ="$ROB_NQ" SEED="$seed" ./experiment_instrumented wtl "$LG" 0.001 8 2)
  emit_check "$out"
  row=$(printf '%s\n' "$out" | awk '/^CSV,/{sub(/^CSV,/,""); print}')
  printf '%s,%s\n' "$seed" "$row" >> reproduced_wtl_robustness.csv
done

# 4) Independent large-table RP/EHB/FH robustness used for reported means.
printf '%s\n' "$LARGE_HEADER" > reproduced_large_robustness.csv
for seed in $LARGE_ROB_SEEDS; do
  for s in up ehb fh; do
    out=$(QUERY_ONLY_INSTRUMENT=1 NQ="$LARGE_ROB_NQ" SEED="$seed" ./experiment_instrumented "$s" "$LARGE_ROB_LG" 0.001 8 2)
    emit_check "$out"
    row=$(printf '%s\n' "$out" | awk -F, '/^CSV,/{print $2","$3","$4","$13","$14","$16","$17","$18","$20","$25}')
    chk=$(printf '%s\n' "$out" | awk -F, '/^CHECKCSV,/{print $5","$6","$7","$8","$10","$11","$12}')
    printf '%s,%s,%s\n' "$seed" "$row" "$chk" >> reproduced_large_robustness.csv
  done
done

# 5) Elastic finite-constant sensitivity.
printf '%s\n' "$HEADER" > reproduced_sensitivity.csv
for x in $SENS; do
  d=${x%:*}; c=${x#*:}
  row_for eh "$SENS_NQ" "$LG" "$d" "$c" 0 >> reproduced_sensitivity.csv
done

# 6) Table-size scaling at delta=0.001.
printf '%s\n' "$HEADER" > reproduced_scaling.csv
for lg in $SCALE_LGS; do
  for s in up eh ehb fh; do
    row_for "$s" "$ROB_NQ" "$lg" 0.001 2 0 >> reproduced_scaling.csv
  done
done

# 7) The theorem-aligned delta=1/1024 point.
printf '%s\n' "$HEADER" > reproduced_theorem_aligned.csv
for s in up dh eh ehb fh; do
  row_for "$s" "$ALIGN_NQ" "$LG" 0.0009765625 2 0 >> reproduced_theorem_aligned.csv
done

# Hard integrity gate across every executed run.
if awk -F, 'NR>1 && ($4+$5+$6+$7)>0 {bad=1} END{exit bad}' reproduced_checks.csv; then :; else
  echo 'fatal correctness counter' >&2
  exit 2
fi

printf '%s\n' 'Reproduction complete: all requested studies generated and all correctness/bound checks passed.'
