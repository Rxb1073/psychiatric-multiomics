#!/bin/bash
ROOT=/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq
declare -A EXP=( [PRJNA319583]=348 [PRJNA314463]=82 [PRJNA574470]=68 )

for D in PRJNA319583 PRJNA314463 PRJNA574470; do
  echo "================ $D ================"
  D2=$ROOT/$D

  # 1. count completed PE pairs (both mates present)
  pairs=$(ls $D2/trimmed/*_val_1.fq.gz 2>/dev/null | wc -l)
  mates=$(ls $D2/trimmed/*_val_2.fq.gz 2>/dev/null | wc -l)
  echo "PE pairs: $pairs/_val_1  $mates/_val_2  (expect ${EXP[$D]})"

  # 2. SE (574470 only)
  if [ "$D" = "PRJNA574470" ]; then
    se=$(ls $D2/trimmed/*_trimmed.fq.gz 2>/dev/null | wc -l)
    echo "SE trimmed: $se  (expect 10)"
  fi

  # 3. empty/truncated outputs (precaution: size check, not just existence)
  tiny=$(find $D2/trimmed -name '*.fq.gz' -size -1M 2>/dev/null | wc -l)
  echo "outputs <1MB (suspicious): $tiny"
  [ "$tiny" -gt 0 ] && find $D2/trimmed -name '*.fq.gz' -size -1M

  # 4. scan logs for failures / non-zero exits
  fails=$(grep -liE 'error|exception|failed|truncated|No such' $D2/logs/trim_*.err 2>/dev/null | wc -l)
  echo "logs with error keywords: $fails"
  [ "$fails" -gt 0 ] && grep -liE 'error|exception|failed|truncated' $D2/logs/trim_*.err 2>/dev/null | head -5

  # 5. confirm each completed task printed its 'done' line
  done_lines=$(grep -l "done:" $D2/logs/trim_*.out 2>/dev/null | wc -l)
  echo "tasks that printed 'done': $done_lines"
  echo ""
done
