#!/bin/bash
set -u
cd /rds/projects/e/elhamsak-dissertation-autism/RNA_seq/_pheno_staging || exit 1

show_matrix () {
  local GSE=$1 SUB=$2 F
  F=${GSE}_series_matrix.txt.gz
  echo "=============== $GSE ==============="
  if [ ! -s "$F" ]; then
    local BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/${SUB}/${GSE}/matrix/"
    echo "  not local; directory listing:"
    curl -s "$BASE" | grep -oE 'GSE[0-9A-Za-z_-]+_series_matrix\.txt\.gz' | sort -u | sed 's/^/    /'
    wget -q "${BASE}${GSE}_series_matrix.txt.gz" -O "$F" 2>/dev/null
  fi
  [ -s "$F" ] || { echo "  no plain-named matrix (use a listed filename above)"; return; }

  echo "  samples (GSM): $(zcat "$F" | grep -m1 '^!Sample_geo_accession' | tr '\t' '\n' | grep -c GSM)"
  echo "  -- characteristic fields -> distinct values --"
  zcat "$F" | grep '^!Sample_characteristics_ch1' | while IFS= read -r line; do
    tag=$(printf '%s' "$line" | cut -f2 | sed 's/^"//' | grep -oE '^[^:]+:' | sed 's/:$//')
    vals=$(printf '%s' "$line" | cut -f2- | tr '\t' '\n' \
           | sed 's/^"//; s/"$//; s/^[^:]*: *//' | sort -u | head -6 | paste -sd'|')
    printf '     %-22s -> %s\n' "${tag:-<none>}" "$vals"
  done
  echo "  -- one SRA relation line (GSM->SRX bridge) --"
  zcat "$F" | grep -i '^!Sample_relation' | head -1 | cut -f2 | sed 's/^/     /'
}

show_matrix GSE102741 GSE102nnn
show_matrix GSE211154 GSE211nnn
show_matrix GSE59288  GSE59nnn
show_matrix GSE51264  GSE51nnn
