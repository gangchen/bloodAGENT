#!/usr/bin/env bash
# Swap data/source/v16/derived/*.dat into data/config/*. Each replaced file is
# backed up alongside its original as <name>.pre-v16.bak (idempotent — re-running
# this script is safe; the backup is only written on the first run).
#
# Usage:
#   data/source/v16/tools/apply_to_config.sh [--apply]
#
# Without --apply this prints the planned actions but does nothing. With --apply
# it performs the swap.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${HERE}/../../../.."  # ..../bloodAGENT
DERIVED="${HERE}/../derived"
CONFIG="${ROOT}/data/config"

APPLY=0
[[ "${1:-}" == "--apply" ]] && APPLY=1

backup_then_install() {
  local src=$1 dst=$2
  if [[ ! -f "$src" ]]; then
    echo "  SKIP (no source): $dst (expected $src)"
    return
  fi
  if [[ -f "$dst" && ! -f "${dst}.pre-v16.bak" ]]; then
    echo "  BAK  $dst -> ${dst}.pre-v16.bak"
    [[ $APPLY -eq 1 ]] && cp "$dst" "${dst}.pre-v16.bak"
  fi
  echo "  COPY $src -> $dst"
  [[ $APPLY -eq 1 ]] && cp "$src" "$dst"
}

append_supplement() {
  local sup=$1 tgt=$2
  if [[ ! -f "$sup" ]]; then
    echo "  SKIP (no supplement): $tgt (expected $sup)"
    return
  fi
  if [[ -f "$tgt" && ! -f "${tgt}.pre-v16.bak" ]]; then
    echo "  BAK  $tgt -> ${tgt}.pre-v16.bak"
    [[ $APPLY -eq 1 ]] && cp "$tgt" "${tgt}.pre-v16.bak"
  fi
  echo "  APPEND $sup -> $tgt"
  [[ $APPLY -eq 1 ]] && cat "$sup" >> "$tgt"
}

echo "=== Master variation_annotation.dat ==="
backup_then_install "${DERIVED}/variation_annotation.v16.dat" "${CONFIG}/variation_annotation.dat"

echo
echo "=== Per-pipeline gt2pt ==="
backup_then_install "${DERIVED}/pipeline/CMR/genotype_to_phenotype_annotation_ICACMR.v16.dat"        "${CONFIG}/CMR/genotype_to_phenotype_annotation_ICACMR.dat"
backup_then_install "${DERIVED}/pipeline/Dragen/genotype_to_phenotype_annotation_Dragen.v16.dat"     "${CONFIG}/Dragen/genotype_to_phenotype_annotation_Dragen.dat"
backup_then_install "${DERIVED}/pipeline/HGDP/genotype_to_phenotype_annotation_HGDP.v16.dat"         "${CONFIG}/HGDP/genotype_to_phenotype_annotation_HGDP.dat"
backup_then_install "${DERIVED}/pipeline/Microarray/genotype_to_phenotype_annotation_Array.v16.dat"  "${CONFIG}/Microarray/genotype_to_phenotype_annotation_Array.dat"
backup_then_install "${DERIVED}/pipeline/ONT/genotype_to_phenotype_annotation_MINIMAP2SNIFFLES.v16.dat" "${CONFIG}/ONT/genotype_to_phenotype_annotation_MINIMAP2SNIFFLES.dat"
backup_then_install "${DERIVED}/pipeline/PacBio/genotype_to_phenotype_annotation_TGSGATK.v16.dat"    "${CONFIG}/PacBio/genotype_to_phenotype_annotation_TGSGATK.dat"
backup_then_install "${DERIVED}/pipeline/PacBio/genotype_to_phenotype_annotation_TGSPBSV.v16.dat"    "${CONFIG}/PacBio/genotype_to_phenotype_annotation_TGSPBSV.dat"

echo
echo "=== exonic_annotation supplement (12 new systems appended) ==="
append_supplement "${DERIVED}/exonic_annotation.hg19.BGStarget.supplement.txt" "${CONFIG}/exonic_annotation.hg19.BGStarget.txt"
append_supplement "${DERIVED}/exonic_annotation.hg38.BGStarget.supplement.txt" "${CONFIG}/exonic_annotation.hg38.BGStarget.txt"

echo
if [[ $APPLY -eq 1 ]]; then
  echo "DONE. Backups are at <file>.pre-v16.bak — to undo:"
  echo "  for f in \$(find data/config -name '*.pre-v16.bak'); do mv \"\$f\" \"\${f%.pre-v16.bak}\"; done"
else
  echo "DRY RUN. Re-run with --apply to actually swap files."
fi
