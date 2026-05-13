#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Build a custom cisTarget database from MACS2 consensus peaks (per compartment)
# -----------------------------------------------------------------------------
# Adapted from the colleague's ArchR-based workflow
#   (scenicplus_from_ArchR_and_Seuratobj/shell_script/generate_customdb_fromarchrbedfile.sh)
# to run on MACS2 narrowPeak outputs produced upstream by Signac in P697.
#
# Rationale: scoring motifs only against peaks observed in our data (vs. the
# ~46 GB public hg38 screen v10 DB) is faster and more specific to the
# chromatin landscape of NK / Treg compartments.
#
# Usage:
#   ./build_customdb_macs2.sh <COMP>
# where <COMP> is one of: NK, Treg
#
# Prerequisites (install once):
#   * aertslab/create_cistarget_databases repo cloned to CREATE_CISTARGET_DB_DIR
#       git clone https://github.com/aertslab/create_cistarget_databases.git
#   * cluster-buster (cbust) compiled & on PATH
#       https://github.com/weng-lab/cluster-buster
#   * Aerts lab motif collection v10nr_clust_public/singletons (~27 GB)
#       https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/singletons/
#   * hg38 reference FASTA (matching chromosome naming of the MACS2 peaks, i.e. 'chr1' style)
#       https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#
# Notes:
#   - This is the "full v10nr_clust_public singletons" motif set (option A,
#     per user choice). To iterate faster later, you can filter MOTIF_LIST to
#     a curated subset before the scoring step.
#   - Padding (+/- 1 kb) matches the colleague's script. Tune via FLANK_BP.
# -----------------------------------------------------------------------------
set -euo pipefail

# ---- Args ----
COMP="${1:-}"
if [[ -z "${COMP}" || ( "${COMP}" != "NK" && "${COMP}" != "Treg" ) ]]; then
  echo "Usage: $0 <NK|Treg>" >&2
  exit 1
fi

# ---- Configurable paths (edit for your environment) ----
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
MACS_DIR="${PROJECT_DIR}/data/outputData/macs2_peaks"
OUT_DIR="${PROJECT_DIR}/data/inputData/cisTarget_dbs_custom"

# Prefer Coder shared scratch (`/scratch/$USER`, fast PVC-backed) for temp
# files; fall back to in-project tmp dir on laptops / non-Coder hosts.
if [[ -n "${USER:-}" && -d "/scratch/${USER}" && -w "/scratch/${USER}" ]]; then
  TMP_DIR="${TMP_DIR:-/scratch/${USER}/scenicplus_customdb_tmp_${COMP}}"
else
  TMP_DIR="${TMP_DIR:-${OUT_DIR}/tmp_${COMP}}"
fi

# External tooling (override via env vars as needed)
CREATE_CISTARGET_DB_DIR="${CREATE_CISTARGET_DB_DIR:-${HOME}/tools/create_cistarget_databases}"
GENOME_FA="${GENOME_FA:-${PROJECT_DIR}/data/inputData/hg38/hg38.fa}"
MOTIF_DIR="${MOTIF_DIR:-${PROJECT_DIR}/data/inputData/aertslab_motif_collection/v10nr_clust_public/singletons}"
MOTIF_ANNOT="${MOTIF_ANNOT:-${PROJECT_DIR}/data/inputData/aertslab_motif_collection/v10nr_clust_public/snapshots/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl}"

# Parameters
FLANK_BP="${FLANK_BP:-1000}"       # BED peak padding (each side)
# Default thread count: explicit env > OMP_NUM_THREADS > nproc > 8.
THREADS="${THREADS:-${OMP_NUM_THREADS:-$(nproc 2>/dev/null || echo 8)}}"
PREFIX="${COMP}_MACS2_peaks_v10nr_hg38"

# Self-documenting startup banner
echo "----------------------------------------------------------------------"
echo " Custom cisTarget DB build for compartment: ${COMP}"
echo "   PROJECT_DIR : ${PROJECT_DIR}"
echo "   TMP_DIR     : ${TMP_DIR}"
echo "   THREADS     : ${THREADS}"
echo "   GENOME_FA   : ${GENOME_FA}"
echo "   MOTIF_DIR   : ${MOTIF_DIR}"
echo "----------------------------------------------------------------------"

# ---- Inputs / outputs ----
PEAKS_NARROW="${MACS_DIR}/${COMP}_MACS2_peaks.narrowPeak"
PEAKS_BED="${TMP_DIR}/${COMP}_peaks.bed"
PADDED_FA="${TMP_DIR}/${PREFIX}.padded.fa"
MOTIF_LIST="${TMP_DIR}/motif_list.txt"

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# ---- 1. Sanity checks ----
echo "[1/5] Checking prerequisites..."
for f in "${PEAKS_NARROW}" "${GENOME_FA}"; do
  [[ -s "${f}" ]] || { echo "ERROR: missing/empty: ${f}" >&2; exit 2; }
done
for d in "${MOTIF_DIR}" "${CREATE_CISTARGET_DB_DIR}"; do
  [[ -d "${d}" ]] || { echo "ERROR: missing dir: ${d}" >&2; exit 2; }
done
command -v cbust >/dev/null || { echo "ERROR: cbust (cluster-buster) not on PATH" >&2; exit 2; }

# Genome FASTA index requirement: create_fasta_with_padded_bg_from_bed.sh needs
# a .fai next to the FASTA. If it's missing AND the FASTA dir isn't writable
# (common when GENOME_FA lives on a read-only NFS share), point the user at a
# scratch-based workaround instead of letting samtools fail cryptically.
if [[ ! -s "${GENOME_FA}.fai" ]]; then
  GENOME_DIR="$(dirname "${GENOME_FA}")"
  if [[ ! -w "${GENOME_DIR}" ]]; then
    cat >&2 <<EOF
ERROR: ${GENOME_FA}.fai is missing and ${GENOME_DIR} is not writable.
       Symlink the FASTA to a writable location and index there, e.g.:

         mkdir -p /scratch/\$USER/ref
         ln -sf "${GENOME_FA}" /scratch/\$USER/ref/hg38.fa
         samtools faidx /scratch/\$USER/ref/hg38.fa
         export GENOME_FA=/scratch/\$USER/ref/hg38.fa

       Then re-run this script.
EOF
    exit 2
  fi
  echo "      .fai not found; running 'samtools faidx ${GENOME_FA}'"
  command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 2; }
  samtools faidx "${GENOME_FA}"
fi

# ---- 2. Convert narrowPeak -> 3-col BED ----
echo "[2/5] Writing 3-col BED from ${PEAKS_NARROW}"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "${PEAKS_NARROW}" \
  | sort -k1,1 -k2,2n -u > "${PEAKS_BED}"
echo "      Peaks: $(wc -l < "${PEAKS_BED}")"

# ---- 3. Padded FASTA ----
echo "[3/5] Extracting padded FASTA (+/- ${FLANK_BP} bp)..."
bash "${CREATE_CISTARGET_DB_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
  "${GENOME_FA}" \
  "${GENOME_FA}.fai" \
  "${PEAKS_BED}" \
  "${PADDED_FA}" \
  "${FLANK_BP}" \
  "yes" \
  "1000"

# ---- 4. Motif list (full v10nr_clust_public singletons) ----
echo "[4/5] Enumerating motifs in ${MOTIF_DIR}"
find "${MOTIF_DIR}" -maxdepth 1 -name '*.cb' -printf '%f\n' \
  | sed 's/\.cb$//' > "${MOTIF_LIST}"
echo "      Motifs: $(wc -l < "${MOTIF_LIST}")"

# ---- 5. Score motifs vs regions, then rank ----
echo "[5/5] Scoring motifs against regions (this is the slow step)..."
python "${CREATE_CISTARGET_DB_DIR}/create_cistarget_motif_databases.py" \
  -f "${PADDED_FA}" \
  -M "${MOTIF_DIR}" \
  -m "${MOTIF_LIST}" \
  -o "${OUT_DIR}/${PREFIX}" \
  -t "${THREADS}"

SCORES_FEATHER="${OUT_DIR}/${PREFIX}.regions_vs_motifs.scores.feather"
RANKS_FEATHER="${OUT_DIR}/${PREFIX}.regions_vs_motifs.rankings.feather"

[[ -s "${SCORES_FEATHER}" ]] || { echo "ERROR: scores feather not produced" >&2; exit 3; }

echo "      Converting scores -> rankings..."
python "${CREATE_CISTARGET_DB_DIR}/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py" \
  -i "${SCORES_FEATHER}" \
  -o "${RANKS_FEATHER}" \
  -s 555

echo ""
echo "DONE. Point the SCENIC+ snakemake config at:"
echo "  ctx_db_fname : ${RANKS_FEATHER}"
echo "  dem_db_fname : ${SCORES_FEATHER}"
echo "  path_to_motif_annotations : ${MOTIF_ANNOT}"
