#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Run pycisTopic MALLET LDA topic modeling in the terminal (NK compartment).
#
# Why terminal (not notebook)?
#   LDA with MALLET on a large peak x cell matrix is memory-intensive and often
#   crashes Jupyter kernels mid-convergence. Running in a standalone terminal
#   with nohup gives us resumable logs, clean process isolation, and avoids
#   silent kernel crashes.
#
# Hyperparameters adopted from the colleague's
#   scenicplus_from_ArchR_and_Seuratobj/shell_script/run_all_topics_CLI.sh
# with the expanded topic sweep (up to 60) to let evaluate_models pick from
# a richer range than the notebook's [2..30].
#
# Usage (from an interactive shell where the scenicplus conda env is active):
#   conda activate scenicplus_env
#   nohup bash code/scenicplus_lda/run_lda_NK.sh > scenicplus_NK_lda.log 2>&1 &
#   tail -f scenicplus_NK_lda.log
# -----------------------------------------------------------------------------
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COMP="NK"
OUT_DIR="${PROJECT_DIR}/data/outputData/scenicplus_${COMP}"

# Prefer Coder shared scratch; fall back to in-project tmp on laptops.
if [[ -n "${USER:-}" && -d "/scratch/${USER}" && -w "/scratch/${USER}" ]]; then
  TMP_DIR="${TMP_DIR:-/scratch/${USER}/mallet_tmp_${COMP}}"
else
  TMP_DIR="${TMP_DIR:-${OUT_DIR}/mallet_tmp}"
fi

# Default thread count: explicit env > OMP_NUM_THREADS > nproc > 8.
THREADS="${THREADS:-${OMP_NUM_THREADS:-$(nproc 2>/dev/null || echo 8)}}"

INPUT_PKL="${OUT_DIR}/${COMP}_cistopic_obj_preLDA.pkl"
OUTPUT_PKL="${OUT_DIR}/models_all_topics.pkl"
MALLET_BIN="${MALLET_BIN:-${PROJECT_DIR}/data/inputData/Mallet-202108/bin/mallet}"

echo "----------------------------------------------------------------------"
echo " MALLET LDA topic modeling for compartment: ${COMP}"
echo "   TMP_DIR    : ${TMP_DIR}"
echo "   THREADS    : ${THREADS}"
echo "   INPUT_PKL  : ${INPUT_PKL}"
echo "   OUTPUT_PKL : ${OUTPUT_PKL}"
echo "----------------------------------------------------------------------"

mkdir -p "${TMP_DIR}"

[[ -s "${INPUT_PKL}" ]] || {
  echo "ERROR: input cistopic pickle missing: ${INPUT_PKL}"
  echo "Run the pre-LDA serialization cell in P697_scenicplus_NK.ipynb first."
  exit 2
}
[[ -x "${MALLET_BIN}" ]] || {
  echo "ERROR: mallet binary not found/executable: ${MALLET_BIN}"
  exit 2
}

# NOTE: tune -p (CPUs) for local machine; colleague used 32.
# Topic sweep [2, 10, 20, 30, 40, 50, 60] matches colleague; -n 150 convergence
# iterations and -m 400 MALLET internal iterations also match.
pycistopic topic_modeling mallet \
  -i "${INPUT_PKL}" \
  -o "${OUTPUT_PKL}" \
  -t 2 10 20 30 40 50 60 \
  -p "${THREADS}" \
  -n 150 \
  -a 50 \
  -A \
  -e 0.1 \
  -m 400 \
  -s 555 \
  -k \
  -T "${TMP_DIR}" \
  -b "${MALLET_BIN}"

echo "DONE. Topic models written to: ${OUTPUT_PKL}"
