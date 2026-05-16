#!/usr/bin/env bash
#SBATCH --job-name=phyloark_quota
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=phyloark_quota_%j.out
#SBATCH --error=phyloark_quota_%j.err

# Run from the directory where you submitted the job.
cd "${SLURM_SUBMIT_DIR}"

# -----------------------------
# User settings
# -----------------------------

SCRIPT="taxonomy_quota_sample.py"
GENOMES="enterobacter_merged_report.tsv"
TAXDUMP="taxonomy"
TARGET=1000
ALPHA=0.5
OUT_PREFIX="phyloark_sampling/phyloark_stability_alpha05_1000"

# Recommended: set this to the taxid representing your focal root.
# Example possibilities, depending on your intended scope:
#   Enterobacteriaceae family: 543
#   Enterobacterales order: 91347
# Leave blank to let the script infer the LCA of all input genome taxids.
ROOT_TAXID=""

# Optional: restrict internal nodes to major ranks only.
# Leave blank to preserve all induced NCBI taxonomy nodes.
# Example: RANK_KEEP="family,genus,species"
RANK_KEEP=""

# -----------------------------
# Environment
# -----------------------------

module load python

mkdir -p "$(dirname "${OUT_PREFIX}")"

# -----------------------------
# Sanity checks
# -----------------------------

if [[ ! -f "${SCRIPT}" ]]; then
  echo "ERROR: Could not find script: ${SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${GENOMES}" ]]; then
  echo "ERROR: Could not find genome table: ${GENOMES}" >&2
  exit 1
fi

if [[ ! -f "${TAXDUMP}/nodes.dmp" ]]; then
  echo "ERROR: Could not find required taxdump file: ${TAXDUMP}/nodes.dmp" >&2
  exit 1
fi

if [[ ! -f "${TAXDUMP}/names.dmp" ]]; then
  echo "WARNING: Could not find ${TAXDUMP}/names.dmp; output names will be blank." >&2
fi

# -----------------------------
# Build command
# -----------------------------

CMD=(
  python3 "${SCRIPT}"
  --genomes "${GENOMES}"
  --taxdump "${TAXDUMP}"
  --target "${TARGET}"
  --alpha "${ALPHA}"
  --out-prefix "${OUT_PREFIX}"
  --include-lineage-names
)

if [[ -n "${ROOT_TAXID}" ]]; then
  CMD+=(--root-taxid "${ROOT_TAXID}")
fi

if [[ -n "${RANK_KEEP}" ]]; then
  CMD+=(--rank-keep "${RANK_KEEP}")
fi

# -----------------------------
# Run
# -----------------------------

echo "Starting PhyloArk taxonomy quota sampling"
echo "Working directory: $(pwd)"
echo "Genome table: ${GENOMES}"
echo "Taxdump: ${TAXDUMP}"
echo "Target genomes: ${TARGET}"
echo "Alpha: ${ALPHA}"
echo "Root taxid: ${ROOT_TAXID:-inferred LCA}"
echo "Output prefix: ${OUT_PREFIX}"
echo
printf 'Command:'
printf ' %q' "${CMD[@]}"
echo

"${CMD[@]}"

echo
echo "Done. Outputs:"
ls -lh "${OUT_PREFIX}".*