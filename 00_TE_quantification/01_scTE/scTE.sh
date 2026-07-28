#!/usr/bin/env bash

# ==============================================================================
# Run scTE on Drosophila melanogaster w1118 testis sample 1
#
# Input:
#   Cell barcode- and UMI-tagged BAM file
#
# Output:
#   Male_testis_sample1_S44Out.h5ad
# ==============================================================================

#SBATCH --account=JRANZ_LAB
#SBATCH --job-name=Male_testis_sample1_S44
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=18G
#SBATCH --time=05:00:00
#SBATCH --error=error/%x.%A.err
#SBATCH --output=error/%x.%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rongyinl@uci.edu

set -euo pipefail

module load samtools/1.15.1

INPUT_BAM="/dfs5/bio/rongyinl/TE/run1-FCA59_Male_testis_adult_1dWT_Fuller_sample1_S44_w1118-clean.bam"
SCTE_INDEX="/share/crsp/lab/jranz/share/SingleCell_TE/w1118/v1-w1118.exclusive.idx"
OUTPUT_PREFIX="/dfs5/bio/rongyinl/TE/scTE_Files/Male_testis_sample1_S44Out"

echo "Starting scTE analysis"
echo "Input BAM: ${INPUT_BAM}"
echo "scTE index: ${SCTE_INDEX}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "Start time: $(date)"

scTE \
    -i "${INPUT_BAM}" \
    -o "${OUTPUT_PREFIX}" \
    -x "${SCTE_INDEX}" \
    --hdf5 True \
    -CB CB \
    -UMI UB

echo "scTE analysis completed successfully."
echo "End time: $(date)"
