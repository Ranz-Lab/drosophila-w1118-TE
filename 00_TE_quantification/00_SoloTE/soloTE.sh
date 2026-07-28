```bash
#!/usr/bin/env bash

# ==============================================================================
# Run SoloTE on Drosophila melanogaster w1118 testis sample 1
#
# Input:
#   Cell barcode- and UMI-tagged BAM file
#
# Output prefix:
#   FCA59_Sample1_w1118
# ==============================================================================

#SBATCH --account=JRANZ_LAB
#SBATCH --job-name=Male_testis_sample1_SoloTE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=90G
#SBATCH --time=48:00:00
#SBATCH --error=error/%x.%A.err
#SBATCH --output=error/%x.%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rongyinl@uci.edu

set -euo pipefail

# Add required software to PATH.
export PATH="/dfs5/bio/rongyinl/TE/samtools-1.19:${PATH}"
export PATH="/dfs5/bio/rongyinl/TE/Polish:${PATH}"
export PATH="/opt/apps/R/4.2.2/bin:${PATH}"

SOLOTE_PIPELINE="/dfs5/bio/rongyinl/TE/Download/SoloTE/SoloTE_pipeline.py"
INPUT_BAM="/dfs5/bio/rongyinl/TE/run1-FCA59_Male_testis_adult_1dWT_Fuller_sample1_S44_w1118-clean.bam"
TE_ANNOTATION="/share/crsp/lab/jranz/share/SingleCell_TE/w1118/w1118_TE-annotation_onecodetofindthemall.bed"
OUTPUT_DIR="/dfs5/bio/rongyinl/TE/Download"
OUTPUT_PREFIX="FCA59_Sample1_w1118"

echo "Starting SoloTE analysis"
echo "Input BAM: ${INPUT_BAM}"
echo "TE annotation: ${TE_ANNOTATION}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "Threads: ${SLURM_CPUS_PER_TASK}"
echo "Start time: $(date)"

python "${SOLOTE_PIPELINE}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --bam "${INPUT_BAM}" \
    --teannotation "${TE_ANNOTATION}" \
    --outputprefix "${OUTPUT_PREFIX}" \
    --outputdir "${OUTPUT_DIR}"

echo "SoloTE analysis completed successfully."
echo "End time: $(date)"
```
