#!/bin/bash
#SBATCH --job-name=H3K9me2_CUTTag
#SBATCH -A JRANZ_LAB
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=72:00:00
#SBATCH --error=%x.%A.err
#SBATCH --output=%x.%A.out

# ==============================================================================
# H3K9me2 CUT&Tag processing workflow
#
#   1. Download SRA accessions and create paired FASTQ files
#   2. Align paired-end reads to the w1118 genome with Bowtie2
#   3. Retain mapped reads with MAPQ >= 2
#   4. Sort and index replicate BAM files
#   5. Merge the two spermatocyte and two wing replicates
#   6. Call broad H3K9me2 peaks with MACS2
#   7. Generate 50-bp RPGC-normalized bigWigs
#   8. Generate 20-kb RPGC-normalized bigWigs
#   9. Generate 20-kb non-normalized bigWigs
#  10. Generate the 20-kb multiBigwigSummary table
#
# The original analytical parameters and output filenames are retained.
#
# Required file:
#   sra-list.txt
#
# Required pre-existing Bowtie2 index prefix:
#   ragtag_bowtie2_index
#
# Run from:
#   /dfs10/bio/ihariyan/ranz_lab/CUT_TAG
# ==============================================================================

set -euo pipefail
shopt -s nullglob


# ------------------------------------------------------------------------------
# 1. Fixed parameters
# ------------------------------------------------------------------------------

CORES=8
MIN_QUALITY_SCORE=2

CUTTAG_DIR="/dfs10/bio/ihariyan/ranz_lab/CUT_TAG"

REFERENCE_GENOME="/pub/ihariyan/ranz_lab/data/genome/drosophila_melanogaster/w1118/scaffolded/ragtag.scaffold.chr.fasta"
BOWTIE2_INDEX="ragtag_bowtie2_index"

EFFECTIVE_GENOME_SIZE=142573017

DEEPTOOLS_ENV="/data/homezvol0/ihariyan/.conda/envs/deeptools_coverage"

SRA_LIST="sra-list.txt"

H3K9_SPERM_REPLICATES=(
    H3K9me2_spermatocyte_rep1_bowtie2.q2.sorted.bam
    H3K9me2_spermatocyte_rep2_bowtie2.q2.sorted.bam
)

H3K9_WING_REPLICATES=(
    H3K9me2_wing_rep1_bowtie2.q2.sorted.bam
    H3K9me2_wing_rep2_bowtie2.q2.sorted.bam
)


# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

require_file() {
    local file="$1"

    if [[ ! -f "$file" ]]; then
        echo "ERROR: Required file not found: $file" >&2
        exit 1
    fi
}


require_files() {
    local file

    for file in "$@"; do
        require_file "$file"
    done
}


require_directory() {
    local directory="$1"

    if [[ ! -d "$directory" ]]; then
        echo "ERROR: Required directory not found: $directory" >&2
        exit 1
    fi
}


log_step() {
    echo
    echo "=============================================================================="
    echo "$1"
    echo "=============================================================================="
}


# ------------------------------------------------------------------------------
# 3. Validate the working directory and fixed inputs
# ------------------------------------------------------------------------------

require_directory "$CUTTAG_DIR"
cd "$CUTTAG_DIR"

require_file "$SRA_LIST"

echo "Working directory: $PWD"
echo "Reference genome:  $REFERENCE_GENOME"
echo "Bowtie2 index:     $BOWTIE2_INDEX"
echo "CPUs:              $CORES"
echo "MAPQ threshold:    $MIN_QUALITY_SCORE"


# ------------------------------------------------------------------------------
# 4. Load software modules
# ------------------------------------------------------------------------------

module load sra-tools/3.0.0
module load bowtie2/2.5.1
module load samtools
module load macs


# ------------------------------------------------------------------------------
# 5. Download SRA accessions and create paired FASTQ files
# ------------------------------------------------------------------------------

log_step "STEP 1: Downloading SRA accessions"

prefetch --option-file "$SRA_LIST"


log_step "STEP 2: Creating paired FASTQ files"

sra_directories=(SRR*/)

if [[ ${#sra_directories[@]} -eq 0 ]]; then
    echo "ERROR: No SRR directories were found after prefetch." >&2
    exit 1
fi

for sra_directory in "${sra_directories[@]}"; do
    accession="${sra_directory%/}"
    sra_file="${accession}/${accession}.sra"

    require_file "$sra_file"

    echo "FASTQ conversion: $accession"

    fastq-dump \
        --split-files \
        -O "$accession" \
        "$sra_file"
done


# ------------------------------------------------------------------------------
# 6. Align paired-end reads with Bowtie2
#
# The index-building command from the original script is retained below as a
# comment because the workflow used an already-built index.
# ------------------------------------------------------------------------------

log_step "STEP 3: Aligning paired-end reads with Bowtie2"

# bowtie2-build "$REFERENCE_GENOME" "$BOWTIE2_INDEX"

for sra_directory in "${sra_directories[@]}"; do
    accession="${sra_directory%/}"

    fastq_1="${CUTTAG_DIR}/${accession}/${accession}_1.fastq"
    fastq_2="${CUTTAG_DIR}/${accession}/${accession}_2.fastq"

    sam_output="${accession}_bowtie2.sam"
    log_output="${accession}_bowtie2.txt"

    require_files \
        "$fastq_1" \
        "$fastq_2"

    echo "Alignment: $accession"

    bowtie2 \
        --end-to-end \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        --phred33 \
        -I 10 \
        -X 700 \
        -p "$CORES" \
        -x "$BOWTIE2_INDEX" \
        -1 "$fastq_1" \
        -2 "$fastq_2" \
        -S "$sam_output" \
        &> "$log_output"
done


# ------------------------------------------------------------------------------
# 7. Filter mapped reads at MAPQ >= 2
# ------------------------------------------------------------------------------

log_step "STEP 4: Filtering mapped reads at MAPQ >= 2"

alignment_sams=(SRR*_bowtie2.sam)

if [[ ${#alignment_sams[@]} -eq 0 ]]; then
    echo "ERROR: No Bowtie2 SAM files were found." >&2
    exit 1
fi

for sam_file in "${alignment_sams[@]}"; do
    base="${sam_file%_bowtie2.sam}"
    mapped_bam="${base}_bowtie2.q${MIN_QUALITY_SCORE}.mapped.bam"

    echo "Filtering: $sam_file -> $mapped_bam"

    samtools view \
        -bS \
        -q "$MIN_QUALITY_SCORE" \
        -F 0x04 \
        "$sam_file" \
        > "$mapped_bam"
done


# ------------------------------------------------------------------------------
# 8. Sort and index per-replicate BAM files
# ------------------------------------------------------------------------------

log_step "STEP 5: Sorting and indexing per-replicate BAM files"

mapped_bams=(*_bowtie2.q*.mapped.bam)

if [[ ${#mapped_bams[@]} -eq 0 ]]; then
    echo "ERROR: No mapped BAM files were found." >&2
    exit 1
fi

for mapped_bam in "${mapped_bams[@]}"; do
    base="${mapped_bam%.mapped.bam}"
    sorted_bam="${base}.sorted.bam"

    if [[ -f "$sorted_bam" ]]; then
        echo "Skipping sort; output already exists: $sorted_bam"
    else
        echo "Sorting: $mapped_bam -> $sorted_bam"

        samtools sort \
            -@ "$CORES" \
            -o "$sorted_bam" \
            "$mapped_bam"
    fi

    if [[ -f "${sorted_bam}.bai" ]]; then
        echo "Skipping index; output already exists: ${sorted_bam}.bai"
    else
        samtools index "$sorted_bam"
    fi
done


# ------------------------------------------------------------------------------
# 9. Validate the H3K9me2 replicate BAM filenames
# ------------------------------------------------------------------------------

require_files \
    "${H3K9_SPERM_REPLICATES[@]}" \
    "${H3K9_WING_REPLICATES[@]}"


# ------------------------------------------------------------------------------
# 10. Merge H3K9me2 replicates
# ------------------------------------------------------------------------------

log_step "STEP 6: Merging H3K9me2 spermatocyte replicates"

printf '  %s\n' "${H3K9_SPERM_REPLICATES[@]}"

samtools merge \
    -@ "$CORES" \
    H3K9me2_sperm.merged.bam \
    "${H3K9_SPERM_REPLICATES[@]}"

samtools sort \
    -@ "$CORES" \
    -o H3K9me2_sperm.merged.sorted.bam \
    H3K9me2_sperm.merged.bam

samtools index \
    H3K9me2_sperm.merged.sorted.bam


log_step "STEP 7: Merging H3K9me2 wing replicates"

printf '  %s\n' "${H3K9_WING_REPLICATES[@]}"

samtools merge \
    -@ "$CORES" \
    H3K9me2_wing.merged.bam \
    "${H3K9_WING_REPLICATES[@]}"

samtools sort \
    -@ "$CORES" \
    -o H3K9me2_wing.merged.sorted.bam \
    H3K9me2_wing.merged.bam

samtools index \
    H3K9me2_wing.merged.sorted.bam


# ------------------------------------------------------------------------------
# 11. Call broad H3K9me2 peaks with MACS2
# ------------------------------------------------------------------------------

log_step "STEP 8: Calling broad H3K9me2 peaks"

macs2 callpeak \
    -t H3K9me2_sperm.merged.sorted.bam \
    -f BAMPE \
    -g dm \
    -n H3K9me2_sperm \
    --keep-dup all \
    --nomodel \
    --broad \
    --broad-cutoff 0.05

macs2 callpeak \
    -t H3K9me2_wing.merged.sorted.bam \
    -f BAMPE \
    -g dm \
    -n H3K9me2_wing \
    --keep-dup all \
    --nomodel \
    --broad \
    --broad-cutoff 0.05


# ------------------------------------------------------------------------------
# 12. Activate the original deepTools environment
# ------------------------------------------------------------------------------

log_step "STEP 9: Activating the deepTools environment"

module load mamba
eval "$(mamba shell hook --shell bash)"
mamba activate "$DEEPTOOLS_ENV"


# ------------------------------------------------------------------------------
# 13. Generate 50-bp RPGC-normalized bigWigs
# ------------------------------------------------------------------------------

log_step "STEP 10: Generating 50-bp RPGC-normalized H3K9me2 bigWigs"

bamCoverage \
    -b H3K9me2_sperm.merged.sorted.bam \
    -o H3K9me2_sperm.RPGC.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
    --binSize 50 \
    --extendReads \
    -p "$CORES"

bamCoverage \
    -b H3K9me2_wing.merged.sorted.bam \
    -o H3K9me2_wing.RPGC.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
    --binSize 50 \
    --extendReads \
    -p "$CORES"


# ------------------------------------------------------------------------------
# 14. Generate 20-kb RPGC-normalized bigWigs
# ------------------------------------------------------------------------------

log_step "STEP 11: Generating 20-kb RPGC-normalized H3K9me2 bigWigs"

bamCoverage \
    -b H3K9me2_sperm.merged.sorted.bam \
    -o H3K9me2_sperm.20kb.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
    --binSize 20000 \
    --centerReads \
    -p "$CORES"

bamCoverage \
    -b H3K9me2_wing.merged.sorted.bam \
    -o H3K9me2_wing.20kb.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
    --binSize 20000 \
    --centerReads \
    -p "$CORES"


# ------------------------------------------------------------------------------
# 15. Generate 20-kb non-normalized bigWigs
# ------------------------------------------------------------------------------

log_step "STEP 12: Generating 20-kb non-normalized H3K9me2 bigWigs"

bamCoverage \
    -b H3K9me2_sperm.merged.sorted.bam \
    -o H3K9me2_sperm.20kb.nonorm.bw \
    --binSize 20000 \
    --centerReads \
    -p "$CORES"

bamCoverage \
    -b H3K9me2_wing.merged.sorted.bam \
    -o H3K9me2_wing.20kb.nonorm.bw \
    --binSize 20000 \
    --centerReads \
    -p "$CORES"


# ------------------------------------------------------------------------------
# 16. Generate the 20-kb multiBigwigSummary table
# ------------------------------------------------------------------------------

log_step "STEP 13: Generating the 20-kb H3K9me2 summary table"

multiBigwigSummary bins \
    -b \
        H3K9me2_sperm.20kb.nonorm.bw \
        H3K9me2_wing.20kb.nonorm.bw \
    -bs 20000 \
    --outRawCounts H3K9me2_sperm_wing_20kb_nonorm_counts.tab \
    -o H3K9me2_sperm_wing_20kb_nonorm.npz \
    -p "$CORES"


# ------------------------------------------------------------------------------
# 17. Final output summary
# ------------------------------------------------------------------------------

log_step "H3K9me2 CUT&Tag workflow complete"

cat <<'EOF'
Primary downstream files:

Merged BAMs:
  H3K9me2_sperm.merged.sorted.bam
  H3K9me2_wing.merged.sorted.bam

RPGC-normalized bigWigs:
  H3K9me2_sperm.RPGC.bw
  H3K9me2_wing.RPGC.bw

20-kb RPGC-normalized bigWigs:
  H3K9me2_sperm.20kb.bw
  H3K9me2_wing.20kb.bw

20-kb non-normalized bigWigs:
  H3K9me2_sperm.20kb.nonorm.bw
  H3K9me2_wing.20kb.nonorm.bw

20-kb summary:
  H3K9me2_sperm_wing_20kb_nonorm_counts.tab
  H3K9me2_sperm_wing_20kb_nonorm.npz
EOF
