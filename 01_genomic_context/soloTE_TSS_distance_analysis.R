#!/usr/bin/env Rscript

# ==============================================================================
# Distance of expressed SoloTE loci to the nearest gene TSS
#
# This script reproduces the genomic-distance analysis used to test whether
# expressed TE loci shift toward or away from gene promoters across Drosophila
# spermatogenesis.
#
# The script:
#   1. Loads the processed SoloTE Seurat object
#   2. Parses genomic coordinates from SoloTE feature names
#   3. Extracts gene transcription start sites (TSSs) from the w1118 GTF
#   4. Calculates the distance from each TE-locus midpoint to the nearest TSS
#   5. Computes locus-level Z-scores across annotated testis cell types
#   6. Classifies loci as expressed within each cell type when Z >= 1
#   7. Summarizes the fraction of expressed loci across fixed distance bins
#   8. Compares expressed and non-expressed loci using two-sided Wilcoxon
#      rank-sum tests, followed by Benjamini-Hochberg correction
#
# Required input files:
#   FCA-testis-merged-normalized-scaled-soloTE-w1118.rds
#   w1118_v4_sorted.protein_coding.gtf
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
  library(matrixStats)
  library(Seurat)
  library(stringr)
  library(tidyr)
})


# ------------------------------------------------------------------------------
# 2. Input files, output directory, and analysis parameters
# ------------------------------------------------------------------------------

seurat_file <- "FCA-testis-merged-normalized-scaled-soloTE-w1118.rds"
gtf_file <- "w1118_v4_sorted.protein_coding.gtf"

output_dir <- "TSS_distance_results"

z_threshold <- 1
minimum_loci_per_group <- 10
distance_bin_width <- 0.25

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 3. Germline cell types and developmental stages
# ------------------------------------------------------------------------------

mitotic_cell_types <- c(
  "spermatogonium",
  "mid-late proliferating spermatogonia",
  "spermatogonium-spermatocyte transition"
)

meiotic_cell_types <- c(
  "spermatocyte 0",
  "spermatocyte 1",
  "spermatocyte 2",
  "spermatocyte 3",
  "spermatocyte 4",
  "spermatocyte 5",
  "spermatocyte 6",
  "spermatocyte 7a",
  "spermatocyte",
  "late primary spermatocyte"
)

postmeiotic_cell_types <- c(
  "spermatid",
  "early elongation stage spermatid",
  "early-mid elongation-stage spermatid",
  "mid-late elongation-stage spermatid"
)

germline_cell_types <- c(
  mitotic_cell_types,
  meiotic_cell_types,
  postmeiotic_cell_types
)

cell_type_order <- germline_cell_types

stage_lookup <- data.frame(
  cell_type = germline_cell_types,
  stage = c(
    rep("Mitotic", length(mitotic_cell_types)),
    rep("Meiotic", length(meiotic_cell_types)),
    rep("Post-meiotic", length(postmeiotic_cell_types))
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

standardize_chromosome <- function(chromosome) {
  chromosome <- as.character(chromosome)
  chromosome <- sub("^chr", "", chromosome)

  chromosome
}


extract_gtf_attribute <- function(attributes, key) {
  pattern <- paste0(
    "(?:^|;\\s*)",
    key,
    "\\s+\"([^\"]+)\""
  )

  match_matrix <- stringr::str_match(
    attributes,
    pattern
  )

  match_matrix[, 2]
}


parse_soloTE_features <- function(features) {
  normalized_ids <- gsub(
    "_",
    "-",
    features,
    fixed = TRUE
  )

  feature_match <- stringr::str_match(
    normalized_ids,
    paste0(
      "^SoloTE-",
      "(X|2L|2R|3L|3R|4|Y)-",
      "(\\d+)-",
      "(\\d+)-",
      "([^:]+):",
      "([^:]+):",
      "([^:-]+)-",
      "(\\d+(?:\\.\\d+)?)-",
      "([+-])$"
    )
  )

  parsed <- data.frame(
    feature = features,
    chromosome = feature_match[, 2],
    start = suppressWarnings(as.integer(feature_match[, 3])),
    end = suppressWarnings(as.integer(feature_match[, 4])),
    family = feature_match[, 5],
    superfamily = feature_match[, 6],
    class = feature_match[, 7],
    divergence = suppressWarnings(as.numeric(feature_match[, 8])),
    strand = feature_match[, 9],
    stringsAsFactors = FALSE
  )

  parsed <- parsed %>%
    filter(
      !is.na(chromosome),
      !is.na(start),
      !is.na(end)
    ) %>%
    mutate(
      chromosome = standardize_chromosome(chromosome),
      midpoint = floor((start + end) / 2)
    )

  parsed
}


build_gene_tss <- function(gtf_file, chromosomes_to_keep) {
  if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file)
  }

  gtf <- read.delim(
    gtf_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "#",
    check.names = FALSE,
    col.names = c(
      "seqname",
      "source",
      "type",
      "start",
      "end",
      "score",
      "strand",
      "phase",
      "attributes"
    )
  )

  gtf <- gtf %>%
    mutate(
      seqname = standardize_chromosome(seqname),
      gene_id = extract_gtf_attribute(attributes, "gene_id"),
      gene_name = extract_gtf_attribute(attributes, "gene_name")
    ) %>%
    filter(
      seqname %in% chromosomes_to_keep
    )

  genes <- gtf %>%
    filter(
      type == "gene",
      !is.na(gene_id)
    ) %>%
    transmute(
      chromosome = seqname,
      start = as.integer(start),
      end = as.integer(end),
      strand = as.character(strand),
      gene_id = as.character(gene_id),
      gene_name = as.character(gene_name)
    )

  # If explicit gene rows are unavailable, derive gene ranges from all rows
  # carrying a gene_id.
  if (nrow(genes) == 0) {
    genes <- gtf %>%
      filter(
        !is.na(gene_id)
      ) %>%
      group_by(
        seqname,
        strand,
        gene_id
      ) %>%
      summarise(
        start = min(start),
        end = max(end),
        gene_name = dplyr::first(
          gene_name[
            !is.na(gene_name) &
              nzchar(gene_name)
          ],
          default = NA_character_
        ),
        .groups = "drop"
      ) %>%
      transmute(
        chromosome = seqname,
        start = as.integer(start),
        end = as.integer(end),
        strand = as.character(strand),
        gene_id = as.character(gene_id),
        gene_name = as.character(gene_name)
      )
  }

  if (nrow(genes) == 0) {
    stop("No gene coordinates could be extracted from the GTF.")
  }

  genes <- genes %>%
    distinct(
      chromosome,
      start,
      end,
      strand,
      gene_id,
      .keep_all = TRUE
    ) %>%
    mutate(
      gene_name = ifelse(
        is.na(gene_name) | !nzchar(gene_name),
        gene_id,
        gene_name
      ),
      tss = case_when(
        strand == "+" ~ start,
        strand == "-" ~ end,
        TRUE ~ NA_integer_
      )
    ) %>%
    filter(
      !is.na(tss)
    )

  if (nrow(genes) == 0) {
    stop("No stranded gene TSS coordinates could be constructed.")
  }

  genes
}


calculate_nearest_tss <- function(te_metadata, gene_tss) {
  te_midpoints_gr <- GRanges(
    seqnames = te_metadata$chromosome,
    ranges = IRanges(
      start = te_metadata$midpoint,
      width = 1
    ),
    strand = "*"
  )

  tss_gr <- GRanges(
    seqnames = gene_tss$chromosome,
    ranges = IRanges(
      start = gene_tss$tss,
      width = 1
    ),
    strand = gene_tss$strand,
    gene_id = gene_tss$gene_id,
    gene_name = gene_tss$gene_name
  )

  nearest_hits <- distanceToNearest(
    te_midpoints_gr,
    tss_gr,
    ignore.strand = TRUE
  )

  nearest_tss_distance <- rep(
    NA_integer_,
    length(te_midpoints_gr)
  )

  nearest_gene_id <- rep(
    NA_character_,
    length(te_midpoints_gr)
  )

  nearest_gene_name <- rep(
    NA_character_,
    length(te_midpoints_gr)
  )

  nearest_tss_distance[queryHits(nearest_hits)] <-
    mcols(nearest_hits)$distance

  nearest_gene_id[queryHits(nearest_hits)] <-
    as.character(
      mcols(tss_gr)$gene_id[
        subjectHits(nearest_hits)
      ]
    )

  nearest_gene_name[queryHits(nearest_hits)] <-
    as.character(
      mcols(tss_gr)$gene_name[
        subjectHits(nearest_hits)
      ]
    )

  te_metadata %>%
    mutate(
      nearest_gene_id = nearest_gene_id,
      nearest_gene_name = nearest_gene_name,
      nearest_TSS_distance_bp = nearest_tss_distance,
      log10_TSS_distance = log10(
        nearest_TSS_distance_bp + 1
      )
    )
}


assign_stage <- function(cell_type) {
  case_when(
    cell_type %in% mitotic_cell_types ~ "Mitotic",
    cell_type %in% meiotic_cell_types ~ "Meiotic",
    cell_type %in% postmeiotic_cell_types ~ "Post-meiotic",
    TRUE ~ NA_character_
  )
}


run_wilcoxon_test <- function(data_for_cell_type) {
  expressed_distances <- data_for_cell_type %>%
    filter(expressed) %>%
    pull(log10_TSS_distance)

  nonexpressed_distances <- data_for_cell_type %>%
    filter(!expressed) %>%
    pull(log10_TSS_distance)

  expressed_distances <- expressed_distances[
    is.finite(expressed_distances)
  ]

  nonexpressed_distances <- nonexpressed_distances[
    is.finite(nonexpressed_distances)
  ]

  n_expressed <- length(expressed_distances)
  n_nonexpressed <- length(nonexpressed_distances)

  median_expressed <- if (n_expressed > 0) {
    median(expressed_distances)
  } else {
    NA_real_
  }

  median_nonexpressed <- if (n_nonexpressed > 0) {
    median(nonexpressed_distances)
  } else {
    NA_real_
  }

  median_all <- median(
    c(
      expressed_distances,
      nonexpressed_distances
    ),
    na.rm = TRUE
  )

  median_shift <- median_expressed - median_nonexpressed

  if (
    n_expressed < minimum_loci_per_group ||
      n_nonexpressed < minimum_loci_per_group
  ) {
    return(
      data.frame(
        n_expressed = n_expressed,
        n_nonexpressed = n_nonexpressed,
        median_log10_distance_expressed = median_expressed,
        median_log10_distance_nonexpressed = median_nonexpressed,
        median_log10_distance_all = median_all,
        median_shift_expressed_minus_nonexpressed = median_shift,
        wilcoxon_W = NA_real_,
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  test_result <- wilcox.test(
    x = expressed_distances,
    y = nonexpressed_distances,
    alternative = "two.sided",
    exact = FALSE
  )

  data.frame(
    n_expressed = n_expressed,
    n_nonexpressed = n_nonexpressed,
    median_log10_distance_expressed = median_expressed,
    median_log10_distance_nonexpressed = median_nonexpressed,
    median_log10_distance_all = median_all,
    median_shift_expressed_minus_nonexpressed = median_shift,
    wilcoxon_W = unname(test_result$statistic),
    p_value = test_result$p.value,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 5. Load the processed SoloTE Seurat object
# ------------------------------------------------------------------------------

if (!file.exists(seurat_file)) {
  stop("Seurat object not found: ", seurat_file)
}

soloTE_object <- readRDS(
  seurat_file
)

metadata_columns <- colnames(
  soloTE_object[[]]
)

annotation_candidates <- c(
  "annotation",
  "newannotation",
  "S_annotation"
)

annotation_column <- annotation_candidates[
  annotation_candidates %in% metadata_columns
][1]

if (
  length(annotation_column) == 0 ||
    is.na(annotation_column)
) {
  stop(
    "No supported cell-type annotation column was found. ",
    "Expected one of: ",
    paste(annotation_candidates, collapse = ", ")
  )
}

message(
  "Using cell-type annotation column: ",
  annotation_column
)


# ------------------------------------------------------------------------------
# 6. Calculate average expression and locus-level Z-scores
# ------------------------------------------------------------------------------

average_expression <- AverageExpression(
  object = soloTE_object,
  assays = "RNA",
  group.by = annotation_column,
  layer = "data",
  verbose = FALSE
)$RNA

soloTE_features <- rownames(average_expression)[
  grepl(
    "^SoloTE[-_]",
    rownames(average_expression)
  )
]

if (length(soloTE_features) == 0) {
  stop("No SoloTE features were found in the Seurat object.")
}

average_TE_expression <- average_expression[
  soloTE_features,
  ,
  drop = FALSE
]

# Standardize each TE locus across annotated cell types.
TE_expression_z <- t(
  scale(
    t(average_TE_expression)
  )
)

TE_expression_z[
  is.na(TE_expression_z)
] <- 0


# ------------------------------------------------------------------------------
# 7. Parse SoloTE coordinates and calculate nearest-TSS distances
# ------------------------------------------------------------------------------

TE_metadata <- parse_soloTE_features(
  rownames(TE_expression_z)
)

if (nrow(TE_metadata) == 0) {
  stop("No SoloTE feature names could be parsed.")
}

n_unparsed <- nrow(TE_expression_z) - nrow(TE_metadata)

if (n_unparsed > 0) {
  warning(
    n_unparsed,
    " SoloTE features could not be parsed and were excluded."
  )
}

TE_expression_z <- TE_expression_z[
  TE_metadata$feature,
  ,
  drop = FALSE
]

gene_tss <- build_gene_tss(
  gtf_file = gtf_file,
  chromosomes_to_keep = unique(
    TE_metadata$chromosome
  )
)

TE_context <- calculate_nearest_tss(
  te_metadata = TE_metadata,
  gene_tss = gene_tss
)

if (any(!is.finite(TE_context$log10_TSS_distance))) {
  warning(
    sum(!is.finite(TE_context$log10_TSS_distance)),
    " TE loci lack a finite nearest-TSS distance."
  )
}

write.csv(
  TE_context,
  file.path(
    output_dir,
    "soloTE_locus_nearest_TSS_distance.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8. Construct the TE-locus x cell-type analysis table
# ------------------------------------------------------------------------------

available_germline_cell_types <- intersect(
  germline_cell_types,
  colnames(TE_expression_z)
)

missing_germline_cell_types <- setdiff(
  germline_cell_types,
  available_germline_cell_types
)

if (length(missing_germline_cell_types) > 0) {
  warning(
    "The following germline cell types were not present and were excluded:\n",
    paste(missing_germline_cell_types, collapse = "\n")
  )
}

if (length(available_germline_cell_types) == 0) {
  stop("No germline cell types were found in the expression matrix.")
}

TE_expression_long <- as.data.frame(
  TE_expression_z[
    ,
    available_germline_cell_types,
    drop = FALSE
  ]
) %>%
  tibble::rownames_to_column("feature") %>%
  pivot_longer(
    cols = -feature,
    names_to = "cell_type",
    values_to = "z_score"
  ) %>%
  left_join(
    TE_context,
    by = "feature"
  ) %>%
  mutate(
    expressed = z_score >= z_threshold,
    stage = assign_stage(cell_type),
    cell_type = factor(
      cell_type,
      levels = cell_type_order
    ),
    stage = factor(
      stage,
      levels = c(
        "Mitotic",
        "Meiotic",
        "Post-meiotic"
      )
    )
  ) %>%
  arrange(
    stage,
    cell_type,
    feature
  )

write.csv(
  TE_expression_long,
  file.path(
    output_dir,
    "soloTE_TSS_distance_by_celltype.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 9. Fixed-width distance-bin summaries
# ------------------------------------------------------------------------------

finite_distances <- TE_expression_long$log10_TSS_distance[
  is.finite(
    TE_expression_long$log10_TSS_distance
  )
]

if (length(finite_distances) == 0) {
  stop("No finite TSS distances were available for binning.")
}

minimum_break <- floor(
  min(finite_distances)
)

maximum_break <- ceiling(
  max(finite_distances)
) + distance_bin_width

distance_breaks <- seq(
  from = minimum_break,
  to = maximum_break,
  by = distance_bin_width
)

TE_expression_binned <- TE_expression_long %>%
  filter(
    is.finite(log10_TSS_distance)
  ) %>%
  mutate(
    distance_bin = cut(
      log10_TSS_distance,
      breaks = distance_breaks,
      include.lowest = TRUE,
      right = TRUE
    )
  )

distance_bin_lookup <- data.frame(
  distance_bin = levels(
    TE_expression_binned$distance_bin
  ),
  bin_start = head(
    distance_breaks,
    -1
  ),
  bin_end = tail(
    distance_breaks,
    -1
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    bin_center = (
      bin_start + bin_end
    ) / 2
  )

distance_bin_summary <- TE_expression_binned %>%
  group_by(
    stage,
    cell_type,
    distance_bin
  ) %>%
  summarise(
    n_loci = n(),
    n_expressed = sum(expressed),
    fraction_expressed = mean(expressed),
    .groups = "drop"
  ) %>%
  mutate(
    distance_bin = as.character(
      distance_bin
    )
  ) %>%
  left_join(
    distance_bin_lookup,
    by = "distance_bin"
  ) %>%
  arrange(
    stage,
    cell_type,
    bin_center
  )

write.csv(
  distance_bin_summary,
  file.path(
    output_dir,
    "soloTE_TSS_distance_binned_expression_fraction.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 10. Wilcoxon tests: expressed versus non-expressed loci
# ------------------------------------------------------------------------------

wilcoxon_results <- TE_expression_long %>%
  filter(
    is.finite(log10_TSS_distance)
  ) %>%
  group_by(
    stage,
    cell_type
  ) %>%
  group_modify(
    ~ run_wilcoxon_test(.x)
  ) %>%
  ungroup() %>%
  mutate(
    p_adjusted_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    significant_FDR_0.05 = (
      !is.na(p_adjusted_BH) &
        p_adjusted_BH < 0.05
    ),
    direction = case_when(
      median_shift_expressed_minus_nonexpressed < 0 ~
        "Expressed loci are closer to TSS",
      median_shift_expressed_minus_nonexpressed > 0 ~
        "Expressed loci are farther from TSS",
      median_shift_expressed_minus_nonexpressed == 0 ~
        "No median shift",
      TRUE ~
        "Not tested"
    ),
    cell_type = factor(
      cell_type,
      levels = cell_type_order
    ),
    stage = factor(
      stage,
      levels = c(
        "Mitotic",
        "Meiotic",
        "Post-meiotic"
      )
    )
  ) %>%
  arrange(
    stage,
    cell_type
  )

write.csv(
  wilcoxon_results,
  file.path(
    output_dir,
    "soloTE_TSS_distance_Wilcoxon_by_celltype.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 11. Save analysis parameters and session information
# ------------------------------------------------------------------------------

analysis_parameters <- data.frame(
  parameter = c(
    "Seurat object",
    "GTF annotation",
    "Expression threshold",
    "Minimum loci per Wilcoxon group",
    "Distance transformation",
    "Distance-bin width"
  ),
  value = c(
    seurat_file,
    gtf_file,
    paste0("Z >= ", z_threshold),
    minimum_loci_per_group,
    "log10(bp to nearest TSS + 1)",
    distance_bin_width
  ),
  stringsAsFactors = FALSE
)

write.csv(
  analysis_parameters,
  file.path(
    output_dir,
    "analysis_parameters.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "sessionInfo_TSS_distance_analysis.txt"
  )
)

message("\nTSS-distance analysis complete.")
message("Results were written to: ", output_dir)
