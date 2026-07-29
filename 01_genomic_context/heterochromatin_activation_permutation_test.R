#!/usr/bin/env Rscript

# ==============================================================================
# Chromosome-stratified permutation test for enrichment of activated TE loci
# in heterochromatin across Drosophila w1118 testis cell types.
#
# Activated locus: per-locus Z-score > 1.
# Heterochromatin: pericentromeric (PCH) + telomeric (Tel).
# Null: shuffle Eu/PCH/Tel labels among loci within each chromosome while
# keeping each cell type's activated loci fixed.
#
# Required input:
#   FCA-testis-merged-normalized-scaled-soloTE-w1118.rds
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(stringr)
  library(tibble)
  library(tidyr)
})

seurat_file <- "FCA-testis-merged-normalized-scaled-soloTE-w1118.rds"
output_dir <- "heterochromatin_permutation_results"

z_threshold <- 1
n_permutations <- 100000L
random_seed <- 2L
fdr_threshold <- 0.01

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Inclusive compartment boundaries used in the manuscript analysis.
regions <- tribble(
  ~chromosome, ~start,    ~end,      ~compartment,
  "2L",        1L,        34496L,     "Tel",
  "2L",        34497L,    21806874L,  "Eu",
  "2L",        21806875L, 23219612L,  "PCH",
  "2R",        1L,        4651836L,   "PCH",
  "2R",        4651837L,  24535303L,  "Eu",
  "2R",        24535304L, 24550778L,  "Tel",
  "3L",        1L,        29781L,     "Tel",
  "3L",        29782L,    22811581L,  "Eu",
  "3L",        22811582L, 27189634L,  "PCH",
  "3R",        1L,        3415214L,   "PCH",
  "3R",        3415215L,  31253423L,  "Eu",
  "3R",        31253424L, 31281498L,  "Tel",
  "X",         1L,        138429L,    "Tel",
  "X",         138430L,   22110469L,  "Eu",
  "X",         22110470L, 22604182L,  "PCH"
)

mitotic_cell_types <- c(
  "spermatogonium",
  "mid-late proliferating spermatogonia",
  "spermatogonium-spermatocyte transition"
)

meiotic_cell_types <- c(
  "spermatocyte 0", "spermatocyte 1", "spermatocyte 2",
  "spermatocyte 3", "spermatocyte 4", "spermatocyte 5",
  "spermatocyte 6", "spermatocyte 7a", "spermatocyte",
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

somatic_cell_types <- c(
  "secretory cell of the male reproductive tract",
  "pigment cell",
  "male gonad associated epithelium",
  "muscle cell",
  "hemocyte",
  "adult fat body",
  "adult tracheal cell",
  "head cyst cell",
  "spermatocyte cyst cell branch b",
  "spermatocyte cyst cell branch a",
  "late cyst cell branch b",
  "late cyst cell branch a",
  "cyst cell branch b",
  "cyst cell branch a",
  "cyst cell intermediate",
  "early cyst cell 2",
  "early cyst cell 1",
  "cyst stem cell"
)

cell_type_order <- c(somatic_cell_types, germline_cell_types)

assign_category <- function(cell_type) {
  case_when(
    cell_type %in% germline_cell_types ~ "Germline",
    cell_type %in% somatic_cell_types ~ "Somatic",
    TRUE ~ NA_character_
  )
}

assign_stage <- function(cell_type) {
  case_when(
    cell_type %in% somatic_cell_types ~ "Somatic",
    cell_type %in% mitotic_cell_types ~ "Mitotic",
    cell_type %in% meiotic_cell_types ~ "Meiotic",
    cell_type %in% postmeiotic_cell_types ~ "Post-meiotic",
    TRUE ~ NA_character_
  )
}

average_by_cell_type <- function(object, annotation_column) {
  tryCatch(
    AverageExpression(
      object,
      assays = "RNA",
      group.by = annotation_column,
      layer = "data",
      verbose = FALSE
    )$RNA,
    error = function(e) {
      AverageExpression(
        object,
        assays = "RNA",
        group.by = annotation_column,
        slot = "data",
        verbose = FALSE
      )$RNA
    }
  )
}

parse_soloTE <- function(features) {
  normalized <- gsub("_", "-", features, fixed = TRUE)

  match_matrix <- str_match(
    normalized,
    paste0(
      "^SoloTE-(X|2L|2R|3L|3R|4|Y)-",
      "(\\d+)-(\\d+)-([^:]+):([^:]+):([^:-]+)-",
      "(\\d+(?:\\.\\d+)?)-([+-])$"
    )
  )

  tibble(
    feature = features,
    chromosome = match_matrix[, 2],
    start = suppressWarnings(as.integer(match_matrix[, 3])),
    end = suppressWarnings(as.integer(match_matrix[, 4])),
    TE_family = match_matrix[, 5],
    TE_superfamily = match_matrix[, 6],
    TE_class = toupper(match_matrix[, 7]),
    divergence = suppressWarnings(as.numeric(match_matrix[, 8])),
    strand = match_matrix[, 9]
  ) %>%
    filter(!is.na(chromosome), !is.na(start), !is.na(end)) %>%
    mutate(midpoint = floor((start + end) / 2))
}

get_compartment <- function(chromosome, midpoint) {
  matched <- regions %>%
    filter(
      .data$chromosome == .env$chromosome,
      .env$midpoint >= start,
      .env$midpoint <= end
    )

  if (nrow(matched) == 1) {
    return(matched$compartment)
  }

  if (nrow(matched) > 1) {
    stop("Locus matched multiple compartments: ", chromosome, ":", midpoint)
  }

  NA_character_
}

# ------------------------------------------------------------------------------
# Load object and calculate per-locus Z-scores
# ------------------------------------------------------------------------------

if (!file.exists(seurat_file)) {
  stop("Input file not found: ", seurat_file)
}

soloTE_object <- readRDS(seurat_file)

annotation_candidates <- c("annotation", "newannotation", "S_annotation")
annotation_column <- annotation_candidates[
  annotation_candidates %in% colnames(soloTE_object[[]])
][1]

if (length(annotation_column) == 0 || is.na(annotation_column)) {
  stop(
    "No annotation column found. Expected one of: ",
    paste(annotation_candidates, collapse = ", ")
  )
}

annotations <- as.character(soloTE_object[[]][[annotation_column]])
names(annotations) <- Cells(soloTE_object)

valid_cells <- (
  !is.na(annotations) &
  nzchar(annotations) &
  annotations != "unannotated" &
  toupper(annotations) != "NA"
)

soloTE_object <- subset(
  soloTE_object,
  cells = names(annotations)[valid_cells]
)

present_annotations <- unique(
  as.character(soloTE_object[[]][[annotation_column]])
)

unknown_cell_types <- setdiff(present_annotations, cell_type_order)

if (length(unknown_cell_types) > 0) {
  stop(
    "Unassigned cell types found:\n",
    paste(sort(unknown_cell_types), collapse = "\n")
  )
}

average_expression <- average_by_cell_type(
  soloTE_object,
  annotation_column
)

TE_features <- rownames(average_expression)[
  grepl("^SoloTE[-_]", rownames(average_expression))
]

if (length(TE_features) == 0) {
  stop("No SoloTE locus features were found.")
}

TE_average <- average_expression[TE_features, , drop = FALSE]
TE_z <- t(scale(t(TE_average)))
TE_z[is.na(TE_z)] <- 0

cell_types <- intersect(cell_type_order, colnames(TE_z))

if (length(cell_types) == 0) {
  stop("No predefined cell types were found in the expression matrix.")
}

TE_z <- TE_z[, cell_types, drop = FALSE]

# ------------------------------------------------------------------------------
# Parse coordinates and assign Eu/PCH/Tel
# ------------------------------------------------------------------------------

TE_metadata <- parse_soloTE(rownames(TE_z))

if (nrow(TE_metadata) == 0) {
  stop("No SoloTE coordinates could be parsed.")
}

TE_metadata$compartment <- mapply(
  get_compartment,
  TE_metadata$chromosome,
  TE_metadata$midpoint,
  USE.NAMES = FALSE
)

write.csv(
  TE_metadata %>% filter(is.na(compartment)),
  file.path(output_dir, "TE_loci_excluded_from_compartment_test.csv"),
  row.names = FALSE
)

TE_metadata <- TE_metadata %>%
  filter(compartment %in% c("Eu", "PCH", "Tel")) %>%
  distinct(feature, .keep_all = TRUE)

TE_z <- TE_z[TE_metadata$feature, , drop = FALSE]

if (!identical(rownames(TE_z), TE_metadata$feature)) {
  stop("TE metadata and expression rows are not aligned.")
}

write.csv(
  TE_metadata,
  file.path(output_dir, "TE_locus_compartment_map.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Activated-locus counts
# ------------------------------------------------------------------------------

activation_matrix <- TE_z > z_threshold
storage.mode(activation_matrix) <- "double"

active_counts <- colSums(activation_matrix)

activation_long <- as.data.frame(TE_z) %>%
  rownames_to_column("feature") %>%
  pivot_longer(
    cols = -feature,
    names_to = "cell_type",
    values_to = "z_score"
  ) %>%
  left_join(TE_metadata, by = "feature") %>%
  mutate(
    activated = z_score > z_threshold,
    broad_category = assign_category(cell_type),
    developmental_stage = assign_stage(cell_type),
    cell_type = factor(cell_type, levels = cell_type_order)
  )

write.csv(
  activation_long %>% filter(activated),
  file.path(output_dir, "activated_TE_loci_by_celltype.csv"),
  row.names = FALSE
)

activated_counts <- activation_long %>%
  filter(activated) %>%
  count(
    broad_category,
    developmental_stage,
    cell_type,
    compartment,
    name = "n_activated_loci"
  ) %>%
  group_by(broad_category, developmental_stage, cell_type) %>%
  mutate(
    n_activated_total = sum(n_activated_loci),
    fraction_of_activated_loci = n_activated_loci / n_activated_total
  ) %>%
  ungroup() %>%
  arrange(cell_type, factor(compartment, levels = c("Eu", "PCH", "Tel")))

write.csv(
  activated_counts,
  file.path(output_dir, "activated_TE_compartment_counts_by_celltype.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Observed heterochromatin fraction
# ------------------------------------------------------------------------------

is_heterochromatic <- TE_metadata$compartment %in% c("PCH", "Tel")

observed_het_counts <- as.numeric(
  crossprod(as.numeric(is_heterochromatic), activation_matrix)
)

observed_het_fraction <- observed_het_counts / active_counts
names(observed_het_fraction) <- colnames(activation_matrix)

# ------------------------------------------------------------------------------
# Chromosome-stratified permutations
# ------------------------------------------------------------------------------

chromosome_groups <- split(
  seq_len(nrow(TE_metadata)),
  TE_metadata$chromosome
)

write.csv(
  TE_metadata %>%
    count(chromosome, compartment, name = "n_loci") %>%
    group_by(chromosome) %>%
    mutate(
      chromosome_total = sum(n_loci),
      chromosome_fraction = n_loci / chromosome_total
    ) %>%
    ungroup(),
  file.path(output_dir, "TE_compartment_counts_by_chromosome.csv"),
  row.names = FALSE
)

set.seed(random_seed)

null_proportions <- matrix(
  NA_real_,
  nrow = n_permutations,
  ncol = length(cell_types),
  dimnames = list(NULL, cell_types)
)

message(
  "Running ",
  format(n_permutations, big.mark = ","),
  " chromosome-stratified permutations..."
)

for (b in seq_len(n_permutations)) {
  permuted_het <- is_heterochromatic

  for (indices in chromosome_groups) {
    permuted_het[indices] <- sample(
      is_heterochromatic[indices],
      length(indices),
      replace = FALSE
    )
  }

  permuted_counts <- as.numeric(
    crossprod(as.numeric(permuted_het), activation_matrix)
  )

  valid <- active_counts > 0
  null_proportions[b, valid] <- permuted_counts[valid] / active_counts[valid]

  if (b %% 10000L == 0L) {
    message("Completed ", format(b, big.mark = ","), " permutations.")
  }
}

# ------------------------------------------------------------------------------
# Two-sided permutation P values and BH correction
# ------------------------------------------------------------------------------

results <- lapply(seq_along(cell_types), function(i) {
  cell_type <- cell_types[i]
  null_values <- null_proportions[, i]
  null_values <- null_values[is.finite(null_values)]
  observed <- observed_het_fraction[cell_type]

  if (!is.finite(observed) || length(null_values) == 0) {
    return(
      tibble(
        cell_type = cell_type,
        n_activated_loci = active_counts[cell_type],
        observed_heterochromatic_fraction = observed,
        null_mean = NA_real_,
        null_median = NA_real_,
        null_sd = NA_real_,
        null_lower_95 = NA_real_,
        null_upper_95 = NA_real_,
        difference_observed_minus_null = NA_real_,
        z_enrichment = NA_real_,
        p_value_two_sided = NA_real_
      )
    )
  }

  null_mean <- mean(null_values)
  null_sd <- sd(null_values)
  null_ci <- quantile(null_values, c(0.025, 0.975), names = FALSE)

  observed_deviation <- abs(observed - null_mean)
  null_deviation <- abs(null_values - null_mean)

  p_two_sided <- (
    1 + sum(null_deviation >= observed_deviation)
  ) / (
    length(null_values) + 1
  )

  tibble(
    cell_type = cell_type,
    n_activated_loci = active_counts[cell_type],
    observed_heterochromatic_fraction = observed,
    null_mean = null_mean,
    null_median = median(null_values),
    null_sd = null_sd,
    null_lower_95 = null_ci[1],
    null_upper_95 = null_ci[2],
    difference_observed_minus_null = observed - null_mean,
    z_enrichment = ifelse(
      is.finite(null_sd) && null_sd > 0,
      (observed - null_mean) / null_sd,
      NA_real_
    ),
    p_value_two_sided = p_two_sided
  )
}) %>%
  bind_rows() %>%
  mutate(
    broad_category = assign_category(cell_type),
    developmental_stage = assign_stage(cell_type),
    p_adjusted_BH = p.adjust(p_value_two_sided, method = "BH"),
    significant_FDR_0.01 = (
      !is.na(p_adjusted_BH) &
      p_adjusted_BH < fdr_threshold
    ),
    direction = case_when(
      significant_FDR_0.01 & z_enrichment > 0 ~
        "Heterochromatin over-enriched",
      significant_FDR_0.01 & z_enrichment < 0 ~
        "Heterochromatin depleted",
      TRUE ~
        "Not significant"
    ),
    cell_type = factor(cell_type, levels = cell_type_order)
  ) %>%
  arrange(cell_type)

write.csv(
  results,
  file.path(
    output_dir,
    "heterochromatin_enrichment_permutation_results.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    activation_z_threshold = z_threshold,
    n_permutations = n_permutations,
    random_seed = random_seed,
    cell_types = cell_types,
    observed_heterochromatic_fraction = observed_het_fraction,
    null_heterochromatic_fractions = null_proportions
  ),
  file.path(
    output_dir,
    "heterochromatin_permutation_null_distributions.rds"
  ),
  compress = "xz"
)

write.csv(
  tibble(
    cell_type = cell_types,
    null_mean = colMeans(null_proportions, na.rm = TRUE),
    null_median = apply(null_proportions, 2, median, na.rm = TRUE),
    null_sd = apply(null_proportions, 2, sd, na.rm = TRUE),
    null_lower_95 = apply(
      null_proportions,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    null_upper_95 = apply(
      null_proportions,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    )
  ),
  file.path(
    output_dir,
    "heterochromatin_permutation_null_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    parameter = c(
      "Seurat object",
      "Activation definition",
      "Heterochromatin definition",
      "Permutation labels",
      "Permutation strata",
      "Number of permutations",
      "Random seed",
      "Test",
      "Multiple-testing correction",
      "FDR threshold",
      "Excluded chromosomes"
    ),
    value = c(
      seurat_file,
      paste0("Per-locus Z-score > ", z_threshold),
      "PCH + Tel",
      "Eu/PCH/Tel",
      "Chromosome",
      n_permutations,
      random_seed,
      "Two-sided relative to the null mean",
      "Benjamini-Hochberg across cell types",
      fdr_threshold,
      "4 and Y"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "analysis_parameters.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "sessionInfo_heterochromatin_permutation_test.txt"
  )
)

message("\nPermutation analysis complete.")
message("Results were written to: ", output_dir)
