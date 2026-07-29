#!/usr/bin/env Rscript

# ==============================================================================
# Protein-coding gene expression across euchromatin and heterochromatin
#
# This script provides the protein-coding gene control for the genomic-context
# analysis. It:
#
#   1. Assigns protein-coding genes to euchromatin (Eu) or pericentromeric
#      heterochromatin (PCH) using gene midpoints
#   2. Calculates replicate-by-cell-type pseudobulk expression
#   3. Z-scores each gene across cell types within each biological replicate
#   4. Aggregates gene-specific Z-scores by replicate, lineage, and compartment
#   5. Fits the mixed-effects model:
#
#        mean_z ~ region * lineage
#        random intercept: replicate
#
#   6. Tests Eu-versus-PCH differences within germline and somatic lineages
#
# Telomeric heterochromatin is excluded from the model because no matched
# protein-coding genes were available in that compartment.
#
# Required input files:
#   FCA-testis-merged-umap-final-soloTE-class-w1118.rds
#   w1118_v4_sorted.protein_coding.gtf
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(GenomicRanges)
  library(IRanges)
  library(Matrix)
  library(nlme)
  library(Seurat)
  library(stringr)
  library(tibble)
  library(tidyr)
})

options(contrasts = c("contr.sum", "contr.poly"))


# ------------------------------------------------------------------------------
# 2. Inputs, outputs, and parameters
# ------------------------------------------------------------------------------

seurat_file <- "FCA-testis-merged-umap-final-soloTE-class-w1118.rds"
gtf_file <- "w1118_v4_sorted.protein_coding.gtf"

output_dir <- "protein_coding_gene_compartment_results"

minimum_expressed_groups <- 3L

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 3. Cell-type definitions
# ------------------------------------------------------------------------------

germline_cell_types <- c(
  "spermatogonium",
  "mid-late proliferating spermatogonia",
  "spermatogonium-spermatocyte transition",
  "spermatocyte 0",
  "spermatocyte 1",
  "spermatocyte 2",
  "spermatocyte 3",
  "spermatocyte 4",
  "spermatocyte 5",
  "spermatocyte 6",
  "spermatocyte 7a",
  "spermatocyte",
  "late primary spermatocyte",
  "spermatid",
  "early elongation stage spermatid",
  "early-mid elongation-stage spermatid",
  "mid-late elongation-stage spermatid"
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

known_cell_types <- c(
  somatic_cell_types,
  germline_cell_types
)


# ------------------------------------------------------------------------------
# 4. Genomic-compartment boundaries
# ------------------------------------------------------------------------------

region_boundaries <- tribble(
  ~chromosome, ~start,    ~end,      ~region,
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


# ------------------------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------------------------

standardize_chromosome <- function(chromosome) {
  sub("^chr", "", as.character(chromosome))
}


extract_gtf_attribute <- function(attributes, key) {
  pattern <- paste0(
    "(?:^|;\\s*)",
    key,
    "\\s+\"([^\"]+)\""
  )

  stringr::str_match(
    attributes,
    pattern
  )[, 2]
}


get_normalized_matrix <- function(seurat_object) {
  if (inherits(seurat_object[["RNA"]], "Assay5")) {
    return(
      GetAssayData(
        seurat_object,
        assay = "RNA",
        layer = "data"
      )
    )
  }

  GetAssayData(
    seurat_object,
    assay = "RNA",
    slot = "data"
  )
}


assign_lineage <- function(cell_type) {
  case_when(
    cell_type %in% germline_cell_types ~ "Germline",
    cell_type %in% somatic_cell_types ~ "Somatic",
    TRUE ~ NA_character_
  )
}


extract_replicate <- function(cell_names) {
  replicate_number <- stringr::str_match(
    cell_names,
    "_([123])$"
  )[, 2]

  if (any(is.na(replicate_number))) {
    stop(
      "Some cell names do not end in _1, _2, or _3. ",
      "Replicate assignment cannot be completed."
    )
  }

  paste0("rep", replicate_number)
}


read_gene_annotation <- function(gtf_file) {
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
      "chromosome",
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

  genes <- gtf %>%
    mutate(
      chromosome = standardize_chromosome(chromosome),
      gene_id_raw = extract_gtf_attribute(attributes, "gene_id"),
      gene_name = extract_gtf_attribute(attributes, "gene_name")
    ) %>%
    filter(
      type == "gene",
      chromosome %in% unique(region_boundaries$chromosome),
      !is.na(gene_id_raw)
    ) %>%
    transmute(
      chromosome,
      start = as.integer(start),
      end = as.integer(end),
      strand = as.character(strand),
      gene_id_raw = as.character(gene_id_raw),
      gene_name = as.character(gene_name)
    ) %>%
    mutate(
      gene_id = ifelse(
        !is.na(gene_name) & nzchar(gene_name),
        gene_name,
        gene_id_raw
      )
    ) %>%
    distinct(
      gene_id_raw,
      chromosome,
      start,
      end,
      strand,
      .keep_all = TRUE
    )

  if (nrow(genes) == 0) {
    stop("No protein-coding gene records were extracted from the GTF.")
  }

  genes
}


assign_gene_regions <- function(genes) {
  region_gr <- makeGRangesFromDataFrame(
    region_boundaries,
    seqnames.field = "chromosome",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )

  gene_gr <- GRanges(
    seqnames = genes$chromosome,
    ranges = IRanges(
      start = genes$start,
      end = genes$end
    ),
    strand = genes$strand
  )

  gene_midpoints <- resize(
    gene_gr,
    width = 1,
    fix = "center"
  )

  hits <- findOverlaps(
    gene_midpoints,
    region_gr,
    ignore.strand = TRUE
  )

  region_assignment <- rep(
    NA_character_,
    length(gene_midpoints)
  )

  region_assignment[queryHits(hits)] <- as.character(
    mcols(region_gr)$region[
      subjectHits(hits)
    ]
  )

  genes %>%
    mutate(
      region = region_assignment
    ) %>%
    filter(
      region %in% c("Eu", "PCH")
    )
}


choose_gene_match_column <- function(
    gene_annotation,
    seurat_features
) {
  match_counts <- c(
    gene_id = sum(
      gene_annotation$gene_id %in% seurat_features
    ),
    gene_id_raw = sum(
      gene_annotation$gene_id_raw %in% seurat_features
    ),
    gene_name = sum(
      gene_annotation$gene_name %in% seurat_features,
      na.rm = TRUE
    )
  )

  selected_column <- names(
    which.max(match_counts)
  )

  message(
    "Gene-feature matching counts: ",
    paste(
      names(match_counts),
      match_counts,
      sep = "=",
      collapse = ", "
    )
  )

  message(
    "Using gene annotation column: ",
    selected_column
  )

  selected_column
}


calculate_group_means <- function(
    expression_matrix,
    metadata,
    feature_names
) {
  group_table <- metadata %>%
    distinct(
      replicate,
      cell_type
    ) %>%
    arrange(
      replicate,
      cell_type
    )

  output_matrix <- matrix(
    NA_real_,
    nrow = length(feature_names),
    ncol = nrow(group_table),
    dimnames = list(
      feature_names,
      paste(
        group_table$replicate,
        group_table$cell_type,
        sep = "|||"
      )
    )
  )

  cell_counts <- integer(
    nrow(group_table)
  )

  for (i in seq_len(nrow(group_table))) {
    group_cells <- metadata$cell[
      metadata$replicate == group_table$replicate[i] &
        metadata$cell_type == group_table$cell_type[i]
    ]

    cell_counts[i] <- length(group_cells)

    output_matrix[, i] <- Matrix::rowMeans(
      expression_matrix[
        feature_names,
        group_cells,
        drop = FALSE
      ]
    )
  }

  group_table$n_cells <- cell_counts

  list(
    expression = output_matrix,
    groups = group_table
  )
}


# ------------------------------------------------------------------------------
# 6. Load and validate the Seurat object
# ------------------------------------------------------------------------------

if (!file.exists(seurat_file)) {
  stop("Seurat object not found: ", seurat_file)
}

seurat_object <- readRDS(
  seurat_file
)

annotation_candidates <- c(
  "newannotation",
  "annotation",
  "S_annotation"
)

annotation_column <- annotation_candidates[
  annotation_candidates %in%
    colnames(seurat_object[[]])
][1]

if (
  length(annotation_column) == 0 ||
    is.na(annotation_column)
) {
  stop(
    "No supported cell-type annotation column was found. Expected one of: ",
    paste(annotation_candidates, collapse = ", ")
  )
}

annotations <- as.character(
  seurat_object[[]][[annotation_column]]
)

names(annotations) <- Cells(
  seurat_object
)

valid_cells <- (
  !is.na(annotations) &
    nzchar(annotations) &
    annotations != "unannotated" &
    toupper(annotations) != "NA"
)

seurat_object <- subset(
  seurat_object,
  cells = names(annotations)[valid_cells]
)

annotations <- as.character(
  seurat_object[[]][[annotation_column]]
)

unknown_cell_types <- sort(
  setdiff(
    unique(annotations),
    known_cell_types
  )
)

if (length(unknown_cell_types) > 0) {
  stop(
    "The following cell types are not assigned to germline or soma:\n",
    paste(unknown_cell_types, collapse = "\n")
  )
}

cell_metadata <- tibble(
  cell = Cells(seurat_object),
  cell_type = annotations,
  replicate = extract_replicate(
    Cells(seurat_object)
  )
) %>%
  mutate(
    lineage = assign_lineage(cell_type)
  )


# ------------------------------------------------------------------------------
# 7. Assign genes to Eu and PCH and match them to the Seurat object
# ------------------------------------------------------------------------------

gene_annotation <- read_gene_annotation(
  gtf_file
) %>%
  assign_gene_regions()

normalized_matrix <- get_normalized_matrix(
  seurat_object
)

match_column <- choose_gene_match_column(
  gene_annotation,
  rownames(normalized_matrix)
)

matched_genes <- gene_annotation %>%
  mutate(
    feature = as.character(
      .data[[match_column]]
    )
  ) %>%
  filter(
    !is.na(feature),
    nzchar(feature),
    feature %in% rownames(normalized_matrix)
  ) %>%
  distinct(
    feature,
    .keep_all = TRUE
  )

if (nrow(matched_genes) == 0) {
  stop("No annotated genes matched the Seurat expression matrix.")
}

write.csv(
  matched_genes,
  file.path(
    output_dir,
    "matched_gene_compartment_annotation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8. Replicate-by-cell-type pseudobulk expression
# ------------------------------------------------------------------------------

pseudobulk_result <- calculate_group_means(
  expression_matrix = normalized_matrix,
  metadata = cell_metadata,
  feature_names = matched_genes$feature
)

pseudobulk_expression <- pseudobulk_result$expression
pseudobulk_groups <- pseudobulk_result$groups

write.csv(
  pseudobulk_groups,
  file.path(
    output_dir,
    "replicate_celltype_cell_counts.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 9. Z-score each gene across cell types within each replicate
# ------------------------------------------------------------------------------

pseudobulk_z <- pseudobulk_expression

for (replicate_name in unique(pseudobulk_groups$replicate)) {
  replicate_columns <- paste(
    pseudobulk_groups$replicate[
      pseudobulk_groups$replicate ==
        replicate_name
    ],
    pseudobulk_groups$cell_type[
      pseudobulk_groups$replicate ==
        replicate_name
    ],
    sep = "|||"
  )

  pseudobulk_z[, replicate_columns] <- t(
    scale(
      t(
        pseudobulk_expression[
          ,
          replicate_columns,
          drop = FALSE
        ]
      )
    )
  )
}

pseudobulk_z[
  is.na(pseudobulk_z)
] <- 0


# ------------------------------------------------------------------------------
# 10. Construct the gene-level long table
# ------------------------------------------------------------------------------

gene_expression_long <- tibble(
  feature = rep(
    rownames(pseudobulk_expression),
    times = ncol(pseudobulk_expression)
  ),
  sample_group = rep(
    colnames(pseudobulk_expression),
    each = nrow(pseudobulk_expression)
  ),
  expression = as.numeric(
    c(pseudobulk_expression)
  ),
  z_score = as.numeric(
    c(pseudobulk_z)
  )
) %>%
  separate(
    sample_group,
    into = c(
      "replicate",
      "cell_type"
    ),
    sep = "\\|\\|\\|",
    remove = FALSE
  ) %>%
  left_join(
    matched_genes %>%
      select(
        feature,
        gene_id,
        gene_id_raw,
        gene_name,
        chromosome,
        start,
        end,
        strand,
        region
      ),
    by = "feature"
  ) %>%
  mutate(
    lineage = assign_lineage(cell_type)
  )

keep_features <- gene_expression_long %>%
  group_by(feature) %>%
  summarise(
    n_expressed_groups = sum(
      expression > 0,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    n_expressed_groups >=
      minimum_expressed_groups
  ) %>%
  pull(feature)

gene_expression_long <- gene_expression_long %>%
  filter(
    feature %in% keep_features,
    region %in% c("Eu", "PCH")
  ) %>%
  mutate(
    replicate = factor(replicate),
    lineage = factor(
      lineage,
      levels = c(
        "Somatic",
        "Germline"
      )
    ),
    region = factor(
      region,
      levels = c(
        "Eu",
        "PCH"
      )
    )
  )

write.csv(
  gene_expression_long,
  file.path(
    output_dir,
    "protein_coding_gene_expression_long.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 11. Aggregate to replicate x lineage x region
# ------------------------------------------------------------------------------

gene_rep_lineage <- gene_expression_long %>%
  group_by(
    replicate,
    lineage,
    region
  ) %>%
  summarise(
    mean_z = mean(
      z_score,
      na.rm = TRUE
    ),
    median_z = median(
      z_score,
      na.rm = TRUE
    ),
    n_observations = n(),
    n_genes = n_distinct(feature),
    .groups = "drop"
  ) %>%
  arrange(
    replicate,
    lineage,
    region
  )

expected_combinations <- expand.grid(
  replicate = levels(
    gene_rep_lineage$replicate
  ),
  lineage = c(
    "Somatic",
    "Germline"
  ),
  region = c(
    "Eu",
    "PCH"
  ),
  stringsAsFactors = FALSE
)

observed_combinations <- gene_rep_lineage %>%
  transmute(
    replicate = as.character(replicate),
    lineage = as.character(lineage),
    region = as.character(region)
  )

missing_combinations <- anti_join(
  expected_combinations,
  observed_combinations,
  by = c(
    "replicate",
    "lineage",
    "region"
  )
)

if (nrow(missing_combinations) > 0) {
  stop(
    "Missing replicate-lineage-region combinations:\n",
    paste(
      apply(
        missing_combinations,
        1,
        paste,
        collapse = " | "
      ),
      collapse = "\n"
    )
  )
}

write.csv(
  gene_rep_lineage,
  file.path(
    output_dir,
    "protein_coding_gene_replicate_lineage_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 12. Fit the mixed-effects model
# ------------------------------------------------------------------------------

gene_model <- nlme::lme(
  fixed = mean_z ~ region * lineage,
  random = ~ 1 | replicate,
  data = gene_rep_lineage,
  method = "ML",
  na.action = na.omit,
  control = nlme::lmeControl(
    opt = "optim",
    maxIter = 1000,
    msMaxIter = 1000
  )
)

saveRDS(
  gene_model,
  file.path(
    output_dir,
    "protein_coding_gene_compartment_model.rds"
  )
)

anova_table <- as.data.frame(
  anova(
    gene_model,
    type = "marginal"
  )
) %>%
  rownames_to_column("term")

write.csv(
  anova_table,
  file.path(
    output_dir,
    "protein_coding_gene_compartment_ANOVA.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 13. Estimated marginal means and contrasts
# ------------------------------------------------------------------------------

region_emmeans <- emmeans(
  gene_model,
  ~ region | lineage
)

region_emmeans_table <- as.data.frame(
  summary(
    region_emmeans,
    infer = c(TRUE, TRUE)
  )
)

write.csv(
  region_emmeans_table,
  file.path(
    output_dir,
    "protein_coding_gene_region_emmeans.csv"
  ),
  row.names = FALSE
)

region_contrasts <- as.data.frame(
  summary(
    pairs(
      region_emmeans,
      adjust = "none"
    ),
    infer = c(TRUE, TRUE),
    adjust = "none"
  )
) %>%
  mutate(
    p_adjusted_BH = p.adjust(
      p.value,
      method = "BH"
    ),
    direction = case_when(
      contrast == "Eu - PCH" & estimate > 0 ~ "Eu > PCH",
      contrast == "Eu - PCH" & estimate < 0 ~ "PCH > Eu",
      TRUE ~ "No difference"
    ),
    significant_FDR_0.05 = (
      !is.na(p_adjusted_BH) &
        p_adjusted_BH < 0.05
    )
  )

write.csv(
  region_contrasts,
  file.path(
    output_dir,
    "protein_coding_gene_Eu_vs_PCH_contrasts.csv"
  ),
  row.names = FALSE
)

lineage_emmeans <- emmeans(
  gene_model,
  ~ lineage | region
)

lineage_contrasts <- as.data.frame(
  summary(
    pairs(
      lineage_emmeans,
      adjust = "none"
    ),
    infer = c(TRUE, TRUE),
    adjust = "none"
  )
) %>%
  mutate(
    p_adjusted_BH = p.adjust(
      p.value,
      method = "BH"
    ),
    direction = case_when(
      contrast == "Somatic - Germline" & estimate > 0 ~
        "Somatic > Germline",
      contrast == "Somatic - Germline" & estimate < 0 ~
        "Germline > Somatic",
      TRUE ~
        "No difference"
    ),
    significant_FDR_0.05 = (
      !is.na(p_adjusted_BH) &
        p_adjusted_BH < 0.05
    )
  )

write.csv(
  lineage_contrasts,
  file.path(
    output_dir,
    "protein_coding_gene_lineage_contrasts.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 14. Save model summary, parameters, and session information
# ------------------------------------------------------------------------------

writeLines(
  c(
    capture.output(
      summary(gene_model)
    ),
    "",
    "Marginal ANOVA:",
    capture.output(
      print(anova_table)
    ),
    "",
    "Eu versus PCH contrasts:",
    capture.output(
      print(region_contrasts)
    )
  ),
  file.path(
    output_dir,
    "protein_coding_gene_compartment_model_summary.txt"
  )
)

analysis_parameters <- data.frame(
  parameter = c(
    "Seurat object",
    "GTF annotation",
    "Gene-compartment assignment",
    "Compartments modelled",
    "Expression filter",
    "Scaling",
    "Model",
    "Random effect",
    "Multiple-testing correction"
  ),
  value = c(
    seurat_file,
    gtf_file,
    "Gene midpoint",
    "Eu and PCH",
    paste0(
      "Expression > 0 in at least ",
      minimum_expressed_groups,
      " replicate-cell-type groups"
    ),
    "Each gene across cell types within each replicate",
    "mean_z ~ region * lineage",
    "Replicate intercept",
    "Benjamini-Hochberg within each contrast set"
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
    "sessionInfo_protein_coding_gene_compartment.txt"
  )
)

message("\nProtein-coding gene compartment analysis complete.")
message("Results were written to: ", output_dir)
