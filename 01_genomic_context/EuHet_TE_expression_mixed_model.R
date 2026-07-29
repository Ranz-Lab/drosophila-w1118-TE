#!/usr/bin/env Rscript

# ==============================================================================
# Euchromatic-versus-heterochromatic TE expression mixed models
#
# This self-contained script:
#
#   1. Loads the processed locus-level SoloTE Seurat object
#   2. Parses genomic coordinates and TE annotations from SoloTE feature names
#   3. Assigns each locus to euchromatin (EU) or heterochromatin (HET)
#   4. Groups loci as DNA transposons or retrotransposons (LINE + LTR)
#   5. Calculates mean scaled expression for EU and HET loci within each cell
#   6. Aggregates expression by replicate, cell type, and partition
#   7. Fits cell-type and compartment mixed-effects models
#   8. Exports Type III ANOVA tables and globally BH-adjusted EU-HET contrasts
#
# Required input:
#   FCA-testis-merged-umap-final-soloTE-locus-w1118.rds
#
# Chromosome 4 and Y loci are classified as HET, matching the original input
# preparation used for this analysis.
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
  library(lme4)
  library(lmerTest)
  library(Matrix)
  library(Seurat)
  library(stringr)
  library(tibble)
  library(tidyr)
})

options(contrasts = c("contr.sum", "contr.poly"))


# ------------------------------------------------------------------------------
# 2. Input, output, and parameters
# ------------------------------------------------------------------------------

seurat_file <- "FCA-testis-merged-umap-final-soloTE-locus-w1118.rds"

output_dir <- "EuHet_TE_mixed_model_results"

minimum_cells_per_replicate_celltype <- 20L
include_chr4_Y_as_HET <- TRUE

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

cell_type_order <- c(
  "secretory cell of the male reproductive tract",
  "mid-late elongation-stage spermatid",
  "early-mid elongation-stage spermatid",
  "early elongation stage spermatid",
  "spermatid",
  "late primary spermatocyte",
  "spermatocyte",
  "spermatocyte 7a",
  "spermatocyte 6",
  "spermatocyte 5",
  "spermatocyte 4",
  "spermatocyte 3",
  "spermatocyte 2",
  "spermatocyte 1",
  "spermatocyte 0",
  "spermatogonium-spermatocyte transition",
  "mid-late proliferating spermatogonia",
  "spermatogonium",
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
  germline_cell_types,
  somatic_cell_types
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

  factor(
    paste0("rep", replicate_number),
    levels = c(
      "rep1",
      "rep2",
      "rep3"
    )
  )
}


assign_compartment <- function(cell_type) {
  case_when(
    cell_type %in% somatic_cell_types ~ "Somatic",
    cell_type %in% germline_cell_types ~ "Germline",
    TRUE ~ NA_character_
  )
}


get_scaled_matrix <- function(seurat_object) {
  if (inherits(seurat_object[["RNA"]], "Assay5")) {
    scaled_matrix <- GetAssayData(
      seurat_object,
      assay = "RNA",
      layer = "scale.data"
    )
  } else {
    scaled_matrix <- GetAssayData(
      seurat_object,
      assay = "RNA",
      slot = "scale.data"
    )
  }

  if (
    is.null(scaled_matrix) ||
      nrow(scaled_matrix) == 0
  ) {
    stop(
      "The RNA scale.data matrix is empty. ",
      "The locus-level object must be scaled before this analysis."
    )
  }

  scaled_matrix
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

  data.frame(
    feature = features,
    chromosome = feature_match[, 2],
    start = suppressWarnings(
      as.integer(feature_match[, 3])
    ),
    end = suppressWarnings(
      as.integer(feature_match[, 4])
    ),
    TE_family = feature_match[, 5],
    TE_superfamily = feature_match[, 6],
    TE_class = toupper(
      feature_match[, 7]
    ),
    divergence = suppressWarnings(
      as.numeric(feature_match[, 8])
    ),
    strand = feature_match[, 9],
    stringsAsFactors = FALSE
  ) %>%
    filter(
      !is.na(chromosome),
      !is.na(start),
      !is.na(end)
    ) %>%
    mutate(
      midpoint = floor(
        (start + end) / 2
      )
    )
}


assign_locus_regions <- function(
    TE_metadata,
    include_chr4_Y
) {
  region_table <- region_boundaries

  if (include_chr4_Y) {
    extra_regions <- TE_metadata %>%
      filter(
        chromosome %in% c(
          "4",
          "Y"
        )
      ) %>%
      group_by(chromosome) %>%
      summarise(
        start = 1L,
        end = max(
          end,
          na.rm = TRUE
        ),
        region = "HET_extra",
        .groups = "drop"
      )

    region_table <- bind_rows(
      region_table,
      extra_regions
    )
  }

  region_gr <- makeGRangesFromDataFrame(
    region_table,
    seqnames.field = "chromosome",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )

  TE_gr <- GRanges(
    seqnames = TE_metadata$chromosome,
    ranges = IRanges(
      start = TE_metadata$midpoint,
      width = 1
    )
  )

  hits <- findOverlaps(
    TE_gr,
    region_gr,
    ignore.strand = TRUE
  )

  region_assignment <- rep(
    NA_character_,
    nrow(TE_metadata)
  )

  region_assignment[queryHits(hits)] <- as.character(
    mcols(region_gr)$region[
      subjectHits(hits)
    ]
  )

  TE_metadata %>%
    mutate(
      region = region_assignment,
      partition = case_when(
        region == "Eu" ~ "EU",
        region %in% c(
          "PCH",
          "Tel",
          "HET_extra"
        ) ~ "HET",
        TRUE ~ NA_character_
      ),
      TE_group = case_when(
        TE_class == "DNA" ~ "DNA",
        TE_class %in% c(
          "LINE",
          "LTR"
        ) ~ "Retro",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(partition),
      !is.na(TE_group)
    ) %>%
    distinct(
      feature,
      .keep_all = TRUE
    )
}


make_model_input <- function(
    scaled_matrix,
    TE_metadata,
    cell_metadata,
    TE_group_name,
    minimum_cells
) {
  metadata_subset <- TE_metadata %>%
    filter(
      TE_group ==
        TE_group_name
    )

  EU_features <- metadata_subset$feature[
    metadata_subset$partition ==
      "EU"
  ]

  HET_features <- metadata_subset$feature[
    metadata_subset$partition ==
      "HET"
  ]

  if (
    length(EU_features) == 0 ||
      length(HET_features) == 0
  ) {
    stop(
      TE_group_name,
      " lacks loci in one of the EU/HET partitions."
    )
  }

  EU_expression <- Matrix::colMeans(
    scaled_matrix[
      EU_features,
      ,
      drop = FALSE
    ]
  )

  HET_expression <- Matrix::colMeans(
    scaled_matrix[
      HET_features,
      ,
      drop = FALSE
    ]
  )

  cell_expression <- cell_metadata %>%
    mutate(
      EU = as.numeric(
        EU_expression[cell]
      ),
      HET = as.numeric(
        HET_expression[cell]
      )
    )

  cell_counts <- cell_expression %>%
    count(
      replicate,
      newannotation,
      compartment,
      name = "n_cells"
    )

  model_input <- cell_expression %>%
    pivot_longer(
      cols = c(
        EU,
        HET
      ),
      names_to = "partition",
      values_to = "scaled_expression"
    ) %>%
    group_by(
      replicate,
      newannotation,
      compartment,
      partition
    ) %>%
    summarise(
      mean_scaled_expr = mean(
        scaled_expression,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    left_join(
      cell_counts,
      by = c(
        "replicate",
        "newannotation",
        "compartment"
      )
    ) %>%
    filter(
      n_cells >= minimum_cells
    ) %>%
    mutate(
      TE_group = TE_group_name,
      n_loci = case_when(
        partition == "EU" ~
          length(EU_features),
        partition == "HET" ~
          length(HET_features),
        TRUE ~
          NA_integer_
      ),
      replicate = factor(
        replicate,
        levels = c(
          "rep1",
          "rep2",
          "rep3"
        )
      ),
      newannotation = factor(
        newannotation,
        levels = cell_type_order
      ),
      compartment = factor(
        compartment,
        levels = c(
          "Somatic",
          "Germline"
        )
      ),
      partition = factor(
        partition,
        levels = c(
          "EU",
          "HET"
        )
      )
    ) %>%
    arrange(
      replicate,
      newannotation,
      partition
    )

  model_input
}


format_anova <- function(anova_result) {
  as.data.frame(
    anova_result
  ) %>%
    rownames_to_column("term")
}


add_partition_interpretation <- function(
    contrast_table
) {
  contrast_table %>%
    mutate(
      p_adjusted_BH = p.adjust(
        p.value,
        method = "BH"
      ),
      direction = case_when(
        contrast == "EU - HET" &
          estimate > 0 ~
          "EU > HET",
        contrast == "EU - HET" &
          estimate < 0 ~
          "HET > EU",
        TRUE ~
          "No difference"
      ),
      significant_FDR_0.05 = (
        !is.na(p_adjusted_BH) &
          p_adjusted_BH < 0.05
      )
    )
}


write_model_diagnostics <- function(
    model,
    output_file,
    model_name
) {
  convergence_messages <-
    model@optinfo$conv$lme4$messages

  convergence_text <- if (
    is.null(convergence_messages)
  ) {
    "None"
  } else {
    paste(
      convergence_messages,
      collapse = "; "
    )
  }

  diagnostic_lines <- c(
    paste0(
      "Model: ",
      model_name
    ),
    "",
    paste0(
      "Number of observations: ",
      nobs(model)
    ),
    paste0(
      "Singular fit: ",
      isSingular(
        model,
        tol = 1e-4
      )
    ),
    paste0(
      "Convergence messages: ",
      convergence_text
    ),
    "",
    "Model summary:",
    capture.output(
      summary(model)
    )
  )

  writeLines(
    diagnostic_lines,
    output_file
  )
}


fit_celltype_model <- function(
    model_input,
    TE_group_name
) {
  model <- lmer(
    mean_scaled_expr ~
      partition * newannotation +
      (1 | replicate) +
      (1 | replicate:newannotation),
    data = model_input,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 100000
      )
    )
  )

  prefix <- file.path(
    output_dir,
    paste0(
      TE_group_name,
      "_celltype"
    )
  )

  saveRDS(
    model,
    paste0(
      prefix,
      "_model.rds"
    )
  )

  write_model_diagnostics(
    model = model,
    output_file = paste0(
      prefix,
      "_model_summary.txt"
    ),
    model_name = paste(
      TE_group_name,
      "EU-HET cell-type model"
    )
  )

  anova_table <- format_anova(
    anova(
      model,
      type = 3,
      ddf = "Satterthwaite"
    )
  )

  write.csv(
    anova_table,
    paste0(
      prefix,
      "_type3_ANOVA.csv"
    ),
    row.names = FALSE
  )

  estimated_means <- emmeans(
    model,
    ~ partition | newannotation,
    lmer.df = "satterthwaite"
  )

  write.csv(
    as.data.frame(
      summary(
        estimated_means,
        infer = c(
          TRUE,
          TRUE
        )
      )
    ),
    paste0(
      prefix,
      "_partition_emmeans.csv"
    ),
    row.names = FALSE
  )

  contrasts <- as.data.frame(
    summary(
      pairs(
        estimated_means,
        adjust = "none"
      ),
      infer = c(
        TRUE,
        TRUE
      ),
      adjust = "none"
    )
  ) %>%
    add_partition_interpretation()

  write.csv(
    contrasts,
    paste0(
      prefix,
      "_EU_vs_HET_globalFDR.csv"
    ),
    row.names = FALSE
  )

  invisible(model)
}


fit_compartment_model <- function(
    model_input,
    TE_group_name
) {
  model <- lmer(
    mean_scaled_expr ~
      partition * compartment +
      (1 | replicate) +
      (1 | newannotation),
    data = model_input,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 100000
      )
    )
  )

  prefix <- file.path(
    output_dir,
    paste0(
      TE_group_name,
      "_compartment"
    )
  )

  saveRDS(
    model,
    paste0(
      prefix,
      "_model.rds"
    )
  )

  write_model_diagnostics(
    model = model,
    output_file = paste0(
      prefix,
      "_model_summary.txt"
    ),
    model_name = paste(
      TE_group_name,
      "EU-HET compartment model"
    )
  )

  anova_table <- format_anova(
    anova(
      model,
      type = 3,
      ddf = "Satterthwaite"
    )
  )

  write.csv(
    anova_table,
    paste0(
      prefix,
      "_type3_ANOVA.csv"
    ),
    row.names = FALSE
  )

  estimated_means <- emmeans(
    model,
    ~ partition | compartment,
    lmer.df = "satterthwaite"
  )

  write.csv(
    as.data.frame(
      summary(
        estimated_means,
        infer = c(
          TRUE,
          TRUE
        )
      )
    ),
    paste0(
      prefix,
      "_partition_emmeans.csv"
    ),
    row.names = FALSE
  )

  contrasts <- as.data.frame(
    summary(
      pairs(
        estimated_means,
        adjust = "none"
      ),
      infer = c(
        TRUE,
        TRUE
      ),
      adjust = "none"
    )
  ) %>%
    add_partition_interpretation()

  write.csv(
    contrasts,
    paste0(
      prefix,
      "_EU_vs_HET_globalFDR.csv"
    ),
    row.names = FALSE
  )

  invisible(model)
}


# ------------------------------------------------------------------------------
# 6. Load and validate the locus-level Seurat object
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
  newannotation = annotations,
  replicate = extract_replicate(
    Cells(seurat_object)
  )
) %>%
  mutate(
    compartment = assign_compartment(
      newannotation
    )
  )


# ------------------------------------------------------------------------------
# 7. Parse TE loci and assign EU/HET partitions
# ------------------------------------------------------------------------------

scaled_matrix <- get_scaled_matrix(
  seurat_object
)

soloTE_features <- rownames(
  scaled_matrix
)[
  grepl(
    "^SoloTE[-_]",
    rownames(
      scaled_matrix
    )
  )
]

if (length(soloTE_features) == 0) {
  stop("No locus-level SoloTE features were found in scale.data.")
}

TE_metadata <- parse_soloTE_features(
  soloTE_features
)

if (nrow(TE_metadata) == 0) {
  stop("No SoloTE locus identifiers could be parsed.")
}

TE_metadata <- assign_locus_regions(
  TE_metadata = TE_metadata,
  include_chr4_Y =
    include_chr4_Y_as_HET
) %>%
  filter(
    feature %in%
      rownames(
        scaled_matrix
      )
  )

if (nrow(TE_metadata) == 0) {
  stop("No parsed TE loci remained after genomic assignment.")
}

scaled_matrix <- scaled_matrix[
  TE_metadata$feature,
  ,
  drop = FALSE
]

write.csv(
  TE_metadata,
  file.path(
    output_dir,
    "TE_locus_EU_HET_annotation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8. Build model inputs for DNA and retrotransposons
# ------------------------------------------------------------------------------

DNA_input <- make_model_input(
  scaled_matrix = scaled_matrix,
  TE_metadata = TE_metadata,
  cell_metadata = cell_metadata,
  TE_group_name = "DNA",
  minimum_cells =
    minimum_cells_per_replicate_celltype
)

Retro_input <- make_model_input(
  scaled_matrix = scaled_matrix,
  TE_metadata = TE_metadata,
  cell_metadata = cell_metadata,
  TE_group_name = "Retro",
  minimum_cells =
    minimum_cells_per_replicate_celltype
)

write.csv(
  DNA_input,
  file.path(
    output_dir,
    "INPUT_LME_EuHet_DNA_locusTE_replicate_celltype.csv"
  ),
  row.names = FALSE
)

write.csv(
  Retro_input,
  file.path(
    output_dir,
    "INPUT_LME_EuHet_Retro_locusTE_replicate_celltype.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 9. Fit the DNA and retrotransposon models
# ------------------------------------------------------------------------------

fit_celltype_model(
  model_input = DNA_input,
  TE_group_name = "DNA"
)

fit_compartment_model(
  model_input = DNA_input,
  TE_group_name = "DNA"
)

fit_celltype_model(
  model_input = Retro_input,
  TE_group_name = "Retro"
)

fit_compartment_model(
  model_input = Retro_input,
  TE_group_name = "Retro"
)


# ------------------------------------------------------------------------------
# 10. Save parameters and session information
# ------------------------------------------------------------------------------

analysis_parameters <- data.frame(
  parameter = c(
    "Seurat object",
    "TE groups",
    "Retrotransposon definition",
    "EU definition",
    "HET definition",
    "Chromosome 4 and Y handling",
    "Minimum cells per replicate-cell type",
    "Cell-type model",
    "Compartment model",
    "Multiple-testing correction"
  ),
  value = c(
    seurat_file,
    "DNA and Retro",
    "LINE + LTR",
    "Eu compartment",
    "PCH + Tel + chromosome 4 + chromosome Y",
    include_chr4_Y_as_HET,
    minimum_cells_per_replicate_celltype,
    paste0(
      "mean_scaled_expr ~ partition * newannotation + ",
      "(1 | replicate) + (1 | replicate:newannotation)"
    ),
    paste0(
      "mean_scaled_expr ~ partition * compartment + ",
      "(1 | replicate) + (1 | newannotation)"
    ),
    "Benjamini-Hochberg across all contrasts within each model"
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
    "sessionInfo_EuHet_TE_mixed_models.txt"
  )
)

message("\nEU-HET TE mixed-model analysis complete.")
message("Results were written to: ", output_dir)
