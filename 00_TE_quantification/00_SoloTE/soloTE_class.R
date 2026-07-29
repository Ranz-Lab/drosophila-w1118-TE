#!/usr/bin/env Rscript

# ==============================================================================
# Process class-level SoloTE output and prepare mixed-model input
#
# This script:
#
#   1. Reads class-level SoloTE matrices for three w1118 testis replicates
#   2. Transfers Fly Cell Atlas cell-type annotations using matched barcodes
#   3. Creates and saves filtered sample-level Seurat objects
#   4. Merges, normalizes, scales, and reduces the three samples
#   5. Adds broad spermatogenesis-stage annotations
#   6. Calculates replicate-by-cell-type mean scaled expression for:
#        DNA, LINE, and LTR
#   7. Writes the input table used by the TE-class mixed-effects analysis
#
# Required input:
#   FCA_testis_10x.rds
#   w1118S1_LBAM2_classtes_MATRIX/
#   w1118S2_LBAM_classtes_MATRIX/
#   w1118S3_LBAM2_classtes_MATRIX/
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(tibble)
  library(tidyr)
})


# ------------------------------------------------------------------------------
# 2. Inputs, outputs, and parameters
# ------------------------------------------------------------------------------

FCA_reference_file <- "FCA_testis_10x.rds"

sample_table <- tibble(
  sample_name = c(
    "sample1",
    "sample2",
    "sample3"
  ),
  replicate = c(
    "rep1",
    "rep2",
    "rep3"
  ),
  replicate_suffix = c(
    "1",
    "2",
    "3"
  ),
  FCA_id = c(
    "FCA59",
    "FCA60",
    "FCA61"
  ),
  matrix_directory = c(
    "w1118S1_LBAM2_classtes_MATRIX",
    "w1118S2_LBAM_classtes_MATRIX",
    "w1118S3_LBAM2_classtes_MATRIX"
  ),
  filtered_object_file = c(
    "Male_testis_sample1_filtered_FCA_soloTE_class.rds",
    "Male_testis_sample2_filtered_FCA_soloTE_class.rds",
    "Male_testis_sample3_filtered_FCA_soloTE_class.rds"
  )
)

normalized_object_file <-
  "FCA-testis-merged-normalized-scaled-soloTE-class-w1118.rds"

final_object_file <-
  "FCA-testis-merged-umap-final-soloTE-class-w1118.rds"

mixed_model_input_file <-
  "INPUT_LME_replicate_celltype_scaled_TEclass.csv"

output_dir <- "SoloTE_class_processing_results"

minimum_cells_feature <- 5L
minimum_cell_features <- 200L
variable_feature_count <- 2000L
minimum_cells_per_replicate_celltype <- 20L

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 3. Cell-type definitions
# ------------------------------------------------------------------------------

somatic_cyst_cell_types <- c(
  "cyst stem cell",
  "early cyst cell 1",
  "early cyst cell 2",
  "cyst cell intermediate",
  "cyst cell branch a",
  "cyst cell branch b",
  "late cyst cell branch a",
  "late cyst cell branch b",
  "spermatocyte cyst cell branch a",
  "spermatocyte cyst cell branch b",
  "head cyst cell"
)

somatic_other_cell_types <- c(
  "secretory cell of the male reproductive tract",
  "pigment cell",
  "male gonad associated epithelium",
  "muscle cell",
  "hemocyte",
  "adult fat body",
  "adult tracheal cell"
)

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

somatic_cell_types <- c(
  somatic_cyst_cell_types,
  somatic_other_cell_types
)

known_cell_types <- c(
  somatic_cell_types,
  germline_cell_types
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


# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

barcode_core <- function(cell_names) {
  sub(
    "-.*$",
    "",
    cell_names
  )
}


assign_compartment <- function(cell_type) {
  case_when(
    cell_type %in% somatic_cell_types ~ "Somatic",
    cell_type %in% germline_cell_types ~ "Germline",
    TRUE ~ NA_character_
  )
}


assign_broad_stage <- function(cell_type) {
  case_when(
    cell_type %in% somatic_cyst_cell_types ~
      "Somatic - Cyst",
    cell_type %in% somatic_other_cell_types ~
      "Somatic - Other",
    cell_type %in% mitotic_cell_types ~
      "Germline - Mitotic",
    cell_type %in% meiotic_cell_types ~
      "Germline - Meiotic",
    cell_type %in% postmeiotic_cell_types ~
      "Germline - Post-meiotic",
    TRUE ~
      NA_character_
  )
}


get_scaled_matrix <- function(seurat_object) {
  if (inherits(seurat_object[["RNA"]], "Assay5")) {
    return(
      GetAssayData(
        seurat_object,
        assay = "RNA",
        layer = "scale.data"
      )
    )
  }

  GetAssayData(
    seurat_object,
    assay = "RNA",
    slot = "scale.data"
  )
}


process_sample <- function(
    matrix_directory,
    FCA_reference,
    FCA_id,
    sample_name,
    replicate,
    replicate_suffix,
    output_file
) {
  if (!dir.exists(matrix_directory)) {
    stop(
      "SoloTE matrix directory not found: ",
      matrix_directory
    )
  }

  FCA_sample <- subset(
    FCA_reference,
    subset = fca_id == FCA_id
  )

  if (ncol(FCA_sample) == 0) {
    stop(
      "No cells were found for FCA ID: ",
      FCA_id
    )
  }

  counts <- Read10X(
    data.dir = matrix_directory,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
  )

  sample_object <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = minimum_cells_feature,
    min.features = minimum_cell_features
  )

  FCA_cores <- barcode_core(
    Cells(FCA_sample)
  )

  target_cores <- barcode_core(
    Cells(sample_object)
  )

  if (anyDuplicated(FCA_cores)) {
    stop(
      "Duplicated barcode cores were found in FCA sample ",
      FCA_id,
      "."
    )
  }

  if (anyDuplicated(target_cores)) {
    stop(
      "Duplicated barcode cores were found in ",
      matrix_directory,
      "."
    )
  }

  annotation_map <- setNames(
    as.character(
      FCA_sample$S_annotation
    ),
    FCA_cores
  )

  matched_cells <- Cells(sample_object)[
    target_cores %in%
      names(annotation_map)
  ]

  if (length(matched_cells) == 0) {
    stop(
      "No matched barcodes were found for ",
      sample_name,
      "."
    )
  }

  sample_object <- subset(
    sample_object,
    cells = matched_cells
  )

  matched_cores <- barcode_core(
    Cells(sample_object)
  )

  transferred_annotation <- annotation_map[
    matched_cores
  ]

  keep_cells <- (
    !is.na(transferred_annotation) &
      nzchar(transferred_annotation) &
      toupper(transferred_annotation) != "NA" &
      transferred_annotation != "unannotated"
  )

  sample_object <- subset(
    sample_object,
    cells = Cells(sample_object)[keep_cells]
  )

  transferred_annotation <- transferred_annotation[
    keep_cells
  ]

  sample_object$annotation <-
    unname(transferred_annotation)

  sample_object$newannotation <-
    unname(transferred_annotation)

  sample_object$replicate <- replicate
  sample_object$FCA_id <- FCA_id

  new_cell_names <- paste0(
    Cells(sample_object),
    "_",
    replicate_suffix
  )

  sample_object <- RenameCells(
    sample_object,
    new.names = new_cell_names
  )

  saveRDS(
    sample_object,
    file.path(
      output_dir,
      output_file
    )
  )

  data.frame(
    sample = sample_name,
    replicate = replicate,
    FCA_id = FCA_id,
    n_raw_cells = ncol(counts),
    n_matched_annotated_cells = ncol(sample_object),
    stringsAsFactors = FALSE
  ) -> sample_summary

  list(
    object = sample_object,
    summary = sample_summary
  )
}


# ------------------------------------------------------------------------------
# 5. Load the FCA reference and process all three samples
# ------------------------------------------------------------------------------

if (!file.exists(FCA_reference_file)) {
  stop(
    "FCA reference object not found: ",
    FCA_reference_file
  )
}

FCA_reference <- readRDS(
  FCA_reference_file
)

required_FCA_columns <- c(
  "fca_id",
  "S_annotation"
)

missing_FCA_columns <- setdiff(
  required_FCA_columns,
  colnames(
    FCA_reference[[]]
  )
)

if (length(missing_FCA_columns) > 0) {
  stop(
    "The FCA object is missing required metadata columns: ",
    paste(
      missing_FCA_columns,
      collapse = ", "
    )
  )
}

processed_samples <- vector(
  mode = "list",
  length = nrow(sample_table)
)

sample_summaries <- vector(
  mode = "list",
  length = nrow(sample_table)
)

for (i in seq_len(nrow(sample_table))) {
  sample_result <- process_sample(
    matrix_directory =
      sample_table$matrix_directory[i],
    FCA_reference = FCA_reference,
    FCA_id = sample_table$FCA_id[i],
    sample_name =
      sample_table$sample_name[i],
    replicate =
      sample_table$replicate[i],
    replicate_suffix =
      sample_table$replicate_suffix[i],
    output_file =
      sample_table$filtered_object_file[i]
  )

  processed_samples[[i]] <-
    sample_result$object

  sample_summaries[[i]] <-
    sample_result$summary
}

sample_summary_table <- bind_rows(
  sample_summaries
)

write.csv(
  sample_summary_table,
  file.path(
    output_dir,
    "class_level_annotation_transfer_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 6. Merge, normalize, scale, and reduce the samples
# ------------------------------------------------------------------------------

Control.samples <- merge(
  x = processed_samples[[1]],
  y = processed_samples[2:3],
  project = "TE.Expression.Class",
  merge.data = TRUE
)

if (
  "JoinLayers" %in%
    getNamespaceExports(
      "SeuratObject"
    )
) {
  Control.samples <- JoinLayers(
    Control.samples
  )
}

Control.samples <- NormalizeData(
  Control.samples,
  verbose = FALSE
)

Control.samples <- FindVariableFeatures(
  Control.samples,
  selection.method = "vst",
  nfeatures = variable_feature_count,
  verbose = FALSE
)

Control.samples <- ScaleData(
  Control.samples,
  features = rownames(Control.samples),
  verbose = FALSE
)

saveRDS(
  Control.samples,
  file.path(
    output_dir,
    normalized_object_file
  )
)

Control.samples <- RunPCA(
  Control.samples,
  features = VariableFeatures(
    Control.samples
  ),
  npcs = 30,
  verbose = FALSE
)

Control.samples <- RunUMAP(
  Control.samples,
  dims = 1:30,
  verbose = FALSE
)

Control.samples <- FindNeighbors(
  Control.samples,
  dims = 1:30,
  verbose = FALSE
)

Control.samples <- FindClusters(
  Control.samples,
  resolution = 0.5,
  verbose = FALSE
)


# ------------------------------------------------------------------------------
# 7. Add broad-stage and compartment metadata
# ------------------------------------------------------------------------------

annotations <- as.character(
  Control.samples$newannotation
)

unknown_cell_types <- sort(
  setdiff(
    unique(annotations),
    known_cell_types
  )
)

if (length(unknown_cell_types) > 0) {
  stop(
    "The following cell types are not assigned to a broad category:\n",
    paste(
      unknown_cell_types,
      collapse = "\n"
    )
  )
}

Control.samples$compartment <-
  assign_compartment(
    annotations
  )

Control.samples$cell_stage <-
  assign_broad_stage(
    annotations
  )

Control.samples$newannotation <- factor(
  annotations,
  levels = cell_type_order
)

Control.samples$compartment <- factor(
  Control.samples$compartment,
  levels = c(
    "Somatic",
    "Germline"
  )
)

Control.samples$cell_stage <- factor(
  Control.samples$cell_stage,
  levels = c(
    "Somatic - Cyst",
    "Somatic - Other",
    "Germline - Mitotic",
    "Germline - Meiotic",
    "Germline - Post-meiotic"
  )
)

saveRDS(
  Control.samples,
  file.path(
    output_dir,
    final_object_file
  )
)

write.csv(
  as.data.frame(
    table(
      Control.samples$replicate,
      Control.samples$newannotation
    )
  ) %>%
    rename(
      replicate = Var1,
      cell_type = Var2,
      n_cells = Freq
    ) %>%
    filter(
      n_cells > 0
    ),
  file.path(
    output_dir,
    "replicate_celltype_cell_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    table(
      Control.samples$cell_stage
    )
  ) %>%
    rename(
      cell_stage = Var1,
      n_cells = Freq
    ),
  file.path(
    output_dir,
    "broad_cell_stage_counts.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8. Prepare replicate-level TE-class mixed-model input
# ------------------------------------------------------------------------------

TE_class_features <- c(
  "SoloTE-DNA",
  "SoloTE-LINE",
  "SoloTE-LTR"
)

scaled_matrix <- get_scaled_matrix(
  Control.samples
)

missing_TE_class_features <- setdiff(
  TE_class_features,
  rownames(
    scaled_matrix
  )
)

if (length(missing_TE_class_features) > 0) {
  stop(
    "The following class-level SoloTE features are absent from scale.data: ",
    paste(
      missing_TE_class_features,
      collapse = ", "
    )
  )
}

cell_metadata <- tibble(
  cell = Cells(Control.samples),
  replicate = as.character(
    Control.samples$replicate
  ),
  newannotation = as.character(
    Control.samples$newannotation
  ),
  compartment = as.character(
    Control.samples$compartment
  )
)

group_counts <- cell_metadata %>%
  count(
    replicate,
    newannotation,
    compartment,
    name = "n_cells"
  )

replicate_celltype_expression <- lapply(
  seq_len(
    nrow(group_counts)
  ),
  function(i) {
    group_cells <- cell_metadata$cell[
      cell_metadata$replicate ==
        group_counts$replicate[i] &
        cell_metadata$newannotation ==
          group_counts$newannotation[i]
    ]

    group_means <- Matrix::rowMeans(
      scaled_matrix[
        TE_class_features,
        group_cells,
        drop = FALSE
      ]
    )

    data.frame(
      replicate =
        group_counts$replicate[i],
      newannotation =
        group_counts$newannotation[i],
      compartment =
        group_counts$compartment[i],
      n_cells =
        group_counts$n_cells[i],
      TE_class = sub(
        "^SoloTE-",
        "",
        names(group_means)
      ),
      mean_expr = as.numeric(
        group_means
      ),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows() %>%
  filter(
    n_cells >=
      minimum_cells_per_replicate_celltype
  ) %>%
  mutate(
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
    TE_class = factor(
      TE_class,
      levels = c(
        "DNA",
        "LINE",
        "LTR"
      )
    )
  ) %>%
  arrange(
    replicate,
    newannotation,
    TE_class
  )

write.csv(
  replicate_celltype_expression,
  file.path(
    output_dir,
    mixed_model_input_file
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 9. Save parameters and session information
# ------------------------------------------------------------------------------

analysis_parameters <- data.frame(
  parameter = c(
    "FCA reference",
    "Sample 1 matrix",
    "Sample 2 matrix",
    "Sample 3 matrix",
    "FCA sample IDs",
    "CreateSeuratObject min.cells",
    "CreateSeuratObject min.features",
    "Variable features",
    "UMAP dimensions",
    "Clustering resolution",
    "Minimum cells per replicate-cell type",
    "TE-class expression source"
  ),
  value = c(
    FCA_reference_file,
    sample_table$matrix_directory[1],
    sample_table$matrix_directory[2],
    sample_table$matrix_directory[3],
    paste(
      sample_table$FCA_id,
      collapse = ", "
    ),
    minimum_cells_feature,
    minimum_cell_features,
    variable_feature_count,
    "1:30",
    0.5,
    minimum_cells_per_replicate_celltype,
    "RNA scale.data"
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
    "sessionInfo_SoloTE_class_processing.txt"
  )
)

message("\nClass-level SoloTE processing complete.")
message("Results were written to: ", output_dir)
