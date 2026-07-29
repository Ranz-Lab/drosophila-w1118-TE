#!/usr/bin/env Rscript

# ==============================================================================
# Process subfamily-level SoloTE output
#
# This script:
#
#   1. Reads subfamily-level SoloTE matrices for three w1118 testis replicates
#   2. Transfers Fly Cell Atlas cell-type annotations using matched barcodes
#   3. Saves the filtered sample-level Seurat objects
#   4. Merges, normalizes, scales, and reduces the three samples
#   5. Saves the merged subfamily-level SoloTE Seurat objects
#   6. Exports average log-normalized TE-subfamily expression by cell type
#
# Required input:
#   FCA_testis_10x.rds
#   w1118S1_LBAM2_subfamilytes_MATRIX/
#   w1118S2_LBAM_subfamilytes_MATRIX/
#   w1118S3_LBAM2_subfamilytes_MATRIX/
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
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
    "w1118S1_LBAM2_subfamilytes_MATRIX",
    "w1118S2_LBAM_subfamilytes_MATRIX",
    "w1118S3_LBAM2_subfamilytes_MATRIX"
  ),
  filtered_object_file = c(
    "Male_testis_sample1_filtered_FCA_soloTE_subfamily.rds",
    "Male_testis_sample2_filtered_FCA_soloTE_subfamily.rds",
    "Male_testis_sample3_filtered_FCA_soloTE_subfamily.rds"
  )
)

normalized_object_file <-
  "FCA-testis-merged-normalized-scaled-soloTE-subfamily-w1118.rds"

final_object_file <-
  "FCA-testis-merged-umap-final-soloTE-subfamily-w1118.rds"

average_expression_file <-
  "SoloTE_subfamily_avg_logNorm_byCellType.csv"

average_expression_long_file <-
  "SoloTE_subfamily_avg_logNorm_byCellType_LONG.csv"

output_dir <- "SoloTE_subfamily_processing_results"

minimum_cells_feature <- 5L
minimum_cell_features <- 200L
variable_feature_count <- 2000L
pca_dimensions <- 30L
clustering_resolution <- 0.5

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 3. Cell-type order
# ------------------------------------------------------------------------------

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


average_expression_by_cell_type <- function(
    seurat_object,
    annotation_column
) {
  tryCatch(
    AverageExpression(
      object = seurat_object,
      assays = "RNA",
      group.by = annotation_column,
      layer = "data",
      verbose = FALSE
    )$RNA,
    error = function(e) {
      AverageExpression(
        object = seurat_object,
        assays = "RNA",
        group.by = annotation_column,
        slot = "data",
        verbose = FALSE
      )$RNA
    }
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

  if (is.list(counts)) {
    if (length(counts) != 1) {
      stop(
        "Read10X returned multiple matrices for ",
        matrix_directory,
        "."
      )
    }

    counts <- counts[[1]]
  }

  raw_cell_count <- ncol(counts)

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
      "No matched FCA barcodes were found for ",
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

  transferred_annotation <- unname(
    annotation_map[
      matched_cores
    ]
  )

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

  if (ncol(sample_object) == 0) {
    stop(
      "No annotated cells remained for ",
      sample_name,
      "."
    )
  }

  sample_object$annotation <-
    transferred_annotation

  sample_object$newannotation <-
    transferred_annotation

  sample_object$replicate <-
    replicate

  sample_object$FCA_id <-
    FCA_id

  sample_object <- RenameCells(
    sample_object,
    new.names = paste0(
      Cells(sample_object),
      "_",
      replicate_suffix
    )
  )

  saveRDS(
    sample_object,
    file.path(
      output_dir,
      output_file
    )
  )

  list(
    object = sample_object,
    summary = data.frame(
      sample = sample_name,
      replicate = replicate,
      FCA_id = FCA_id,
      matrix_directory = matrix_directory,
      n_raw_cells = raw_cell_count,
      n_matched_annotated_cells =
        ncol(sample_object),
      stringsAsFactors = FALSE
    )
  )
}


# ------------------------------------------------------------------------------
# 5. Load the FCA reference
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
    "The FCA reference is missing required metadata columns: ",
    paste(
      missing_FCA_columns,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# 6. Process the three subfamily-level SoloTE matrices
# ------------------------------------------------------------------------------

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
    FCA_reference =
      FCA_reference,
    FCA_id =
      sample_table$FCA_id[i],
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
    "subfamily_annotation_transfer_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. Merge, normalize, scale, and reduce the samples
# ------------------------------------------------------------------------------

Control.samples <- merge(
  x = processed_samples[[1]],
  y = processed_samples[2:3],
  project = "TE.Expression.Subfamily",
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

Control.samples$newannotation <- factor(
  as.character(
    Control.samples$newannotation
  ),
  levels = cell_type_order
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
  npcs = pca_dimensions,
  verbose = FALSE
)

Control.samples <- RunUMAP(
  Control.samples,
  dims = seq_len(
    pca_dimensions
  ),
  verbose = FALSE
)

Control.samples <- FindNeighbors(
  Control.samples,
  dims = seq_len(
    pca_dimensions
  ),
  verbose = FALSE
)

Control.samples <- FindClusters(
  Control.samples,
  resolution = clustering_resolution,
  verbose = FALSE
)

saveRDS(
  Control.samples,
  file.path(
    output_dir,
    final_object_file
  )
)


# ------------------------------------------------------------------------------
# 8. Export average log-normalized subfamily expression by cell type
# ------------------------------------------------------------------------------

average_expression <- average_expression_by_cell_type(
  seurat_object = Control.samples,
  annotation_column = "newannotation"
)

subfamily_features <- rownames(
  average_expression
)[
  grepl(
    "^SoloTE-",
    rownames(
      average_expression
    )
  )
]

if (length(subfamily_features) == 0) {
  stop(
    "No subfamily-level SoloTE features were found in the merged object."
  )
}

subfamily_average <- average_expression[
  subfamily_features,
  ,
  drop = FALSE
]

present_cell_types <- intersect(
  cell_type_order,
  colnames(
    subfamily_average
  )
)

subfamily_average <- subfamily_average[
  ,
  present_cell_types,
  drop = FALSE
]

clean_subfamily_names <- sub(
  "^SoloTE-",
  "",
  rownames(
    subfamily_average
  )
)

if (anyDuplicated(clean_subfamily_names)) {
  stop(
    "Removing the SoloTE- prefix produced duplicated subfamily names."
  )
}

rownames(
  subfamily_average
) <- clean_subfamily_names

write.csv(
  as.data.frame(
    subfamily_average
  ),
  file.path(
    output_dir,
    average_expression_file
  ),
  quote = FALSE
)

subfamily_average_long <- as.data.frame(
  subfamily_average
) %>%
  rownames_to_column(
    "subfamily"
  ) %>%
  pivot_longer(
    cols = -subfamily,
    names_to = "cell_type",
    values_to = "average_log_normalized_expression"
  ) %>%
  mutate(
    cell_type = factor(
      cell_type,
      levels = cell_type_order
    )
  ) %>%
  arrange(
    subfamily,
    cell_type
  )

write.csv(
  subfamily_average_long,
  file.path(
    output_dir,
    average_expression_long_file
  ),
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------------------------
# 9. Save cell counts, parameters, and session information
# ------------------------------------------------------------------------------

cell_count_table <- as.data.frame(
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
  )

write.csv(
  cell_count_table,
  file.path(
    output_dir,
    "replicate_celltype_cell_counts.csv"
  ),
  row.names = FALSE
)

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
    "PCA and UMAP dimensions",
    "Clustering resolution",
    "Average-expression source"
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
    paste0(
      "1:",
      pca_dimensions
    ),
    clustering_resolution,
    "RNA data layer: log-normalized expression"
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
    "sessionInfo_SoloTE_subfamily_processing.txt"
  )
)

message("\nSubfamily-level SoloTE processing complete.")
message("Results were written to: ", output_dir)
