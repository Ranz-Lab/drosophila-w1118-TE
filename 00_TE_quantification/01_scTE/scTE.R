#!/usr/bin/env Rscript

# ==============================================================================
# Process scTE results from Drosophila w1118 testis snRNA-seq
#
# This script:
#   1. Converts scTE .h5ad outputs to Seurat objects
#   2. Transfers Fly Cell Atlas testis annotations by cell barcode
#   3. Retains cells with valid annotations
#   4. Merges three biological replicates
#   5. Performs normalization, dimensionality reduction, and clustering
#   6. Saves the annotated Seurat object and UMAP
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(ggplot2)
})

set.seed(1234)


# ------------------------------------------------------------------------------
# 2. File paths and parameters
# ------------------------------------------------------------------------------

input_dir  <- "path/to/input"
output_dir <- "path/to/output"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Fly Cell Atlas Seurat object containing testis cell-type annotations
fca_object_file <- file.path(input_dir, "FCA_testis_10x.rds")

# scTE output files
sample_table <- data.frame(
  sample_id = c(
    "testis_sample1",
    "testis_sample2",
    "testis_sample3"
  ),
  h5ad_file = c(
    "Male_testis_sample1_S44Out5.h5ad",
    "Male_testis_sample2_S48Out.h5ad",
    "Male_testis_sample3_S60Out.h5ad"
  ),
  stringsAsFactors = FALSE
)

# Metadata column containing the FCA cell-type annotation
annotation_column <- "S_annotation"

# Optional: subset the FCA object by sample before matching barcodes.
#
# Leave this as NULL when cell barcodes are unique across the FCA object.
# If the same barcode occurs in multiple FCA samples with different annotations,
# specify the FCA sample metadata column and add the corresponding sample values
# to sample_table.
fca_sample_column <- NULL

sample_table$fca_sample_value <- NA_character_

# Analysis parameters
n_variable_features <- 2000
n_pcs               <- 30
cluster_resolution  <- 0.5


# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

# Remove the suffix beginning with the first hyphen from a 10x cell barcode.
strip_barcode_suffix <- function(cell_ids) {
  sub("-.*$", "", cell_ids)
}


# Construct a barcode-to-annotation lookup table from the FCA object.
create_annotation_lookup <- function(
    fca_object,
    annotation_column,
    sample_column = NULL,
    sample_value = NULL
) {
  fca_metadata <- fca_object[[]]

  if (!annotation_column %in% colnames(fca_metadata)) {
    stop(
      "Annotation column '",
      annotation_column,
      "' was not found in the FCA metadata."
    )
  }

  # Optionally restrict the FCA metadata to one biological sample.
  if (!is.null(sample_column)) {
    if (!sample_column %in% colnames(fca_metadata)) {
      stop(
        "FCA sample column '",
        sample_column,
        "' was not found in the FCA metadata."
      )
    }

    if (is.null(sample_value) || is.na(sample_value)) {
      stop(
        "An FCA sample value must be supplied when fca_sample_column is used."
      )
    }

    fca_metadata <- fca_metadata[
      fca_metadata[[sample_column]] == sample_value,
      ,
      drop = FALSE
    ]

    if (nrow(fca_metadata) == 0) {
      stop(
        "No FCA cells were found for ",
        sample_column,
        " = ",
        sample_value,
        "."
      )
    }
  }

  annotation_lookup <- data.frame(
    barcode = strip_barcode_suffix(rownames(fca_metadata)),
    annotation = as.character(fca_metadata[[annotation_column]]),
    stringsAsFactors = FALSE
  )

  # Remove missing or empty annotations.
  valid_annotation <- (
    !is.na(annotation_lookup$annotation) &
      nzchar(annotation_lookup$annotation) &
      toupper(annotation_lookup$annotation) != "NA"
  )

  annotation_lookup <- annotation_lookup[
    valid_annotation,
    ,
    drop = FALSE
  ]

  # Check whether duplicated barcodes have conflicting annotations.
  annotations_by_barcode <- split(
    annotation_lookup$annotation,
    annotation_lookup$barcode
  )

  conflicting_barcodes <- names(
    Filter(
      function(x) length(unique(x)) > 1,
      annotations_by_barcode
    )
  )

  if (length(conflicting_barcodes) > 0) {
    example_barcodes <- paste(
      head(conflicting_barcodes, 10),
      collapse = ", "
    )

    stop(
      "FCA barcodes were associated with conflicting annotations. ",
      "Subset the FCA object by sample using fca_sample_column. ",
      "Example conflicting barcodes: ",
      example_barcodes
    )
  }

  # Duplicated barcodes with the same annotation can be safely collapsed.
  annotation_lookup <- annotation_lookup[
    !duplicated(annotation_lookup$barcode),
    ,
    drop = FALSE
  ]

  annotation_lookup
}


# Convert one scTE .h5ad file and return an annotated Seurat object.
process_scTE_sample <- function(
    sample_id,
    h5ad_file,
    fca_object,
    annotation_column,
    input_dir,
    fca_sample_column = NULL,
    fca_sample_value = NULL
) {
  message("\nProcessing ", sample_id, "...")

  h5ad_path <- file.path(input_dir, h5ad_file)

  if (!file.exists(h5ad_path)) {
    stop("Input file not found: ", h5ad_path)
  }

  h5seurat_path <- sub(
    "\\.h5ad$",
    ".h5seurat",
    h5ad_path,
    ignore.case = TRUE
  )

  # Convert the scTE output to H5Seurat format.
  Convert(
    source = h5ad_path,
    dest = "h5seurat",
    overwrite = TRUE
  )

  if (!file.exists(h5seurat_path)) {
    stop("H5Seurat conversion failed for: ", h5ad_path)
  }

  scTE_object <- LoadH5Seurat(h5seurat_path)

  original_cell_count <- ncol(scTE_object)

  # Store sample metadata.
  scTE_object$sample_id <- sample_id
  scTE_object$strain    <- "w1118"

  # Standardize barcodes for matching.
  scTE_object$CellName <- strip_barcode_suffix(Cells(scTE_object))

  if (anyDuplicated(scTE_object$CellName)) {
    stop(
      "Duplicated stripped cell barcodes were detected within ",
      sample_id,
      "."
    )
  }

  annotation_lookup <- create_annotation_lookup(
    fca_object = fca_object,
    annotation_column = annotation_column,
    sample_column = fca_sample_column,
    sample_value = fca_sample_value
  )

  annotation_vector <- setNames(
    annotation_lookup$annotation,
    annotation_lookup$barcode
  )

  # Transfer annotations by barcode rather than relying on row order.
  scTE_object$annotation <- unname(
    annotation_vector[scTE_object$CellName]
  )

  cells_to_keep <- (
    !is.na(scTE_object$annotation) &
      nzchar(scTE_object$annotation) &
      toupper(scTE_object$annotation) != "NA"
  )

  scTE_object <- subset(
    scTE_object,
    cells = Cells(scTE_object)[cells_to_keep]
  )

  if (ncol(scTE_object) == 0) {
    stop("No annotated cells were retained for ", sample_id, ".")
  }

  Idents(scTE_object) <- "annotation"

  message(
    sample_id,
    ": retained ",
    format(ncol(scTE_object), big.mark = ","),
    " of ",
    format(original_cell_count, big.mark = ","),
    " cells with FCA annotations."
  )

  scTE_object
}


# ------------------------------------------------------------------------------
# 4. Load the Fly Cell Atlas object
# ------------------------------------------------------------------------------

if (!file.exists(fca_object_file)) {
  stop("FCA Seurat object not found: ", fca_object_file)
}

FCA_testis_10x <- readRDS(fca_object_file)

message(
  "Loaded FCA testis object containing ",
  format(ncol(FCA_testis_10x), big.mark = ","),
  " cells."
)


# ------------------------------------------------------------------------------
# 5. Convert, annotate, and filter each scTE sample
# ------------------------------------------------------------------------------

sample_objects <- vector(
  mode = "list",
  length = nrow(sample_table)
)

names(sample_objects) <- sample_table$sample_id

for (i in seq_len(nrow(sample_table))) {
  sample_objects[[i]] <- process_scTE_sample(
    sample_id = sample_table$sample_id[i],
    h5ad_file = sample_table$h5ad_file[i],
    fca_object = FCA_testis_10x,
    annotation_column = annotation_column,
    input_dir = input_dir,
    fca_sample_column = fca_sample_column,
    fca_sample_value = sample_table$fca_sample_value[i]
  )
}


# ------------------------------------------------------------------------------
# 6. Merge the three biological replicates
# ------------------------------------------------------------------------------

testis_scTE <- merge(
  x = sample_objects[[1]],
  y = sample_objects[-1],
  add.cell.ids = names(sample_objects),
  project = "w1118_testis_scTE",
  merge.data = FALSE
)

Idents(testis_scTE) <- "annotation"

message(
  "\nMerged object contains ",
  format(ncol(testis_scTE), big.mark = ","),
  " cells and ",
  format(nrow(testis_scTE), big.mark = ","),
  " features."
)


# ------------------------------------------------------------------------------
# 7. Normalize and perform dimensionality reduction
# ------------------------------------------------------------------------------

testis_scTE <- NormalizeData(
  object = testis_scTE,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

testis_scTE <- FindVariableFeatures(
  object = testis_scTE,
  selection.method = "vst",
  nfeatures = n_variable_features,
  verbose = FALSE
)

testis_scTE <- ScaleData(
  object = testis_scTE,
  features = VariableFeatures(testis_scTE),
  verbose = FALSE
)

testis_scTE <- RunPCA(
  object = testis_scTE,
  features = VariableFeatures(testis_scTE),
  npcs = n_pcs,
  seed.use = 1234,
  verbose = FALSE
)

dimensions_to_use <- seq_len(
  min(n_pcs, ncol(Embeddings(testis_scTE, reduction = "pca")))
)

testis_scTE <- FindNeighbors(
  object = testis_scTE,
  dims = dimensions_to_use,
  verbose = FALSE
)

testis_scTE <- FindClusters(
  object = testis_scTE,
  resolution = cluster_resolution,
  random.seed = 1234,
  verbose = FALSE
)

testis_scTE <- RunUMAP(
  object = testis_scTE,
  dims = dimensions_to_use,
  seed.use = 1234,
  verbose = FALSE
)


# ------------------------------------------------------------------------------
# 8. Generate and save the annotated UMAP
# ------------------------------------------------------------------------------

umap_plot <- DimPlot(
  object = testis_scTE,
  reduction = "umap",
  group.by = "annotation",
  label = TRUE,
  label.size = 3,
  repel = TRUE
) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      override.aes = list(size = 3)
    )
  ) +
  theme_classic()

print(umap_plot)

ggsave(
  filename = file.path(
    output_dir,
    "FCA-testis-merged-scTE-annotated-UMAP-w1118.pdf"
  ),
  plot = umap_plot,
  width = 10,
  height = 7
)


# ------------------------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------------------------

output_rds <- file.path(
  output_dir,
  "FCA-testis-merged-scTE-annotated-normalized-w1118.rds"
)

saveRDS(
  object = testis_scTE,
  file = output_rds
)

annotation_counts <- as.data.frame(
  table(
    sample_id = testis_scTE$sample_id,
    annotation = testis_scTE$annotation
  )
)

annotation_counts <- annotation_counts[
  annotation_counts$Freq > 0,
  ,
  drop = FALSE
]

write.csv(
  annotation_counts,
  file = file.path(
    output_dir,
    "FCA-testis-scTE-cell-counts-by-sample-and-annotation.csv"
  ),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

message("\nSaved annotated Seurat object to:")
message(output_rds)
message("\nProcessing complete.")
