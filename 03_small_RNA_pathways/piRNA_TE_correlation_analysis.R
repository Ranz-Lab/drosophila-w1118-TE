```r
#!/usr/bin/env Rscript

# ==============================================================================
# Process SoloTE locus-level expression from Drosophila w1118 testis snRNA-seq
#
# This script:
#   1. Creates Seurat objects from SoloTE count matrices
#   2. Transfers Fly Cell Atlas testis annotations by matching cell barcodes
#   3. Filters cells without valid annotations
#   4. Merges and normalizes three biological samples
#   5. Aggregates locus-level TE expression across cell types
#   6. Classifies TE loci as germline-biased, soma-biased, or similar
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(matrixStats)
  library(Seurat)
  library(stringr)
})

set.seed(1234)


# ------------------------------------------------------------------------------
# 2. Input files and analysis parameters
# ------------------------------------------------------------------------------

fca_file <- "FCA_testis_10x.rds"

sample_table <- data.frame(
  sample_id = c(
    "Male_testis_sample1",
    "Male_testis_sample2",
    "Male_testis_sample3"
  ),
  fca_id = c(
    "FCA59",
    "FCA60",
    "FCA61"
  ),
  matrix_directory = c(
    "w1118S1_LBAM2_locustes_MATRIX",
    "w1118S2_LBAM_locustes_MATRIX",
    "w1118S3_LBAM2_locustes_MATRIX"
  ),
  output_rds = c(
    "Male_testis_sample1_filtered_FCA_soloTE_locus.rds",
    "Male_testis_sample2_filtered_FCA_soloTE_locus.rds",
    "Male_testis_sample3_filtered_FCA_soloTE_locus.rds"
  ),
  stringsAsFactors = FALSE
)

minimum_cells_per_feature <- 5
minimum_features_per_cell <- 200
bias_threshold <- 1


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
  "head cyst cell",
  "adult tracheal cell",
  "adult fat body",
  "hemocyte",
  "muscle cell",
  "male gonad associated epithelium",
  "pigment cell",
  "secretory cell of the male reproductive tract"
)


# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

strip_barcode_suffix <- function(cell_ids) {
  sub("-.*$", "", cell_ids)
}


transfer_fca_annotations <- function(
    solo_object,
    fca_object,
    fca_sample_id,
    sample_id
) {
  if (!"fca_id" %in% colnames(fca_object[[]])) {
    stop("The FCA object does not contain an 'fca_id' metadata column.")
  }

  if (!"S_annotation" %in% colnames(fca_object[[]])) {
    stop("The FCA object does not contain an 'S_annotation' metadata column.")
  }

  fca_cells <- Cells(fca_object)[fca_object$fca_id == fca_sample_id]

  if (length(fca_cells) == 0) {
    stop("No FCA cells were found for ", fca_sample_id, ".")
  }

  fca_sample <- subset(
    x = fca_object,
    cells = fca_cells
  )

  fca_metadata <- fca_sample[[]]

  annotation_lookup <- data.frame(
    barcode = strip_barcode_suffix(rownames(fca_metadata)),
    annotation = as.character(fca_metadata$S_annotation),
    stringsAsFactors = FALSE
  )

  valid_fca_annotations <- (
    !is.na(annotation_lookup$annotation) &
      nzchar(annotation_lookup$annotation) &
      toupper(annotation_lookup$annotation) != "NA"
  )

  annotation_lookup <- annotation_lookup[
    valid_fca_annotations,
    ,
    drop = FALSE
  ]

  # Confirm that duplicated stripped barcodes do not have conflicting
  # annotations within the FCA sample.
  barcode_conflicts <- annotation_lookup %>%
    group_by(barcode) %>%
    summarise(
      n_annotations = n_distinct(annotation),
      .groups = "drop"
    ) %>%
    filter(n_annotations > 1)

  if (nrow(barcode_conflicts) > 0) {
    stop(
      "Conflicting FCA annotations were found for stripped barcodes in ",
      fca_sample_id,
      "."
    )
  }

  annotation_lookup <- annotation_lookup %>%
    distinct(barcode, .keep_all = TRUE)

  solo_object$barcode <- strip_barcode_suffix(Cells(solo_object))
  solo_object$sample_id <- sample_id
  solo_object$fca_id <- fca_sample_id

  if (anyDuplicated(solo_object$barcode)) {
    stop(
      "Duplicated stripped barcodes were detected in ",
      sample_id,
      "."
    )
  }

  annotation_vector <- setNames(
    annotation_lookup$annotation,
    annotation_lookup$barcode
  )

  solo_object$annotation <- unname(
    annotation_vector[solo_object$barcode]
  )

  valid_target_annotations <- (
    !is.na(solo_object$annotation) &
      nzchar(solo_object$annotation) &
      toupper(solo_object$annotation) != "NA"
  )

  cells_to_keep <- Cells(solo_object)[valid_target_annotations]

  solo_object <- subset(
    x = solo_object,
    cells = cells_to_keep
  )

  if (ncol(solo_object) == 0) {
    stop("No annotated cells were retained for ", sample_id, ".")
  }

  Idents(solo_object) <- "annotation"

  solo_object
}


process_soloTE_sample <- function(
    sample_id,
    fca_sample_id,
    matrix_directory,
    output_rds,
    fca_object
) {
  message("\nProcessing ", sample_id, "...")

  if (!dir.exists(matrix_directory)) {
    stop("SoloTE matrix directory not found: ", matrix_directory)
  }

  solo_counts <- Read10X(
    data.dir = matrix_directory,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
  )

  solo_object <- CreateSeuratObject(
    counts = solo_counts,
    project = sample_id,
    min.cells = minimum_cells_per_feature,
    min.features = minimum_features_per_cell
  )

  original_cell_count <- ncol(solo_object)

  solo_object <- transfer_fca_annotations(
    solo_object = solo_object,
    fca_object = fca_object,
    fca_sample_id = fca_sample_id,
    sample_id = sample_id
  )

  saveRDS(
    object = solo_object,
    file = output_rds
  )

  message(
    sample_id,
    ": retained ",
    format(ncol(solo_object), big.mark = ","),
    " of ",
    format(original_cell_count, big.mark = ","),
    " cells with FCA annotations."
  )

  solo_object
}


parse_soloTE_features <- function(features) {
  feature_match <- str_match(
    features,
    paste0(
      "^SoloTE[-_]",
      "(X|2L|2R|3L|3R|4|Y)[-_]",
      "(\\d+)[-_]",
      "(\\d+)[-_]",
      "([^:]+):",
      "([^:]+):",
      "([^:-]+)[-_]",
      "(\\d+(?:\\.\\d+)?)[-_]",
      "([+-])$"
    )
  )

  data.frame(
    feature = features,
    chromosome = feature_match[, 2],
    start = as.numeric(feature_match[, 3]),
    end = as.numeric(feature_match[, 4]),
    family = feature_match[, 5],
    superfamily = feature_match[, 6],
    class = feature_match[, 7],
    score = as.numeric(feature_match[, 8]),
    strand = feature_match[, 9],
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 5. Load the Fly Cell Atlas testis object
# ------------------------------------------------------------------------------

if (!file.exists(fca_file)) {
  stop("FCA Seurat object not found: ", fca_file)
}

FCA_testis_10x <- readRDS(fca_file)

message(
  "Loaded FCA testis object containing ",
  format(ncol(FCA_testis_10x), big.mark = ","),
  " cells."
)


# ------------------------------------------------------------------------------
# 6. Process each SoloTE sample
# ------------------------------------------------------------------------------

sample_objects <- vector(
  mode = "list",
  length = nrow(sample_table)
)

names(sample_objects) <- sample_table$sample_id

for (i in seq_len(nrow(sample_table))) {
  sample_objects[[i]] <- process_soloTE_sample(
    sample_id = sample_table$sample_id[i],
    fca_sample_id = sample_table$fca_id[i],
    matrix_directory = sample_table$matrix_directory[i],
    output_rds = sample_table$output_rds[i],
    fca_object = FCA_testis_10x
  )
}


# ------------------------------------------------------------------------------
# 7. Merge and normalize the three samples
# ------------------------------------------------------------------------------

w1118_testis_soloTE <- merge(
  x = sample_objects[[1]],
  y = sample_objects[-1],
  add.cell.ids = sample_table$fca_id,
  project = "TE.Expression.Locus",
  merge.data = FALSE
)

# Seurat v5 stores merged count matrices in separate layers.
if (inherits(w1118_testis_soloTE[["RNA"]], "Assay5")) {
  w1118_testis_soloTE <- JoinLayers(w1118_testis_soloTE)
}

w1118_testis_soloTE <- NormalizeData(
  object = w1118_testis_soloTE,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

w1118_testis_soloTE <- FindVariableFeatures(
  object = w1118_testis_soloTE,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

w1118_testis_soloTE <- ScaleData(
  object = w1118_testis_soloTE,
  features = rownames(w1118_testis_soloTE),
  verbose = FALSE
)

Idents(w1118_testis_soloTE) <- "annotation"


# ------------------------------------------------------------------------------
# 8. Assign cells to germline and somatic categories
# ------------------------------------------------------------------------------

w1118_testis_soloTE$cell_category <- case_when(
  w1118_testis_soloTE$annotation %in% germline_cell_types ~ "Germline",
  w1118_testis_soloTE$annotation %in% somatic_cell_types ~ "Somatic",
  TRUE ~ NA_character_
)

unassigned_annotations <- sort(
  unique(
    w1118_testis_soloTE$annotation[
      is.na(w1118_testis_soloTE$cell_category)
    ]
  )
)

if (length(unassigned_annotations) > 0) {
  message(
    "\nThe following annotations were not assigned to germline or soma and ",
    "will be excluded from the comparison:\n",
    paste(unassigned_annotations, collapse = "\n")
  )
}

saveRDS(
  object = w1118_testis_soloTE,
  file = "FCA-testis-merged-normalized-scaled-soloTE-w1118.rds"
)


# ------------------------------------------------------------------------------
# 9. Save cell counts by sample and annotation
# ------------------------------------------------------------------------------

cell_counts <- w1118_testis_soloTE[[]] %>%
  count(
    sample_id,
    fca_id,
    annotation,
    cell_category,
    name = "n_cells"
  ) %>%
  arrange(
    sample_id,
    cell_category,
    annotation
  )

write.csv(
  cell_counts,
  file = "FCA-testis-soloTE-cell-counts-by-sample-and-annotation.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 10. Calculate average TE expression across annotated cell types
# ------------------------------------------------------------------------------

average_expression <- AverageExpression(
  object = w1118_testis_soloTE,
  assays = "RNA",
  group.by = "annotation",
  layer = "data",
  verbose = FALSE
)$RNA

soloTE_features <- rownames(average_expression)[
  grepl("^SoloTE[-_]", rownames(average_expression))
]

if (length(soloTE_features) == 0) {
  stop("No SoloTE locus-level features were identified in the merged object.")
}

average_expression <- average_expression[
  soloTE_features,
  ,
  drop = FALSE
]

# Standardize each TE locus across cell types.
average_expression_z <- t(
  scale(
    t(average_expression)
  )
)

average_expression_z[is.na(average_expression_z)] <- 0

saveRDS(
  object = average_expression,
  file = "FCA-testis-soloTE-average-expression-by-cell-type-w1118.rds"
)

saveRDS(
  object = average_expression_z,
  file = "FCA-testis-soloTE-cell-type-zscores-w1118.rds"
)


# ------------------------------------------------------------------------------
# 11. Compare locus-level expression in germline and somatic cell types
# ------------------------------------------------------------------------------

germline_present <- intersect(
  germline_cell_types,
  colnames(average_expression_z)
)

somatic_present <- intersect(
  somatic_cell_types,
  colnames(average_expression_z)
)

if (length(germline_present) == 0) {
  stop("No germline cell types were found in the averaged expression matrix.")
}

if (length(somatic_present) == 0) {
  stop("No somatic cell types were found in the averaged expression matrix.")
}

germline_expression <- rowMaxs(
  average_expression_z[
    ,
    germline_present,
    drop = FALSE
  ]
)

somatic_expression <- rowMaxs(
  average_expression_z[
    ,
    somatic_present,
    drop = FALSE
  ]
)

soloTE_metadata <- parse_soloTE_features(
  rownames(average_expression_z)
)

germline_soma_results <- soloTE_metadata %>%
  mutate(
    germline_max_z = germline_expression[feature],
    somatic_max_z = somatic_expression[feature],
    germline_minus_somatic = germline_max_z - somatic_max_z,
    expression_bias = case_when(
      germline_minus_somatic > bias_threshold ~ "Germline-biased",
      germline_minus_somatic < -bias_threshold ~ "Soma-biased",
      TRUE ~ "Similar"
    )
  ) %>%
  arrange(
    expression_bias,
    desc(abs(germline_minus_somatic))
  )

write.csv(
  germline_soma_results,
  file = "FCA-testis-soloTE-germline-vs-soma-w1118.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 12. Summarize expression bias by TE class
# ------------------------------------------------------------------------------

bias_summary <- germline_soma_results %>%
  mutate(
    class = ifelse(
      is.na(class) | !nzchar(class),
      "Unparsed",
      class
    )
  ) %>%
  count(
    class,
    expression_bias,
    name = "n_loci"
  ) %>%
  group_by(class) %>%
  mutate(
    total_loci = sum(n_loci),
    proportion = n_loci / total_loci
  ) %>%
  ungroup() %>%
  arrange(
    class,
    expression_bias
  )

write.csv(
  bias_summary,
  file = "FCA-testis-soloTE-germline-vs-soma-summary-by-class-w1118.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 13. Save session information
# ------------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  con = "sessionInfo_soloTE_processing.txt"
)

message("\nSoloTE processing complete.")
message(
  "Merged object: ",
  "FCA-testis-merged-normalized-scaled-soloTE-w1118.rds"
)
message(
  "Germline-versus-soma results: ",
  "FCA-testis-soloTE-germline-vs-soma-w1118.csv"
)
```
