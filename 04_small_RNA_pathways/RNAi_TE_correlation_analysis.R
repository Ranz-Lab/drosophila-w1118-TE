#!/usr/bin/env Rscript

# ==============================================================================
# RNAi pathway gene and TE expression analysis
#
# This script generates the data underlying the RNAi pathway analysis:
#
#   A. RNAi pathway gene expression across annotated testis cell types
#   B. Within-cell-type Spearman correlations between RNAi gene expression and
#      TE-class expression loads across individual nuclei
#   C. Between-cell-type Spearman correlations between RNAi gene expression and
#      TE-class expression loads within germline and somatic compartments
#   D. Cell-type medians and correlations for loqs versus DNA, LINE, and LTR
#      expression loads
#
#
# Required input files:
#   FCA-testis-merged-umap-final-soloTE-subfamily-w1118.rds
#   TE-class_mapping-soloTE.xlsx
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Packages and parameters
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(readxl)
  library(Seurat)
  library(stringr)
  library(tibble)
  library(tidyr)
})

seurat_file <- "FCA-testis-merged-umap-final-soloTE-subfamily-w1118.rds"
te_mapping_file <- "TE-class_mapping-soloTE.xlsx"
output_dir <- "RNAi_TE_results"

rna_i_genes <- c("Dcr-2", "AGO2", "r2d2", "loqs", "Hen1")
te_classes <- c("DNA", "LINE", "LTR")

minimum_nonzero_fraction_gene <- 0.01
minimum_nonzero_fraction_te <- 0.01
minimum_observations <- 3
loqs_lower_quantile_filter <- 0.10

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 2. Cell-type definitions
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

known_cell_types <- c(germline_cell_types, somatic_cell_types)

assign_compartment <- function(cell_type) {
  case_when(
    cell_type %in% germline_cell_types ~ "Germline",
    cell_type %in% somatic_cell_types ~ "Somatic",
    TRUE ~ NA_character_
  )
}

assign_stage <- function(cell_type) {
  case_when(
    cell_type %in% mitotic_cell_types ~ "Mitotic",
    cell_type %in% meiotic_cell_types ~ "Meiotic",
    cell_type %in% postmeiotic_cell_types ~ "Post-meiotic",
    cell_type %in% somatic_cell_types ~ "Somatic",
    TRUE ~ NA_character_
  )
}


# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

get_normalized_matrix <- function(seurat_object) {
  if (inherits(seurat_object[["RNA"]], "Assay5")) {
    return(GetAssayData(seurat_object, assay = "RNA", layer = "data"))
  }

  GetAssayData(seurat_object, assay = "RNA", slot = "data")
}

average_expression_by_celltype <- function(seurat_object, annotation_column) {
  tryCatch(
    AverageExpression(
      seurat_object,
      assays = "RNA",
      group.by = annotation_column,
      layer = "data",
      verbose = FALSE
    )$RNA,
    error = function(e) {
      AverageExpression(
        seurat_object,
        assays = "RNA",
        group.by = annotation_column,
        slot = "data",
        verbose = FALSE
      )$RNA
    }
  )
}

match_features_case_insensitive <- function(requested, available) {
  index <- match(tolower(requested), tolower(available))
  matched <- available[index[!is.na(index)]]
  names(matched) <- requested[!is.na(index)]
  matched
}

spearman_test <- function(x, y, minimum_n = 3) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (
    length(x) < minimum_n ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(
      data.frame(
        n = length(x),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  result <- suppressWarnings(
    cor.test(x, y, method = "spearman", exact = FALSE)
  )

  data.frame(
    n = length(x),
    rho = unname(result$estimate),
    p_value = result$p.value
  )
}


# ------------------------------------------------------------------------------
# 4. Load and validate the Seurat object
# ------------------------------------------------------------------------------

if (!file.exists(seurat_file)) {
  stop("Seurat object not found: ", seurat_file)
}

Control.samples <- readRDS(seurat_file)

annotation_candidates <- c("newannotation", "annotation", "S_annotation")
annotation_column <- annotation_candidates[
  annotation_candidates %in% colnames(Control.samples[[]])
][1]

if (is.na(annotation_column) || length(annotation_column) == 0) {
  stop(
    "No supported annotation column was found. Expected one of: ",
    paste(annotation_candidates, collapse = ", ")
  )
}

annotations <- as.character(Control.samples[[]][[annotation_column]])
names(annotations) <- Cells(Control.samples)

valid_cells <- (
  !is.na(annotations) &
  nzchar(annotations) &
  annotations != "unannotated"
)

Control.samples <- subset(
  Control.samples,
  cells = names(annotations)[valid_cells]
)

annotations <- as.character(Control.samples[[]][[annotation_column]])
names(annotations) <- Cells(Control.samples)

unknown_cell_types <- setdiff(unique(annotations), known_cell_types)

if (length(unknown_cell_types) > 0) {
  stop(
    "The following cell types are not assigned to a compartment:\n",
    paste(sort(unknown_cell_types), collapse = "\n")
  )
}

normalized_matrix <- get_normalized_matrix(Control.samples)

if (!identical(colnames(normalized_matrix), Cells(Control.samples))) {
  stop("Expression-matrix columns are not aligned with Seurat cell names.")
}


# ------------------------------------------------------------------------------
# 5. Match RNAi genes and map TE families to TE classes
# ------------------------------------------------------------------------------

matched_genes <- match_features_case_insensitive(
  rna_i_genes,
  rownames(normalized_matrix)
)

missing_genes <- setdiff(rna_i_genes, names(matched_genes))

if (length(missing_genes) > 0) {
  stop(
    "The following RNAi genes were not found: ",
    paste(missing_genes, collapse = ", ")
  )
}

if (!file.exists(te_mapping_file)) {
  stop("TE mapping file not found: ", te_mapping_file)
}

mapping_raw <- read_excel(te_mapping_file)
mapping_names <- colnames(mapping_raw)

superfamily_column <- mapping_names[
  grepl("^superfamily$", mapping_names, ignore.case = TRUE)
][1]

class_column <- mapping_names[
  grepl("^(TE_)?class$", mapping_names, ignore.case = TRUE)
][1]

if (is.na(superfamily_column) || is.na(class_column)) {
  stop(
    "The TE mapping file must contain 'superfamily' and ",
    "'class' or 'TE_class' columns."
  )
}

te_mapping <- mapping_raw %>%
  transmute(
    TE_family = trimws(as.character(.data[[superfamily_column]])),
    TE_class = toupper(trimws(as.character(.data[[class_column]]))),
    mapping_key = tolower(TE_family)
  ) %>%
  filter(
    !is.na(TE_family),
    nzchar(TE_family),
    TE_class %in% te_classes
  ) %>%
  distinct(mapping_key, TE_class, .keep_all = TRUE)

mapping_conflicts <- te_mapping %>%
  count(mapping_key) %>%
  filter(n > 1)

if (nrow(mapping_conflicts) > 0) {
  stop("Some TE families map to more than one TE class.")
}

te_features <- rownames(normalized_matrix)[
  str_detect(rownames(normalized_matrix), "^(SoloTE[-_]|TE[-_])")
]

te_metadata <- tibble(
  TE_feature = te_features,
  TE_family = str_remove(te_features, "^(SoloTE[-_]|TE[-_])"),
  mapping_key = tolower(trimws(TE_family))
) %>%
  left_join(
    te_mapping %>% select(mapping_key, TE_class),
    by = "mapping_key"
  ) %>%
  filter(!is.na(TE_class)) %>%
  distinct(TE_feature, .keep_all = TRUE) %>%
  select(TE_feature, TE_family, TE_class)

if (nrow(te_metadata) == 0) {
  stop("No TE features mapped to DNA, LINE, or LTR.")
}

write.csv(
  te_metadata,
  file.path(output_dir, "TE_family_class_mapping_used.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    TE_feature = setdiff(te_features, te_metadata$TE_feature)
  ),
  file.path(output_dir, "unmapped_TE_features.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 6. Build per-cell RNAi expression and TE-class loads
# ------------------------------------------------------------------------------

cell_table <- tibble(
  cell = Cells(Control.samples),
  cell_type = annotations[Cells(Control.samples)]
) %>%
  mutate(
    compartment = assign_compartment(cell_type),
    stage = assign_stage(cell_type)
  )

gene_expression_per_cell <- as.data.frame(
  t(
    as.matrix(
      normalized_matrix[
        unname(matched_genes),
        ,
        drop = FALSE
      ]
    )
  )
)

colnames(gene_expression_per_cell) <- names(matched_genes)

gene_expression_per_cell <- gene_expression_per_cell %>%
  rownames_to_column("cell")

linear_expression <- expm1(normalized_matrix)

class_to_features <- split(
  te_metadata$TE_feature,
  te_metadata$TE_class
)

te_loads_per_cell <- lapply(te_classes, function(te_class) {
  features <- intersect(
    class_to_features[[te_class]],
    rownames(linear_expression)
  )

  if (length(features) == 0) {
    stop("No mapped TE features were found for class: ", te_class)
  }

  log1p(
    Matrix::colSums(
      linear_expression[
        features,
        ,
        drop = FALSE
      ]
    )
  )
})

names(te_loads_per_cell) <- te_classes

te_loads_per_cell <- as.data.frame(te_loads_per_cell) %>%
  mutate(cell = colnames(normalized_matrix))

per_cell_data <- cell_table %>%
  inner_join(gene_expression_per_cell, by = "cell") %>%
  inner_join(te_loads_per_cell, by = "cell")


# ------------------------------------------------------------------------------
# 7. RNAi pathway expression across cell types
# ------------------------------------------------------------------------------

rna_i_expression_by_celltype <- per_cell_data %>%
  pivot_longer(
    cols = all_of(rna_i_genes),
    names_to = "RNAi_gene",
    values_to = "expression"
  ) %>%
  group_by(
    cell_type,
    compartment,
    stage,
    RNAi_gene
  ) %>%
  summarise(
    n_cells = n(),
    average_expression = mean(expression),
    median_expression = median(expression),
    fraction_expressing = mean(expression > 0),
    .groups = "drop"
  ) %>%
  arrange(RNAi_gene, compartment, cell_type)

write.csv(
  rna_i_expression_by_celltype,
  file.path(output_dir, "RNAi_gene_expression_by_celltype.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8. Within-cell-type RNAi gene versus TE-class correlations
# ------------------------------------------------------------------------------

within_results <- list()
result_index <- 1

for (cell_type_name in sort(unique(per_cell_data$cell_type))) {
  cell_type_data <- per_cell_data %>%
    filter(cell_type == cell_type_name)

  for (te_class in te_classes) {
    te_values <- cell_type_data[[te_class]]
    te_nonzero_fraction <- mean(te_values > 0)

    for (gene_name in rna_i_genes) {
      gene_values <- cell_type_data[[gene_name]]
      gene_nonzero_fraction <- mean(gene_values > 0)

      if (
        gene_nonzero_fraction < minimum_nonzero_fraction_gene ||
        te_nonzero_fraction < minimum_nonzero_fraction_te
      ) {
        test_result <- data.frame(
          n = nrow(cell_type_data),
          rho = NA_real_,
          p_value = NA_real_
        )
      } else {
        test_result <- spearman_test(
          gene_values,
          te_values,
          minimum_n = minimum_observations
        )
      }

      within_results[[result_index]] <- data.frame(
        cell_type = cell_type_name,
        compartment = unique(cell_type_data$compartment),
        stage = unique(cell_type_data$stage),
        RNAi_gene = gene_name,
        TE_class = te_class,
        n_cells = test_result$n,
        fraction_expressing_gene = gene_nonzero_fraction,
        fraction_nonzero_TE_load = te_nonzero_fraction,
        rho = test_result$rho,
        p_value = test_result$p_value,
        stringsAsFactors = FALSE
      )

      result_index <- result_index + 1
    }
  }
}

within_correlations <- bind_rows(within_results) %>%
  group_by(cell_type, TE_class) %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH"),
    significant_FDR_0.05 = (
      !is.na(p_adjusted_BH) &
      p_adjusted_BH < 0.05
    )
  ) %>%
  ungroup() %>%
  arrange(compartment, stage, cell_type, TE_class, RNAi_gene)

write.csv(
  within_correlations,
  file.path(
    output_dir,
    "within_celltype_RNAi_TEclass_correlations.csv"
  ),
  row.names = FALSE
)

within_summary_by_compartment <- within_correlations %>%
  group_by(compartment, RNAi_gene, TE_class) %>%
  summarise(
    mean_rho = mean(rho, na.rm = TRUE),
    median_rho = median(rho, na.rm = TRUE),
    n_cell_types = sum(!is.na(rho)),
    fraction_negative = mean(rho < 0, na.rm = TRUE),
    fraction_significant_FDR_0.05 = mean(
      significant_FDR_0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(compartment, RNAi_gene, TE_class)

write.csv(
  within_summary_by_compartment,
  file.path(
    output_dir,
    "within_celltype_RNAi_TEclass_summary_by_compartment.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 9. Between-cell-type correlations within germline and soma
# ------------------------------------------------------------------------------

average_expression <- average_expression_by_celltype(
  Control.samples,
  annotation_column
)

gene_average <- average_expression[
  unname(matched_genes),
  ,
  drop = FALSE
]

rownames(gene_average) <- names(matched_genes)

te_average <- average_expression[
  te_metadata$TE_feature,
  ,
  drop = FALSE
]

te_class_load_by_celltype <- sapply(te_classes, function(te_class) {
  features <- te_metadata$TE_feature[
    te_metadata$TE_class == te_class
  ]

  features <- intersect(features, rownames(te_average))

  if (length(features) == 0) {
    return(rep(NA_real_, ncol(te_average)))
  }

  colSums(
    as.matrix(
      te_average[
        features,
        ,
        drop = FALSE
      ]
    ),
    na.rm = TRUE
  )
})

rownames(te_class_load_by_celltype) <- colnames(te_average)

between_results <- list()
result_index <- 1

for (compartment_name in c("Germline", "Somatic")) {
  requested_cell_types <- if (compartment_name == "Germline") {
    germline_cell_types
  } else {
    somatic_cell_types
  }

  cell_types_present <- Reduce(
    intersect,
    list(
      requested_cell_types,
      colnames(gene_average),
      rownames(te_class_load_by_celltype)
    )
  )

  if (length(cell_types_present) < minimum_observations) {
    stop(
      compartment_name,
      " has fewer than ",
      minimum_observations,
      " available cell types."
    )
  }

  for (te_class in te_classes) {
    te_values <- te_class_load_by_celltype[
      cell_types_present,
      te_class
    ]

    for (gene_name in rna_i_genes) {
      gene_values <- as.numeric(
        gene_average[
          gene_name,
          cell_types_present,
          drop = TRUE
        ]
      )

      test_result <- spearman_test(
        gene_values,
        te_values,
        minimum_n = minimum_observations
      )

      between_results[[result_index]] <- data.frame(
        compartment = compartment_name,
        RNAi_gene = gene_name,
        TE_class = te_class,
        n_cell_types = test_result$n,
        rho = test_result$rho,
        p_value = test_result$p_value,
        stringsAsFactors = FALSE
      )

      result_index <- result_index + 1
    }
  }
}

between_correlations <- bind_rows(between_results) %>%
  group_by(compartment, TE_class) %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH"),
    significant_FDR_0.05 = (
      !is.na(p_adjusted_BH) &
      p_adjusted_BH < 0.05
    )
  ) %>%
  ungroup() %>%
  arrange(compartment, TE_class, RNAi_gene)

write.csv(
  between_correlations,
  file.path(
    output_dir,
    "between_celltype_RNAi_TEclass_correlations.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 10. loqs versus TE-class load across cell types
# ------------------------------------------------------------------------------

loqs_long <- per_cell_data %>%
  select(
    cell,
    cell_type,
    compartment,
    stage,
    loqs,
    all_of(te_classes)
  ) %>%
  pivot_longer(
    cols = all_of(te_classes),
    names_to = "TE_class",
    values_to = "TE_load"
  ) %>%
  filter(
    is.finite(loqs),
    is.finite(TE_load),
    loqs > 0,
    TE_load > 0
  )

loqs_cutoffs <- loqs_long %>%
  group_by(TE_class) %>%
  summarise(
    loqs_threshold = quantile(
      loqs,
      probs = loqs_lower_quantile_filter,
      na.rm = TRUE
    ),
    TE_load_threshold = quantile(
      TE_load,
      probs = loqs_lower_quantile_filter,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

loqs_filtered <- loqs_long %>%
  inner_join(loqs_cutoffs, by = "TE_class") %>%
  filter(
    loqs >= loqs_threshold,
    TE_load >= TE_load_threshold
  )

loqs_celltype_summary <- loqs_filtered %>%
  group_by(
    TE_class,
    compartment,
    stage,
    cell_type
  ) %>%
  summarise(
    n_cells = n(),
    median_loqs_expression = median(loqs),
    median_TE_class_load = median(TE_load),
    .groups = "drop"
  ) %>%
  arrange(TE_class, compartment, stage, cell_type)

write.csv(
  loqs_celltype_summary,
  file.path(
    output_dir,
    "loqs_TEclass_celltype_medians.csv"
  ),
  row.names = FALSE
)

loqs_correlation_results <- list()

for (i in seq_along(te_classes)) {
  te_class <- te_classes[i]

  class_data <- loqs_celltype_summary %>%
    filter(TE_class == te_class)

  test_result <- spearman_test(
    class_data$median_loqs_expression,
    class_data$median_TE_class_load,
    minimum_n = minimum_observations
  )

  loqs_correlation_results[[i]] <- data.frame(
    TE_class = te_class,
    n_cell_types = test_result$n,
    rho = test_result$rho,
    p_value = test_result$p_value,
    stringsAsFactors = FALSE
  )
}

loqs_correlations <- bind_rows(loqs_correlation_results) %>%
  mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH"),
    significant_FDR_0.05 = (
      !is.na(p_adjusted_BH) &
      p_adjusted_BH < 0.05
    )
  ) %>%
  arrange(TE_class)

write.csv(
  loqs_correlations,
  file.path(
    output_dir,
    "loqs_TEclass_correlations_across_celltypes.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 11. Save parameters and session information
# ------------------------------------------------------------------------------

analysis_parameters <- data.frame(
  parameter = c(
    "Seurat object",
    "TE mapping file",
    "RNAi genes",
    "TE classes",
    "Minimum nonzero fraction for RNAi genes",
    "Minimum nonzero fraction for TE-class loads",
    "Minimum observations per correlation",
    "loqs lower-quantile filter",
    "Within-cell-type FDR scope",
    "Between-cell-type FDR scope"
  ),
  value = c(
    seurat_file,
    te_mapping_file,
    paste(rna_i_genes, collapse = ", "),
    paste(te_classes, collapse = ", "),
    minimum_nonzero_fraction_gene,
    minimum_nonzero_fraction_te,
    minimum_observations,
    loqs_lower_quantile_filter,
    "BH across RNAi genes within each cell type and TE class",
    "BH across RNAi genes within each compartment and TE class"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  analysis_parameters,
  file.path(output_dir, "analysis_parameters.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_RNAi_TE_analysis.txt")
)

message("\nRNAi-TE analysis complete.")
message("Results were written to: ", output_dir)
