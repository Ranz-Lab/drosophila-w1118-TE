```r
#!/usr/bin/env Rscript

# ==============================================================================
# Mixed-effects analysis of euchromatic and heterochromatic TE expression
#
# This script compares mean scaled expression from euchromatic (EU) and
# heterochromatic (HET) TE loci in Drosophila w1118 testis cell types.
#
# Separate models are fitted for:
#   1. DNA transposons
#   2. Retrotransposons
#
# The Retro input is assumed to combine LINE and LTR loci, as in the original
# analysis.
#
# For each TE group, the script fits:
#
#   Cell-type model:
#     mean_scaled_expr ~ partition * newannotation +
#       (1 | replicate) +
#       (1 | replicate:newannotation)
#
#   Compartment model:
#     mean_scaled_expr ~ partition * compartment +
#       (1 | replicate) +
#       (1 | newannotation)
#
# The script exports Type III ANOVA tables, estimated marginal means,
# EU-versus-HET contrasts with global Benjamini-Hochberg correction, model
# objects, model summaries, diagnostics, and session information.
#
# Required input files:
#   INPUT_LME_EuHet_DNA_locusTE_replicate_celltype.csv
#   INPUT_LME_EuHet_REtro_locusTE_replicate_celltype.csv
#
# Required input columns:
#   replicate
#   newannotation
#   partition
#   mean_scaled_expr
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(lme4)
  library(lmerTest)
})

# Sum-to-zero contrasts are required for interpretable Type III tests.
options(contrasts = c("contr.sum", "contr.poly"))


# ------------------------------------------------------------------------------
# 2. Input files and output directory
# ------------------------------------------------------------------------------

input_table <- data.frame(
  TE_group = c(
    "DNA",
    "Retro"
  ),
  input_file = c(
    "INPUT_LME_EuHet_DNA_locusTE_replicate_celltype.csv",
    "INPUT_LME_EuHet_REtro_locusTE_replicate_celltype.csv"
  ),
  stringsAsFactors = FALSE
)

output_dir <- "EuHet_mixed_model_results"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------------------------
# 3. Define testis cell types and compartments
# ------------------------------------------------------------------------------

Cell_Types_Order <- c(
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
  germline_cell_types,
  somatic_cell_types
)


# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

prepare_input_data <- function(input_file, TE_group) {

  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  df <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c(
    "replicate",
    "newannotation",
    "partition",
    "mean_scaled_expr"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )

  if (length(missing_columns) > 0) {
    stop(
      TE_group,
      " input is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (nrow(df) == 0) {
    stop(TE_group, " input contains no rows.")
  }

  df$mean_scaled_expr <- suppressWarnings(
    as.numeric(df$mean_scaled_expr)
  )

  if (any(!is.finite(df$mean_scaled_expr))) {
    stop(
      TE_group,
      " mean_scaled_expr contains missing, infinite, or non-numeric values."
    )
  }

  unknown_cell_types <- sort(
    setdiff(
      unique(df$newannotation),
      known_cell_types
    )
  )

  if (length(unknown_cell_types) > 0) {
    stop(
      TE_group,
      " input contains cell types that are not assigned to a compartment:\n",
      paste(unknown_cell_types, collapse = "\n")
    )
  }

  unknown_partitions <- sort(
    setdiff(
      unique(df$partition),
      c("EU", "HET")
    )
  )

  if (length(unknown_partitions) > 0) {
    stop(
      TE_group,
      " input contains unexpected partition values: ",
      paste(unknown_partitions, collapse = ", ")
    )
  }

  duplicate_rows <- df %>%
    count(
      replicate,
      newannotation,
      partition,
      name = "n_rows"
    ) %>%
    filter(n_rows > 1)

  if (nrow(duplicate_rows) > 0) {
    stop(
      TE_group,
      " input contains multiple rows for the same ",
      "replicate-cell type-partition combination."
    )
  }

  partition_counts <- df %>%
    distinct(
      replicate,
      newannotation,
      partition
    ) %>%
    count(
      replicate,
      newannotation,
      name = "n_partitions"
    ) %>%
    filter(n_partitions != 2)

  if (nrow(partition_counts) > 0) {
    stop(
      TE_group,
      " input contains replicate-cell type combinations without both ",
      "EU and HET observations."
    )
  }

  df <- df %>%
    mutate(
      compartment = case_when(
        newannotation %in% somatic_cell_types ~ "Somatic",
        newannotation %in% germline_cell_types ~ "Germline",
        TRUE ~ NA_character_
      ),
      replicate = factor(replicate),
      partition = factor(
        partition,
        levels = c("EU", "HET")
      ),
      compartment = factor(
        compartment,
        levels = c("Somatic", "Germline")
      ),
      newannotation = factor(
        newannotation,
        levels = Cell_Types_Order
      )
    )

  if (nlevels(df$replicate) < 2) {
    stop(
      TE_group,
      " analysis requires at least two biological replicates."
    )
  }

  df
}


format_anova_table <- function(anova_result) {

  anova_table <- as.data.frame(anova_result)

  anova_table$term <- rownames(anova_table)

  anova_table <- anova_table[
    ,
    c(
      "term",
      setdiff(
        colnames(anova_table),
        "term"
      )
    ),
    drop = FALSE
  ]

  rownames(anova_table) <- NULL

  anova_table
}


add_partition_interpretation <- function(contrast_table) {

  contrast_table %>%
    mutate(
      padj_fdr_global = p.adjust(
        p.value,
        method = "BH"
      ),
      direction = case_when(
        contrast == "EU - HET" & estimate > 0 ~ "EU > HET",
        contrast == "EU - HET" & estimate < 0 ~ "HET > EU",
        TRUE ~ "No difference"
      ),
      significant_fdr_0.05 = padj_fdr_global < 0.05
    )
}


write_model_diagnostics <- function(
    model,
    output_file,
    model_name
) {

  convergence_messages <- model@optinfo$conv$lme4$messages

  if (is.null(convergence_messages)) {
    convergence_text <- "None"
  } else {
    convergence_text <- paste(
      convergence_messages,
      collapse = "; "
    )
  }

  singular_fit <- isSingular(
    model,
    tol = 1e-4
  )

  random_effects <- as.data.frame(
    VarCorr(model)
  )

  diagnostic_lines <- c(
    paste0("Model: ", model_name),
    "",
    paste0("Number of observations: ", nobs(model)),
    paste0("AIC: ", AIC(model)),
    paste0("BIC: ", BIC(model)),
    paste0("Log likelihood: ", logLik(model)),
    paste0("Singular fit: ", singular_fit),
    paste0("Convergence messages: ", convergence_text),
    "",
    "Random-effect variance estimates:",
    capture.output(print(random_effects)),
    "",
    "Model summary:",
    capture.output(summary(model))
  )

  writeLines(
    diagnostic_lines,
    output_file
  )
}


fit_celltype_model <- function(
    df,
    TE_group,
    output_dir
) {

  message(
    "Fitting cell-type model for ",
    TE_group,
    "..."
  )

  model <- lmer(
    mean_scaled_expr ~ partition * newannotation +
      (1 | replicate) +
      (1 | replicate:newannotation),
    data = df,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 100000
      )
    )
  )

  model_prefix <- file.path(
    output_dir,
    paste0(
      TE_group,
      "_celltype"
    )
  )

  saveRDS(
    model,
    paste0(
      model_prefix,
      "_model.rds"
    )
  )

  write_model_diagnostics(
    model = model,
    output_file = paste0(
      model_prefix,
      "_model_summary.txt"
    ),
    model_name = paste(
      TE_group,
      "EU-HET cell-type model"
    )
  )

  anova_result <- anova(
    model,
    type = 3,
    ddf = "Satterthwaite"
  )

  anova_table <- format_anova_table(
    anova_result
  )

  write.csv(
    anova_table,
    paste0(
      model_prefix,
      "_type3_ANOVA.csv"
    ),
    row.names = FALSE
  )

  emm_partition <- emmeans(
    model,
    ~ partition | newannotation,
    lmer.df = "satterthwaite"
  )

  emm_table <- as.data.frame(
    summary(
      emm_partition,
      infer = c(TRUE, TRUE)
    )
  )

  write.csv(
    emm_table,
    paste0(
      model_prefix,
      "_partition_emmeans.csv"
    ),
    row.names = FALSE
  )

  contrast_table <- as.data.frame(
    summary(
      pairs(
        emm_partition,
        adjust = "none"
      ),
      infer = c(TRUE, TRUE),
      adjust = "none"
    )
  ) %>%
    add_partition_interpretation()

  write.csv(
    contrast_table,
    paste0(
      model_prefix,
      "_EU_vs_HET_globalFDR.csv"
    ),
    row.names = FALSE
  )

  invisible(
    list(
      model = model,
      anova = anova_table,
      emmeans = emm_table,
      contrasts = contrast_table
    )
  )
}


fit_compartment_model <- function(
    df,
    TE_group,
    output_dir
) {

  message(
    "Fitting compartment model for ",
    TE_group,
    "..."
  )

  model <- lmer(
    mean_scaled_expr ~ partition * compartment +
      (1 | replicate) +
      (1 | newannotation),
    data = df,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(
        maxfun = 100000
      )
    )
  )

  model_prefix <- file.path(
    output_dir,
    paste0(
      TE_group,
      "_compartment"
    )
  )

  saveRDS(
    model,
    paste0(
      model_prefix,
      "_model.rds"
    )
  )

  write_model_diagnostics(
    model = model,
    output_file = paste0(
      model_prefix,
      "_model_summary.txt"
    ),
    model_name = paste(
      TE_group,
      "EU-HET compartment model"
    )
  )

  anova_result <- anova(
    model,
    type = 3,
    ddf = "Satterthwaite"
  )

  anova_table <- format_anova_table(
    anova_result
  )

  write.csv(
    anova_table,
    paste0(
      model_prefix,
      "_type3_ANOVA.csv"
    ),
    row.names = FALSE
  )

  emm_partition <- emmeans(
    model,
    ~ partition | compartment,
    lmer.df = "satterthwaite"
  )

  emm_table <- as.data.frame(
    summary(
      emm_partition,
      infer = c(TRUE, TRUE)
    )
  )

  write.csv(
    emm_table,
    paste0(
      model_prefix,
      "_partition_emmeans.csv"
    ),
    row.names = FALSE
  )

  contrast_table <- as.data.frame(
    summary(
      pairs(
        emm_partition,
        adjust = "none"
      ),
      infer = c(TRUE, TRUE),
      adjust = "none"
    )
  ) %>%
    add_partition_interpretation()

  write.csv(
    contrast_table,
    paste0(
      model_prefix,
      "_EU_vs_HET_globalFDR.csv"
    ),
    row.names = FALSE
  )

  invisible(
    list(
      model = model,
      anova = anova_table,
      emmeans = emm_table,
      contrasts = contrast_table
    )
  )
}


# ------------------------------------------------------------------------------
# 5. Run DNA and retrotransposon analyses
# ------------------------------------------------------------------------------

analysis_results <- vector(
  mode = "list",
  length = nrow(input_table)
)

names(analysis_results) <- input_table$TE_group

for (i in seq_len(nrow(input_table))) {

  TE_group <- input_table$TE_group[i]
  input_file <- input_table$input_file[i]

  message(
    "\nProcessing ",
    TE_group,
    " input..."
  )

  analysis_data <- prepare_input_data(
    input_file = input_file,
    TE_group = TE_group
  )

  write.csv(
    analysis_data,
    file.path(
      output_dir,
      paste0(
        TE_group,
        "_validated_input.csv"
      )
    ),
    row.names = FALSE
  )

  celltype_results <- fit_celltype_model(
    df = analysis_data,
    TE_group = TE_group,
    output_dir = output_dir
  )

  compartment_results <- fit_compartment_model(
    df = analysis_data,
    TE_group = TE_group,
    output_dir = output_dir
  )

  analysis_results[[TE_group]] <- list(
    celltype = celltype_results,
    compartment = compartment_results
  )
}


# ------------------------------------------------------------------------------
# 6. Save session information
# ------------------------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    output_dir,
    "sessionInfo_EuHet_mixed_models.txt"
  )
)

message("\nEU-HET mixed-model analyses complete.")
message("Results were written to: ", output_dir)
```
