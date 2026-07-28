#!/usr/bin/env Rscript

# ==============================================================================
# TE class expression across germline and somatic compartments
#
# This script:
#   1. Assigns annotated testis cell types to germline or somatic compartments
#   2. Collapses replicate-level TE-class expression across cell types
#   3. Fits a linear mixed-effects model:
#
#        mean_expr ~ TE_class * compartment + (1 | replicate)
#
#   4. Tests fixed effects using a Type III ANOVA
#   5. Performs pairwise comparisons:
#        - TE classes within each compartment
#        - compartments within each TE class
#   6. Applies Benjamini-Hochberg correction across each set of comparisons
#
# Input:
#   INPUT_LME_replicate_celltype_scaled_TEclass.csv
#
# Required input columns:
#   replicate, newannotation, TE_class, mean_expr
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
# 2. Input and output files
# ------------------------------------------------------------------------------

input_file <- "INPUT_LME_replicate_celltype_scaled_TEclass.csv"

collapsed_output_file <-
  "INPUT_LME_replicate_compartment_scaled_TEclass.csv"

anova_output_file <-
  "LME_ANOVA_TEclass_compartment.csv"

emmeans_te_output_file <-
  "LME_emmeans_TEclass_within_compartment.csv"

pairwise_te_output_file <-
  "LME_pairwise_TEclass_within_compartment_globalFDR.csv"

emmeans_compartment_output_file <-
  "LME_emmeans_compartment_within_TEclass.csv"

pairwise_compartment_output_file <-
  "LME_pairwise_compartment_within_TEclass_globalFDR.csv"

model_output_file <-
  "LME_TEclass_compartment_model.rds"

model_summary_output_file <-
  "LME_TEclass_compartment_model_summary.txt"

session_info_output_file <-
  "sessionInfo_TEclass_compartment_LME.txt"


# ------------------------------------------------------------------------------
# 3. Define testis cell-type compartments
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
  germline_cell_types,
  somatic_cell_types
)

te_class_levels <- c(
  "DNA",
  "LINE",
  "LTR"
)

compartment_levels <- c(
  "Somatic",
  "Germline"
)


# ------------------------------------------------------------------------------
# 4. Load and validate input data
# ------------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

df_rep_long <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "replicate",
  "newannotation",
  "TE_class",
  "mean_expr"
)

missing_columns <- setdiff(
  required_columns,
  colnames(df_rep_long)
)

if (length(missing_columns) > 0) {
  stop(
    "The input file is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(df_rep_long) == 0) {
  stop("The input file contains no rows.")
}

df_rep_long$mean_expr <- suppressWarnings(
  as.numeric(df_rep_long$mean_expr)
)

if (any(!is.finite(df_rep_long$mean_expr))) {
  stop(
    "The mean_expr column contains missing, infinite, or non-numeric values."
  )
}

unknown_cell_types <- sort(
  setdiff(
    unique(df_rep_long$newannotation),
    known_cell_types
  )
)

if (length(unknown_cell_types) > 0) {
  stop(
    "The following cell types are not assigned to a compartment:\n",
    paste(unknown_cell_types, collapse = "\n")
  )
}

unknown_te_classes <- sort(
  setdiff(
    unique(df_rep_long$TE_class),
    te_class_levels
  )
)

if (length(unknown_te_classes) > 0) {
  stop(
    "Unexpected TE classes found in the input data: ",
    paste(unknown_te_classes, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 5. Assign cell types to compartments
# ------------------------------------------------------------------------------

df_rep_long <- df_rep_long %>%
  mutate(
    compartment = case_when(
      newannotation %in% somatic_cell_types ~ "Somatic",
      newannotation %in% germline_cell_types ~ "Germline",
      TRUE ~ NA_character_
    ),
    replicate = factor(replicate),
    compartment = factor(
      compartment,
      levels = compartment_levels
    ),
    TE_class = factor(
      TE_class,
      levels = te_class_levels
    )
  )

if (nlevels(df_rep_long$replicate) < 2) {
  stop("At least two biological replicates are required.")
}


# ------------------------------------------------------------------------------
# 6. Collapse to replicate x compartment x TE class
# ------------------------------------------------------------------------------

df_compartment <- df_rep_long %>%
  group_by(
    replicate,
    compartment,
    TE_class
  ) %>%
  summarise(
    mean_expr = mean(mean_expr),
    n_cell_types = n_distinct(newannotation),
    .groups = "drop"
  ) %>%
  arrange(
    replicate,
    compartment,
    TE_class
  )

expected_combinations <- expand.grid(
  replicate = levels(df_rep_long$replicate),
  compartment = compartment_levels,
  TE_class = te_class_levels,
  stringsAsFactors = FALSE
)

observed_combinations <- df_compartment %>%
  transmute(
    replicate = as.character(replicate),
    compartment = as.character(compartment),
    TE_class = as.character(TE_class)
  )

missing_combinations <- anti_join(
  expected_combinations,
  observed_combinations,
  by = c(
    "replicate",
    "compartment",
    "TE_class"
  )
)

if (nrow(missing_combinations) > 0) {
  stop(
    "The collapsed dataset has missing replicate-compartment-TE class ",
    "combinations:\n",
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
  df_compartment,
  collapsed_output_file,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. Fit the linear mixed-effects model
# ------------------------------------------------------------------------------

m_lme_compartment <- lmer(
  mean_expr ~ TE_class * compartment + (1 | replicate),
  data = df_compartment,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 100000)
  )
)

saveRDS(
  m_lme_compartment,
  model_output_file
)

convergence_messages <-
  m_lme_compartment@optinfo$conv$lme4$messages

if (!is.null(convergence_messages)) {
  warning(
    "Model convergence message: ",
    paste(convergence_messages, collapse = "; ")
  )
}

singular_fit <- isSingular(
  m_lme_compartment,
  tol = 1e-4
)

if (singular_fit) {
  warning(
    "The fitted model is singular; inspect the estimated replicate variance."
  )
}


# ------------------------------------------------------------------------------
# 8. Type III ANOVA
# ------------------------------------------------------------------------------

anova_res <- anova(
  m_lme_compartment,
  type = 3,
  ddf = "Satterthwaite"
)

anova_table <- as.data.frame(anova_res) %>%
  tibble::rownames_to_column("term")

write.csv(
  anova_table,
  anova_output_file,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 9. Compare TE classes within each compartment
# ------------------------------------------------------------------------------

emm_te_by_compartment <- emmeans(
  m_lme_compartment,
  ~ TE_class | compartment
)

emm_te_table <- as.data.frame(
  emm_te_by_compartment
)

write.csv(
  emm_te_table,
  emmeans_te_output_file,
  row.names = FALSE
)

pairwise_te_by_compartment <- pairs(
  emm_te_by_compartment,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    padj_fdr_global = p.adjust(
      p.value,
      method = "BH"
    ),
    interpretation = case_when(
      contrast == "DNA - LINE" & estimate < 0 ~ "LINE > DNA",
      contrast == "DNA - LINE" & estimate > 0 ~ "DNA > LINE",
      contrast == "DNA - LTR"  & estimate < 0 ~ "LTR > DNA",
      contrast == "DNA - LTR"  & estimate > 0 ~ "DNA > LTR",
      contrast == "LINE - LTR" & estimate < 0 ~ "LTR > LINE",
      contrast == "LINE - LTR" & estimate > 0 ~ "LINE > LTR",
      TRUE ~ "No difference"
    ),
    significant_fdr_0.05 = padj_fdr_global < 0.05
  )

write.csv(
  pairwise_te_by_compartment,
  pairwise_te_output_file,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 10. Compare compartments within each TE class
# ------------------------------------------------------------------------------

emm_compartment_by_te <- emmeans(
  m_lme_compartment,
  ~ compartment | TE_class
)

emm_compartment_table <- as.data.frame(
  emm_compartment_by_te
)

write.csv(
  emm_compartment_table,
  emmeans_compartment_output_file,
  row.names = FALSE
)

pairwise_compartment_by_te <- pairs(
  emm_compartment_by_te,
  adjust = "none"
) %>%
  as.data.frame() %>%
  mutate(
    padj_fdr_global = p.adjust(
      p.value,
      method = "BH"
    ),
    interpretation = case_when(
      contrast == "Somatic - Germline" & estimate < 0 ~
        "Germline > Somatic",
      contrast == "Somatic - Germline" & estimate > 0 ~
        "Somatic > Germline",
      TRUE ~ "No difference"
    ),
    significant_fdr_0.05 = padj_fdr_global < 0.05
  )

write.csv(
  pairwise_compartment_by_te,
  pairwise_compartment_output_file,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 11. Save model summary and session information
# ------------------------------------------------------------------------------

summary_lines <- c(
  capture.output(summary(m_lme_compartment)),
  "",
  paste0("Singular fit: ", singular_fit),
  "",
  "Type III ANOVA:",
  capture.output(print(anova_table)),
  "",
  "TE-class comparisons within compartments:",
  capture.output(print(pairwise_te_by_compartment)),
  "",
  "Compartment comparisons within TE classes:",
  capture.output(print(pairwise_compartment_by_te))
)

writeLines(
  summary_lines,
  model_summary_output_file
)

writeLines(
  capture.output(sessionInfo()),
  session_info_output_file
)

message("Mixed-model analysis complete.")
message("Model: ", model_output_file)
message("ANOVA: ", anova_output_file)
message("TE-class contrasts: ", pairwise_te_output_file)
message("Compartment contrasts: ", pairwise_compartment_output_file)
