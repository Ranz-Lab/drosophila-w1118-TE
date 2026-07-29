#!/usr/bin/env Rscript

# Meiotic hdWGCNA and TE-gene co-modularity analysis
#
# Inputs:
#   FCA-testis-merged-normalized-scaled-soloTE-w1118.rds
#   w1118_v4_sorted.protein_coding.gtf
#
# Outputs are written to meiotic_hdWGCNA_results/

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(hdWGCNA)
  library(IRanges)
  library(matrixStats)
  library(Seurat)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(WGCNA)
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

seurat_file <- "FCA-testis-merged-normalized-scaled-soloTE-w1118.rds"
gtf_file <- "w1118_v4_sorted.protein_coding.gtf"
output_dir <- "meiotic_hdWGCNA_results"

wgcna_name <- "germline_testis_TE"
module_prefix <- "Meiotic-M"
feature_fraction <- 0.05
metacell_k <- 25
metacell_max_shared <- 10
module_activity_z <- 0.5
neighbor_distance <- 5000L
genomic_bin_size <- 250000L
n_permutations <- 10000L
minimum_test_pairs <- 5L
seed <- 2L
n_threads <- 8L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(12345)
enableWGCNAThreads(nThreads = n_threads)

mitotic <- c(
  "spermatogonium",
  "mid-late proliferating spermatogonia",
  "spermatogonium-spermatocyte transition"
)

meiotic <- c(
  "spermatocyte 0", "spermatocyte 1", "spermatocyte 2",
  "spermatocyte 3", "spermatocyte 4", "spermatocyte 5",
  "spermatocyte 6", "spermatocyte 7a", "spermatocyte",
  "late primary spermatocyte"
)

postmeiotic <- c(
  "spermatid",
  "early elongation stage spermatid",
  "early-mid elongation-stage spermatid",
  "mid-late elongation-stage spermatid"
)

somatic <- c(
  "secretory cell of the male reproductive tract",
  "pigment cell", "male gonad associated epithelium", "muscle cell",
  "hemocyte", "adult fat body", "adult tracheal cell", "head cyst cell",
  "spermatocyte cyst cell branch b", "spermatocyte cyst cell branch a",
  "late cyst cell branch b", "late cyst cell branch a",
  "cyst cell branch b", "cyst cell branch a", "cyst cell intermediate",
  "early cyst cell 2", "early cyst cell 1", "cyst stem cell"
)

all_cell_types <- c(mitotic, meiotic, postmeiotic, somatic)
stage_levels <- c("somatic", "mitotic", "meiotic", "post-meiotic")

stage_for <- function(x) {
  case_when(
    x %in% somatic ~ "somatic",
    x %in% mitotic ~ "mitotic",
    x %in% meiotic ~ "meiotic",
    x %in% postmeiotic ~ "post-meiotic",
    TRUE ~ NA_character_
  )
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

get_data <- function(object) {
  if (inherits(object[["RNA"]], "Assay5")) {
    GetAssayData(object, assay = "RNA", layer = "data")
  } else {
    GetAssayData(object, assay = "RNA", slot = "data")
  }
}

average_by_celltype <- function(object, annotation_col) {
  tryCatch(
    AverageExpression(
      object, assays = "RNA", group.by = annotation_col,
      layer = "data", verbose = FALSE
    )$RNA,
    error = function(e) {
      AverageExpression(
        object, assays = "RNA", group.by = annotation_col,
        slot = "data", verbose = FALSE
      )$RNA
    }
  )
}

parse_te <- function(features) {
  ids <- gsub("_", "-", features, fixed = TRUE)
  m <- str_match(
    ids,
    paste0(
      "^SoloTE-(X|2L|2R|3L|3R|4|Y)-(\\d+)-(\\d+)-",
      "([^:]+):([^:]+):([^:-]+)-(\\d+(?:\\.\\d+)?)-([+-])$"
    )
  )

  tibble(
    feature = features,
    chromosome = m[, 2],
    start = suppressWarnings(as.integer(m[, 3])),
    end = suppressWarnings(as.integer(m[, 4])),
    TE_family = m[, 5],
    TE_superfamily = m[, 6],
    TE_class = toupper(m[, 7]),
    divergence = suppressWarnings(as.numeric(m[, 8])),
    TE_strand = m[, 9]
  ) %>%
    filter(!is.na(chromosome), !is.na(start), !is.na(end))
}

gtf_attr <- function(x, key) {
  str_match(x, paste0("(?:^|;\\s*)", key, "\\s+\"([^\"]+)\""))[, 2]
}

read_genes <- function(filename, chromosomes) {
  if (!file.exists(filename)) stop("GTF file not found: ", filename)

  gtf <- read.delim(
    filename, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
    quote = "", comment.char = "#", check.names = FALSE,
    col.names = c(
      "seqname", "source", "type", "start", "end",
      "score", "strand", "phase", "attributes"
    )
  ) %>%
    mutate(
      seqname = sub("^chr", "", seqname),
      gene_id = gtf_attr(attributes, "gene_id"),
      gene_name = gtf_attr(attributes, "gene_name")
    ) %>%
    filter(seqname %in% chromosomes)

  genes <- gtf %>%
    filter(type == "gene", !is.na(gene_id)) %>%
    transmute(
      chromosome = seqname, start = as.integer(start), end = as.integer(end),
      strand = as.character(strand), gene_id = as.character(gene_id),
      gene_name = as.character(gene_name)
    )

  if (!nrow(genes)) {
    genes <- gtf %>%
      filter(!is.na(gene_id)) %>%
      group_by(seqname, strand, gene_id) %>%
      summarise(
        start = min(start), end = max(end),
        gene_name = {
          x <- gene_name[!is.na(gene_name) & nzchar(gene_name)]
          if (length(x)) x[1] else NA_character_
        },
        .groups = "drop"
      ) %>%
      transmute(
        chromosome = seqname, start = as.integer(start), end = as.integer(end),
        strand = as.character(strand), gene_id = as.character(gene_id),
        gene_name = as.character(gene_name)
      )
  }

  if (!nrow(genes)) stop("No gene coordinates could be extracted from the GTF.")

  genes %>%
    mutate(
      gene_name = ifelse(
        is.na(gene_name) | !nzchar(gene_name), gene_id, gene_name
      )
    ) %>%
    distinct(chromosome, start, end, strand, gene_id, .keep_all = TRUE)
}

expression_quartile <- function(x) {
  out <- rep(NA_integer_, length(x))
  keep <- which(is.finite(x))
  if (length(keep)) out[keep] <- ntile(x[keep], 4)
  out
}

shuffle_te_modules <- function(x) {
  x %>%
    group_by(TE_family, expression_quartile) %>%
    mutate(TE_module_permuted = sample(TE_module, n(), replace = FALSE)) %>%
    ungroup()
}

permutation_test <- function(x, B = n_permutations, test_seed = seed) {
  x <- x %>%
    filter(
      !is.na(TE_module), !is.na(gene_module),
      !is.na(TE_family), !is.na(expression_quartile)
    )

  if (nrow(x) < minimum_test_pairs) {
    return(list(
      summary = tibble(
        n_pairs = nrow(x), observed_fraction = NA_real_,
        expected_fraction = NA_real_, enrichment = NA_real_,
        null_lower_95 = NA_real_, null_upper_95 = NA_real_,
        p_value = NA_real_
      ),
      null = numeric()
    ))
  }

  observed <- mean(x$TE_module == x$gene_module)
  set.seed(test_seed)

  null <- replicate(B, {
    xp <- shuffle_te_modules(x)
    mean(xp$TE_module_permuted == xp$gene_module)
  })

  interval <- quantile(null, c(0.025, 0.975), na.rm = TRUE)
  expected <- mean(null)

  list(
    summary = tibble(
      n_pairs = nrow(x),
      observed_fraction = observed,
      expected_fraction = expected,
      enrichment = ifelse(expected > 0, observed / expected, NA_real_),
      null_lower_95 = unname(interval[1]),
      null_upper_95 = unname(interval[2]),
      p_value = (1 + sum(null >= observed)) / (B + 1)
    ),
    null = null
  )
}

grouped_tests <- function(x, group_col, seed_offset) {
  groups <- sort(unique(as.character(x[[group_col]])))
  groups <- groups[!is.na(groups) & nzchar(groups)]

  results <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    dat <- x[as.character(x[[group_col]]) == groups[i], , drop = FALSE]
    test <- permutation_test(dat, test_seed = seed + seed_offset + i)
    results[[i]] <- bind_cols(tibble(group = groups[i]), test$summary)
  }

  bind_rows(results) %>%
    mutate(
      p_adjusted_BH = p.adjust(p_value, method = "BH"),
      significant_FDR_0.05 = !is.na(p_adjusted_BH) & p_adjusted_BH < 0.05
    )
}

# ------------------------------------------------------------------------------
# Prepare Seurat object
# ------------------------------------------------------------------------------

if (!file.exists(seurat_file)) stop("Seurat object not found: ", seurat_file)
object <- readRDS(seurat_file)

annotation_candidates <- c("newannotation", "annotation", "S_annotation")
annotation_col <- annotation_candidates[
  annotation_candidates %in% colnames(object[[]])
][1]
if (!length(annotation_col) || is.na(annotation_col)) {
  stop("Expected one of these annotation columns: ",
       paste(annotation_candidates, collapse = ", "))
}

annotations <- as.character(object[[]][[annotation_col]])
names(annotations) <- Cells(object)
keep_cells <- !is.na(annotations) & nzchar(annotations) & annotations != "unannotated"
object <- subset(object, cells = names(annotations)[keep_cells])

annotations <- as.character(object[[]][[annotation_col]])
unknown <- sort(setdiff(unique(annotations), all_cell_types))
if (length(unknown)) {
  stop("Unassigned cell types:\n", paste(unknown, collapse = "\n"))
}

object$germline_stage <- factor(stage_for(annotations), levels = stage_levels)
if (!any(object$germline_stage == "meiotic")) stop("No meiotic cells were found.")

if (!"pca" %in% Reductions(object)) {
  if (!length(VariableFeatures(object))) {
    object <- FindVariableFeatures(
      object, selection.method = "vst", nfeatures = 2000, verbose = FALSE
    )
  }
  object <- ScaleData(
    object, features = VariableFeatures(object), verbose = FALSE
  )
  object <- RunPCA(
    object, features = VariableFeatures(object), npcs = 30, verbose = FALSE
  )
}

# ------------------------------------------------------------------------------
# Meiotic hdWGCNA
# ------------------------------------------------------------------------------

object <- SetupForWGCNA(
  object,
  gene_select = "fraction",
  fraction = feature_fraction,
  wgcna_name = wgcna_name
)

object <- MetacellsByGroups(
  seurat_obj = object,
  group.by = "germline_stage",
  reduction = "pca",
  k = metacell_k,
  max_shared = metacell_max_shared,
  ident.group = "germline_stage"
)

object <- NormalizeMetacells(object)

object <- SetDatExpr(
  object,
  group_name = "meiotic",
  group.by = "germline_stage",
  assay = "RNA",
  layer = "data"
)

object <- TestSoftPowers(object, networkType = "signed")
write.csv(
  GetPowerTable(object),
  file.path(output_dir, "meiotic_soft_power_table.csv"),
  row.names = FALSE
)

object <- ConstructNetwork(object, tom_name = "meiotic")
object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
object <- ModuleEigengenes(object, group.by.vars = NULL)
object <- ModuleConnectivity(
  object, group.by = "germline_stage", group_name = "meiotic"
)
object <- ResetModuleNames(object, new_name = module_prefix)

saveRDS(object, file.path(output_dir, "meiotic_hdWGCNA_object.rds"))

modules <- GetModules(object) %>%
  mutate(
    feature = as.character(gene_name),
    module = as.character(module)
  ) %>%
  filter(module != "grey") %>%
  distinct(feature, .keep_all = TRUE)

write.csv(
  modules,
  file.path(output_dir, "meiotic_modules.csv"),
  row.names = FALSE
)

write.csv(
  GetHubGenes(object, n_hubs = 10),
  file.path(output_dir, "meiotic_hub_features_top10.csv"),
  row.names = FALSE
)

top200 <- GetHubGenes(object, n_hubs = 200) %>%
  mutate(feature = as.character(gene_name)) %>%
  filter(!grepl("^SoloTE[-_]", feature))

write.csv(
  top200,
  file.path(output_dir, "meiotic_top200_coding_genes_by_kME.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Module activity and composition
# ------------------------------------------------------------------------------

MEs <- GetMEs(object, harmonized = TRUE)
common_cells <- intersect(rownames(MEs), Cells(object))
module_names <- unique(modules$module)
ME_cols <- intersect(module_names, colnames(MEs))
if (!length(ME_cols)) {
  ME_cols <- colnames(MEs)[grepl(module_prefix, colnames(MEs), fixed = TRUE)]
}
if (!length(ME_cols)) stop("No renamed meiotic eigengene columns were found.")

ME_z <- scale(as.matrix(MEs[common_cells, ME_cols, drop = FALSE]))
ME_z[is.na(ME_z)] <- 0

ME_long <- as.data.frame(ME_z) %>%
  rownames_to_column("cell") %>%
  pivot_longer(-cell, names_to = "module", values_to = "module_z") %>%
  left_join(
    object[[]] %>%
      rownames_to_column("cell") %>%
      transmute(
        cell,
        cell_type = as.character(.data[[annotation_col]]),
        germline_stage = as.character(germline_stage)
      ),
    by = "cell"
  )

activity_stage <- ME_long %>%
  group_by(germline_stage, module) %>%
  summarise(
    n_cells = n(),
    mean_module_z = mean(module_z),
    fraction_active = mean(module_z > module_activity_z),
    .groups = "drop"
  ) %>%
  mutate(germline_stage = factor(germline_stage, levels = stage_levels)) %>%
  arrange(germline_stage, module)

activity_meiotic <- ME_long %>%
  filter(cell_type %in% meiotic) %>%
  group_by(cell_type, module) %>%
  summarise(
    n_cells = n(),
    mean_module_z = mean(module_z),
    fraction_active = mean(module_z > module_activity_z),
    .groups = "drop"
  ) %>%
  mutate(cell_type = factor(cell_type, levels = meiotic)) %>%
  arrange(cell_type, module)

write.csv(
  activity_stage,
  file.path(output_dir, "meiotic_module_activity_by_stage.csv"),
  row.names = FALSE
)
write.csv(
  activity_meiotic,
  file.path(output_dir, "meiotic_module_activity_by_celltype.csv"),
  row.names = FALSE
)

TE_modules <- parse_te(modules$feature[grepl("^SoloTE[-_]", modules$feature)])

composition <- modules %>%
  transmute(
    feature, module,
    feature_type = ifelse(grepl("^SoloTE[-_]", feature), "TE", "Gene")
  ) %>%
  left_join(
    TE_modules %>% select(feature, TE_family, TE_superfamily, TE_class),
    by = "feature"
  ) %>%
  mutate(feature_class = ifelse(feature_type == "Gene", "Gene", TE_class)) %>%
  count(module, feature_type, feature_class, name = "n_features") %>%
  arrange(module, feature_type, feature_class)

write.csv(
  composition,
  file.path(output_dir, "meiotic_module_composition.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Neighboring TE-gene pairs
# ------------------------------------------------------------------------------

avg <- average_by_celltype(object, annotation_col)
TE_features <- intersect(TE_modules$feature, rownames(avg))
if (!length(TE_features)) stop("No module-assigned TE features were found.")

TE_avg <- avg[TE_features, , drop = FALSE]
TE_z <- t(scale(t(TE_avg)))
TE_z[is.na(TE_z)] <- 0
TE_max_z <- rowMaxs(TE_z)
names(TE_max_z) <- rownames(TE_z)

TEs <- parse_te(TE_features) %>%
  mutate(
    TE_maxZ_all = TE_max_z[feature],
    expression_quartile = expression_quartile(TE_maxZ_all)
  )

genes <- read_genes(gtf_file, unique(TEs$chromosome))

TE_gr <- GRanges(
  TEs$chromosome,
  IRanges(TEs$start, TEs$end),
  strand = TEs$TE_strand
)

gene_gr <- GRanges(
  genes$chromosome,
  IRanges(genes$start, genes$end),
  strand = genes$strand,
  gene_id = genes$gene_id,
  gene_name = genes$gene_name
)

hits <- findOverlaps(TE_gr, gene_gr, ignore.strand = TRUE)
best_gene <- rep(NA_integer_, length(TE_gr))

if (length(hits)) {
  overlap_width <- width(
    pintersect(
      ranges(TE_gr[queryHits(hits)]),
      ranges(gene_gr[subjectHits(hits)])
    )
  )
  best <- tibble(
    TE_index = queryHits(hits),
    gene_index = subjectHits(hits),
    overlap_width
  ) %>%
    group_by(TE_index) %>%
    slice_max(overlap_width, n = 1, with_ties = FALSE) %>%
    ungroup()
  best_gene[best$TE_index] <- best$gene_index
}

nearest <- distanceToNearest(
  resize(TE_gr, width = 1, fix = "center"),
  gene_gr,
  ignore.strand = TRUE
)

nearest_index <- rep(NA_integer_, length(TE_gr))
nearest_distance <- rep(NA_integer_, length(TE_gr))
nearest_index[queryHits(nearest)] <- subjectHits(nearest)
nearest_distance[queryHits(nearest)] <- mcols(nearest)$distance

genic <- !is.na(best_gene)
use_nearest <- (
  !genic & is.finite(nearest_distance) & nearest_distance <= neighbor_distance
)
best_gene[use_nearest] <- nearest_index[use_nearest]

host_gene_id <- host_gene_name <- rep(NA_character_, length(TE_gr))
valid_gene <- !is.na(best_gene)
host_gene_id[valid_gene] <- mcols(gene_gr)$gene_id[best_gene[valid_gene]]
host_gene_name[valid_gene] <- mcols(gene_gr)$gene_name[best_gene[valid_gene]]

orientation <- rep(NA_character_, length(TE_gr))
if (any(genic)) {
  idx <- which(genic)
  te_strand <- as.character(strand(TE_gr[idx]))
  gene_strand <- as.character(strand(gene_gr[best_gene[idx]]))
  valid_strand <- te_strand %in% c("+", "-") & gene_strand %in% c("+", "-")
  orientation[idx[valid_strand]] <- ifelse(
    te_strand[valid_strand] == gene_strand[valid_strand],
    "sense", "antisense"
  )
}

pairs <- TEs %>%
  mutate(
    host_gene_id,
    host_gene_name,
    neighbor_type = case_when(
      genic ~ "genic_overlap",
      use_nearest ~ "nearest_gene_within_5kb",
      TRUE ~ NA_character_
    ),
    neighbor_distance_bp = ifelse(genic, 0, nearest_distance),
    orientation
  ) %>%
  filter(!is.na(neighbor_type))

module_map <- modules %>% select(feature, module)
TE_map <- module_map %>%
  filter(grepl("^SoloTE[-_]", feature)) %>%
  rename(TE_module = module)
gene_map <- module_map %>%
  filter(!grepl("^SoloTE[-_]", feature)) %>%
  rename(gene_feature = feature, gene_module = module)

pairs <- pairs %>%
  left_join(TE_map, by = "feature") %>%
  left_join(gene_map, by = c("host_gene_name" = "gene_feature")) %>%
  rename(module_by_name = gene_module) %>%
  left_join(gene_map, by = c("host_gene_id" = "gene_feature")) %>%
  rename(module_by_id = gene_module) %>%
  mutate(
    gene_module = coalesce(module_by_name, module_by_id),
    gene_match = case_when(
      !is.na(module_by_name) ~ "gene_name",
      !is.na(module_by_id) ~ "gene_id",
      TRUE ~ NA_character_
    ),
    same_module = TE_module == gene_module
  ) %>%
  select(-module_by_name, -module_by_id) %>%
  filter(!is.na(TE_module), !is.na(gene_module)) %>%
  distinct(feature, .keep_all = TRUE)

if (!nrow(pairs)) stop("No neighboring TE-gene pairs had both module labels.")

write.csv(
  pairs,
  file.path(output_dir, "meiotic_neighboring_TE_gene_pairs.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Same-module enrichment tests
# ------------------------------------------------------------------------------

global_test <- permutation_test(pairs)

global_summary <- global_test$summary %>%
  mutate(
    analysis = "All neighboring TE-gene pairs",
    p_adjusted_BH = p_value,
    significant_FDR_0.05 = p_adjusted_BH < 0.05
  ) %>%
  select(analysis, everything())

write.csv(
  global_summary,
  file.path(output_dir, "meiotic_comodularity_global.csv"),
  row.names = FALSE
)
write.csv(
  tibble(
    permutation = seq_along(global_test$null),
    same_module_fraction = global_test$null
  ),
  file.path(output_dir, "meiotic_comodularity_global_null.csv"),
  row.names = FALSE
)

class_results <- grouped_tests(
  pairs %>% filter(TE_class %in% c("DNA", "LINE", "LTR")),
  "TE_class",
  100L
) %>%
  rename(TE_class = group)

write.csv(
  class_results,
  file.path(output_dir, "meiotic_comodularity_by_TE_class.csv"),
  row.names = FALSE
)

family_counts <- pairs %>%
  count(TE_family, name = "n_pairs") %>%
  filter(n_pairs >= minimum_test_pairs)

family_results <- grouped_tests(
  pairs %>% semi_join(family_counts, by = "TE_family"),
  "TE_family",
  200L
) %>%
  rename(TE_family = group) %>%
  left_join(pairs %>% distinct(TE_family, TE_class), by = "TE_family") %>%
  arrange(p_adjusted_BH, TE_class, TE_family)

write.csv(
  family_results,
  file.path(output_dir, "meiotic_comodularity_by_TE_family.csv"),
  row.names = FALSE
)

orientation_pairs <- pairs %>%
  filter(
    neighbor_type == "genic_overlap",
    orientation %in% c("sense", "antisense")
  )

orientation_results <- grouped_tests(
  orientation_pairs,
  "orientation",
  300L
) %>%
  rename(orientation = group)

write.csv(
  orientation_results,
  file.path(output_dir, "meiotic_comodularity_by_orientation.csv"),
  row.names = FALSE
)

orientation_difference <- tibble(
  n_sense = sum(orientation_pairs$orientation == "sense"),
  n_antisense = sum(orientation_pairs$orientation == "antisense"),
  observed_sense_fraction = NA_real_,
  observed_antisense_fraction = NA_real_,
  observed_difference = NA_real_,
  null_mean_difference = NA_real_,
  null_lower_95 = NA_real_,
  null_upper_95 = NA_real_,
  p_value_two_sided = NA_real_
)

if (
  orientation_difference$n_sense >= minimum_test_pairs &&
  orientation_difference$n_antisense >= minimum_test_pairs
) {
  sense_obs <- mean(
    orientation_pairs$same_module[orientation_pairs$orientation == "sense"]
  )
  antisense_obs <- mean(
    orientation_pairs$same_module[orientation_pairs$orientation == "antisense"]
  )
  observed_diff <- sense_obs - antisense_obs

  set.seed(seed + 400L)
  null_diff <- replicate(n_permutations, {
    xp <- shuffle_te_modules(orientation_pairs)
    mean(
      xp$TE_module_permuted[xp$orientation == "sense"] ==
      xp$gene_module[xp$orientation == "sense"]
    ) -
    mean(
      xp$TE_module_permuted[xp$orientation == "antisense"] ==
      xp$gene_module[xp$orientation == "antisense"]
    )
  })

  interval <- quantile(null_diff, c(0.025, 0.975))

  orientation_difference <- tibble(
    n_sense = orientation_difference$n_sense,
    n_antisense = orientation_difference$n_antisense,
    observed_sense_fraction = sense_obs,
    observed_antisense_fraction = antisense_obs,
    observed_difference = observed_diff,
    null_mean_difference = mean(null_diff),
    null_lower_95 = unname(interval[1]),
    null_upper_95 = unname(interval[2]),
    p_value_two_sided = (
      1 + sum(abs(null_diff) >= abs(observed_diff))
    ) / (n_permutations + 1)
  )
}

write.csv(
  orientation_difference,
  file.path(output_dir, "meiotic_orientation_difference_test.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Genomic density of module-assigned features
# ------------------------------------------------------------------------------

TE_density_features <- parse_te(
  modules$feature[grepl("^SoloTE[-_]", modules$feature)]
) %>%
  left_join(module_map, by = "feature")

gene_modules <- module_map %>%
  filter(!grepl("^SoloTE[-_]", feature))

module_genes <- genes %>%
  left_join(
    gene_modules %>% rename(module_by_name = module),
    by = c("gene_name" = "feature")
  ) %>%
  left_join(
    gene_modules %>% rename(module_by_id = module),
    by = c("gene_id" = "feature")
  ) %>%
  mutate(module = coalesce(module_by_name, module_by_id)) %>%
  filter(!is.na(module)) %>%
  distinct(gene_id, module, .keep_all = TRUE)

TE_density_gr <- GRanges(
  TE_density_features$chromosome,
  IRanges(TE_density_features$start, TE_density_features$end),
  module = TE_density_features$module
)

gene_density_gr <- GRanges(
  module_genes$chromosome,
  IRanges(module_genes$start, module_genes$end),
  module = module_genes$module
)

chromosome_lengths <- bind_rows(
  TE_density_features %>%
    group_by(chromosome) %>%
    summarise(max_end = max(end), .groups = "drop"),
  genes %>%
    group_by(chromosome) %>%
    summarise(max_end = max(end), .groups = "drop")
) %>%
  group_by(chromosome) %>%
  summarise(length = max(max_end), .groups = "drop") %>%
  filter(chromosome %in% c("2L", "2R", "3L", "3R", "4", "X", "Y"))

make_bins <- function(chromosome, length, bin_size) {
  bin_start <- seq(1, length, by = bin_size)
  bin_end <- pmin(bin_start + bin_size - 1, length)
  GRanges(
    chromosome,
    IRanges(bin_start, bin_end),
    bin_id = paste0(chromosome, ":", bin_start, "-", bin_end),
    bin_midpoint = (bin_start + bin_end) / 2
  )
}

bins <- do.call(
  c,
  lapply(seq_len(nrow(chromosome_lengths)), function(i) {
    make_bins(
      chromosome_lengths$chromosome[i],
      chromosome_lengths$length[i],
      genomic_bin_size
    )
  })
)

count_bins <- function(features, feature_type) {
  hits <- findOverlaps(bins, features, ignore.strand = TRUE)
  tibble(
    chromosome = as.character(seqnames(bins)[queryHits(hits)]),
    bin_id = mcols(bins)$bin_id[queryHits(hits)],
    bin_start = start(bins)[queryHits(hits)],
    bin_end = end(bins)[queryHits(hits)],
    bin_midpoint = mcols(bins)$bin_midpoint[queryHits(hits)],
    module = mcols(features)$module[subjectHits(hits)],
    feature_type
  ) %>%
    count(
      chromosome, bin_id, bin_start, bin_end, bin_midpoint,
      module, feature_type, name = "n_features"
    ) %>%
    mutate(density_per_Mb = n_features / (genomic_bin_size / 1e6))
}

density_table <- bind_rows(
  count_bins(gene_density_gr, "Gene"),
  count_bins(TE_density_gr, "TE")
) %>%
  arrange(feature_type, module, chromosome, bin_start)

write.csv(
  density_table,
  file.path(output_dir, "meiotic_module_genomic_density_250kb.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Reproducibility information
# ------------------------------------------------------------------------------

parameters <- tibble(
  parameter = c(
    "Seurat object", "GTF annotation", "Feature detection fraction",
    "Metacell k", "Metacell max_shared", "Network type",
    "Network compartment", "Module prefix", "Neighbor definition",
    "Permutation strata", "Permutations", "Genomic bin size"
  ),
  value = c(
    seurat_file, gtf_file, feature_fraction,
    metacell_k, metacell_max_shared, "signed",
    "meiotic", module_prefix,
    paste0("genic overlap or nearest gene <= ", neighbor_distance, " bp"),
    "TE family and expression quartile",
    n_permutations, genomic_bin_size
  )
)

write.csv(
  parameters,
  file.path(output_dir, "analysis_parameters.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_meiotic_hdWGCNA.txt")
)

message("Meiotic hdWGCNA analysis complete.")
message("Results written to: ", output_dir)
