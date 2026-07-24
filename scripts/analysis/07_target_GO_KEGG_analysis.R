###############################################################################
# RBZ STEP16 MASTER CLEAN PIPELINE
# Downstream target, overlap, shared-target, network, pre-hub, and enrichment
# analysis for locked miRNAs
#
# This script consolidates the revised workflow used in the project:
#   16A  Prepare locked miRNA input
#   16B  Collect predicted targets with multiMiR
#   16C  Collect validated targets with multiMiR
#   16D  Harmonize gene symbols
#   16E  Build consensus predicted targets (>= 2 databases)
#   16F  Build overlap-supported targets
#   16G  Identify shared targets + Venn diagram
#   16H  Prepare miRNA-target network files
#   16I  Prepare PPI input and pre-hub candidates
#   16K  Run GO/KEGG enrichment (clean final version)
#
# Notes:
# - Step 16J (final PPI hub genes from external STRING/Cytoscape exports) is
#   not included here because it depends on files generated outside R.
# - This script preserves the folder structure used in the project.
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

###############################################################################
# 0) USER SETTINGS
###############################################################################
PROJECT_ROOT <- getwd()
INPUT_FILE_RAW <- file.path(PROJECT_ROOT, "Locked_miRNAs.csv")

# Main behavior switches
RUN_16A <- TRUE
RUN_16B_16C <- TRUE
RUN_16D <- TRUE
RUN_16E <- TRUE
RUN_16F <- TRUE
RUN_16G <- TRUE
RUN_16G_VENN <- TRUE
RUN_16H <- TRUE
RUN_16I <- TRUE
RUN_16K <- TRUE

# Step 16B/C multiMiR settings
USE_SUMMARY_RESULT    <- FALSE
USE_TIBBLE_OUTPUT     <- FALSE
LIMIT_PER_QUERY       <- NULL
PREDICTED_SITE_MODE   <- "all"
PREDICTED_CUTOFF      <- NULL
PREDICTED_CUTOFF_TYPE <- "p"

# Step 16E consensus rule
MIN_DATABASE_SUPPORT <- 2

# Step 16K enrichment settings
GO_PVALUE_CUTOFF   <- 0.05
GO_QVALUE_CUTOFF   <- 0.20
KEGG_PVALUE_CUTOFF <- 0.05
KEGG_QVALUE_CUTOFF <- 0.20
TOP_N_PLOT         <- 15
HUB_TOP_N_FALLBACK <- 10

###############################################################################
# 1) OUTPUT FOLDER STRUCTURE
###############################################################################
dirs_to_create <- c(
  "00_Input",
  "01_Target_Prediction",
  "01_Target_Prediction/Raw_By_Database",
  "01_Target_Prediction/Per_miRNA",
  "01_Target_Prediction/Tables",
  "02_Target_Validation",
  "02_Target_Validation/Raw_Validated",
  "02_Target_Validation/Tables",
  "03_Harmonized_Targets",
  "04_Consensus_Predicted_Targets",
  "05_Overlap_Targets",
  "05_Overlap_Targets/Per_miRNA",
  "05_Overlap_Targets/Tables",
  "06_Shared_Targets",
  "06_Shared_Targets/Plots",
  "07_miRNA_Target_Network",
  "07_miRNA_Target_Network/Cytoscape_Files",
  "07_miRNA_Target_Network/Edge_Lists",
  "07_miRNA_Target_Network/Node_Attributes",
  "07_miRNA_Target_Network/Plots",
  "08_PPI_Network",
  "08_PPI_Network/Edge_Lists",
  "08_PPI_Network/Node_Attributes",
  "08_PPI_Network/Hub_Genes",
  "08_PPI_Network/Plots",
  "09_Enrichment",
  "09_Enrichment/GO_BP",
  "09_Enrichment/GO_MF",
  "09_Enrichment/GO_CC",
  "09_Enrichment/KEGG",
  "09_Enrichment/Tables",
  "10_Integrated_Figures",
  "11_Methods_and_Logs",
  "12_Scripts"
)
for (d in dirs_to_create) dir.create(file.path(PROJECT_ROOT, d), recursive = TRUE, showWarnings = FALSE)

MASTER_LOG <- file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Run_Log_STEP16_MASTER.txt")

###############################################################################
# 2) HELPER FUNCTIONS
###############################################################################
log_message <- function(..., log_file = MASTER_LOG) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  if (!is.null(log_file)) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    cat(msg, "\n", file = log_file, append = TRUE)
  }
}

safe_write_csv <- function(df, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (is.null(df)) df <- data.frame()
  write.csv(df, file, row.names = FALSE, quote = TRUE)
}

safe_write_txt_vector <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  writeLines(x, con = file, sep = "\n")
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null", "na", "n/a")] <- NA_character_
  x
}

clean_mirna_names <- function(x) {
  x <- trimws(x)
  x <- gsub('"', "", x)
  x <- gsub("^[[:space:]]+|[[:space:]]+$", "", x)
  x
}

sanitize_filename <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

collapse_unique_sorted <- function(x, sep = "; ") {
  x <- clean_text(x)
  x <- unique(x[!is.na(x)])
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

get_first_existing_column <- function(df, candidates) {
  hits <- intersect(candidates, colnames(df))
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

extract_target_gene_column <- function(df) {
  get_first_existing_column(df, c(
    "target_symbol", "target_gene", "target_entrez", "target_ensembl",
    "target_id", "target_alias", "gene_symbol", "target"
  ))
}
extract_mirna_column <- function(df) {
  get_first_existing_column(df, c("mature_mirna_id", "mature_mirna_acc", "mirna", "mature_mirna"))
}
extract_database_column <- function(df) {
  get_first_existing_column(df, c("database", "table", "source"))
}
extract_score_column <- function(df) {
  get_first_existing_column(df, c("score", "predicted_score", "binding", "probability", "context"))
}
extract_validation_method_column <- function(df) {
  get_first_existing_column(df, c("experiment", "support_type", "validation_method"))
}
extract_pubmed_column <- function(df) {
  get_first_existing_column(df, c("pubmed_id", "pubmed", "pmid"))
}

empty_predicted_df <- function() {
  data.frame(
    miRNA = character(), gene_symbol = character(), database = character(),
    score = character(), source_type = character(), stringsAsFactors = FALSE
  )
}
empty_validated_df <- function() {
  data.frame(
    miRNA = character(), gene_symbol = character(), database = character(),
    validation_method = character(), evidence_strength = character(),
    PMID = character(), source_type = character(), stringsAsFactors = FALSE
  )
}

harmonize_gene_symbol <- function(x) {
  x <- clean_text(x)
  x <- toupper(x)
  x <- gsub("^\\s+|\\s+$", "", x)
  x <- gsub("[[:space:]]+", "", x)
  x
}

is_obviously_invalid_gene <- function(x) {
  x0 <- clean_text(x)
  x1 <- harmonize_gene_symbol(x0)
  invalid_set <- c(
    NA, "", "-", ".", "...", "NAN", "NONE", "UNKNOWN", "UNCHARACTERIZED",
    "UNDEFINED", "MIRNA", "MICRORNA", "TARGET", "GENE", "NA_CHARACTER_"
  )
  bad_basic <- is.na(x1) | x1 %in% invalid_set
  bad_numeric_only <- grepl("^[0-9._-]+$", x1)
  bad_phrase <- grepl("[;,:/\\\\]", x0)
  bad_basic | bad_numeric_only | bad_phrase
}

consensus_class_from_n <- function(n) {
  ifelse(is.na(n), NA_character_,
         ifelse(n >= 4, "Very strong",
                ifelse(n == 3, "Strong",
                       ifelse(n == 2, "Moderate", "Low"))))
}

make_indicator_name <- function(mir) paste0("targeted_by_", gsub("[^A-Za-z0-9_]", "_", mir))

parse_gene_ratio_numeric <- function(x) {
  x <- as.character(x)
  sapply(x, function(v) {
    if (is.na(v) || !grepl("/", v, fixed = TRUE)) return(NA_real_)
    parts <- strsplit(v, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NA_real_)
    num <- suppressWarnings(as.numeric(parts[1]))
    den <- suppressWarnings(as.numeric(parts[2]))
    if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
    num / den
  })
}

standardize_enrichment_table <- function(res_df, source_gene_set, analysis_type) {
  if (is.null(res_df) || nrow(res_df) == 0) {
    return(data.frame(
      term_id = character(), term_name = character(), gene_ratio = character(),
      bg_ratio = character(), p_value = numeric(), adjusted_p_value = numeric(),
      qvalue = numeric(), gene_count = numeric(), gene_list = character(),
      source_gene_set = character(), analysis_type = character(),
      stringsAsFactors = FALSE
    ))
  }
  needed <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue",
              "p.adjust", "qvalue", "Count", "geneID")
  for (nm in needed) if (!nm %in% colnames(res_df)) res_df[[nm]] <- NA
  data.frame(
    term_id = as.character(res_df$ID),
    term_name = as.character(res_df$Description),
    gene_ratio = as.character(res_df$GeneRatio),
    bg_ratio = as.character(res_df$BgRatio),
    p_value = as.numeric(res_df$pvalue),
    adjusted_p_value = as.numeric(res_df$p.adjust),
    qvalue = as.numeric(res_df$qvalue),
    gene_count = as.numeric(res_df$Count),
    gene_list = as.character(res_df$geneID),
    source_gene_set = source_gene_set,
    analysis_type = analysis_type,
    stringsAsFactors = FALSE
  )
}

prepare_plot_df <- function(df, top_n = 15) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  plot_df <- df
  plot_df <- plot_df[!is.na(plot_df$adjusted_p_value), , drop = FALSE]
  if (nrow(plot_df) == 0) return(NULL)
  plot_df$gene_ratio_numeric <- parse_gene_ratio_numeric(plot_df$gene_ratio)
  plot_df$neglog10_padj <- -log10(plot_df$adjusted_p_value)
  plot_df <- plot_df[order(plot_df$adjusted_p_value, -plot_df$gene_count), , drop = FALSE]
  plot_df <- utils::head(plot_df, top_n)
  if (nrow(plot_df) == 0) return(NULL)
  plot_df$term_name <- factor(plot_df$term_name, levels = rev(plot_df$term_name))
  plot_df
}

save_colorful_dotplot <- function(df, out_prefix, title_text, top_n = 15) {
  plot_df <- prepare_plot_df(df, top_n)
  if (is.null(plot_df)) return(invisible(NULL))
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = gene_ratio_numeric, y = term_name, size = gene_count, color = neglog10_padj)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient(low = "#2C7BB6", high = "#D7191C",
                                  name = expression(-log[10]("adj. p"))) +
    ggplot2::labs(title = title_text, x = "Gene ratio", y = NULL, size = "Gene count") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   axis.text.y = ggplot2::element_text(size = 10),
                   legend.title = ggplot2::element_text(face = "bold"))
  grDevices::pdf(paste0(out_prefix, ".pdf"), width = 10, height = 7); print(p); grDevices::dev.off()
  grDevices::png(paste0(out_prefix, ".png"), width = 3000, height = 2100, res = 300); print(p); grDevices::dev.off()
}

save_colorful_barplot <- function(df, out_prefix, title_text, top_n = 15) {
  plot_df <- prepare_plot_df(df, top_n)
  if (is.null(plot_df)) return(invisible(NULL))
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = gene_count, y = term_name, fill = neglog10_padj)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient(low = "#ABD9E9", high = "#D73027",
                                 name = expression(-log[10]("adj. p"))) +
    ggplot2::labs(title = title_text, x = "Gene count", y = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   axis.text.y = ggplot2::element_text(size = 10),
                   legend.title = ggplot2::element_text(face = "bold"))
  grDevices::pdf(paste0(out_prefix, ".pdf"), width = 10, height = 7); print(p); grDevices::dev.off()
  grDevices::png(paste0(out_prefix, ".png"), width = 3000, height = 2100, res = 300); print(p); grDevices::dev.off()
}

extract_unique_genes <- function(file, gene_col = "gene_symbol", top_n = NULL) {
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!gene_col %in% colnames(df)) stop("Gene column '", gene_col, "' not found in: ", file)
  genes <- toupper(clean_text(df[[gene_col]]))
  genes <- genes[!is.na(genes)]
  genes <- unique(genes)
  if (!is.null(top_n)) genes <- utils::head(genes, top_n)
  genes
}

log_message("STEP16 master pipeline started.")

###############################################################################
# 3) PACKAGE SETUP
###############################################################################
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran_pkgs <- c("dplyr", "tibble", "tidyr", "ggplot2", "VennDiagram")
bioc_pkgs <- c("multiMiR", "clusterProfiler", "org.Hs.eg.db")
for (pkg in cran_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
for (pkg in bioc_pkgs) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
  library(multiMiR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

###############################################################################
# STEP 16A
###############################################################################
if (RUN_16A) {
  log_message("Running Step 16A")
  INPUT_FILE <- INPUT_FILE_RAW
  locked_raw <- read.csv(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(locked_raw)[1] <- "miRNA"
  locked_raw$miRNA_original <- locked_raw$miRNA
  locked_raw$miRNA <- clean_mirna_names(locked_raw$miRNA)
  locked_clean <- locked_raw[nzchar(locked_raw$miRNA), , drop = FALSE]
  dup_flag <- duplicated(locked_clean$miRNA)
  duplicates_table <- locked_clean[dup_flag, c("miRNA_original", "miRNA"), drop = FALSE]
  locked_clean <- locked_clean[!duplicated(locked_clean$miRNA), , drop = FALSE]
  valid_pattern <- "^hsa-miR-[A-Za-z0-9.-]+-(3p|5p)$"
  locked_clean$valid_mature_miRNA_name <- grepl(valid_pattern, locked_clean$miRNA)
  invalid_table <- locked_clean[!locked_clean$valid_mature_miRNA_name, , drop = FALSE]
  locked_final <- locked_clean[order(locked_clean$miRNA), c("miRNA"), drop = FALSE]
  rownames(locked_final) <- NULL
  safe_write_csv(locked_final, file.path(PROJECT_ROOT, "00_Input", "Locked_miRNAs.csv"))
  safe_write_csv(locked_clean, file.path(PROJECT_ROOT, "00_Input", "Locked_miRNAs_Audit.csv"))
  safe_write_csv(duplicates_table, file.path(PROJECT_ROOT, "00_Input", "Locked_miRNAs_Duplicates_Removed.csv"))
  safe_write_csv(invalid_table, file.path(PROJECT_ROOT, "00_Input", "Locked_miRNAs_Invalid_Names.csv"))
  summary_tbl <- data.frame(
    metric = c("n_input_rows_raw", "n_nonempty_rows_after_cleaning", "n_duplicates_removed",
               "n_final_unique_miRNAs", "n_invalid_miRNA_names"),
    value = c(nrow(locked_raw), nrow(locked_clean), nrow(duplicates_table),
              nrow(locked_final), nrow(invalid_table)),
    stringsAsFactors = FALSE
  )
  safe_write_csv(summary_tbl, file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Locked_miRNAs_Summary.csv"))
  sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16A.txt")); print(sessionInfo()); sink()
}

###############################################################################
# STEP 16B + 16C
###############################################################################
if (RUN_16B_16C) {
  log_message("Running Step 16B + 16C")
  INPUT_FILE <- file.path(PROJECT_ROOT, "00_Input", "Locked_miRNAs.csv")
  locked_mirs <- read.csv(INPUT_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  locked_mirs$miRNA <- clean_text(locked_mirs$miRNA)
  locked_mirs <- locked_mirs[nzchar(locked_mirs$miRNA), , drop = FALSE]
  locked_mirs <- locked_mirs[!duplicated(locked_mirs$miRNA), , drop = FALSE]

  safe_write_csv(as.data.frame(tryCatch(multimir_dbInfo(), error = function(e) NULL)),
                 file.path(PROJECT_ROOT, "11_Methods_and_Logs", "multiMiR_DB_Info.csv"))
  safe_write_csv(as.data.frame(tryCatch(multimir_dbTables(), error = function(e) NULL)),
                 file.path(PROJECT_ROOT, "11_Methods_and_Logs", "multiMiR_DB_Tables.csv"))
  safe_write_csv(as.data.frame(tryCatch(multimir_dbInfoVersions(), error = function(e) NULL)),
                 file.path(PROJECT_ROOT, "11_Methods_and_Logs", "multiMiR_DB_Versions.csv"))

  run_multimir_query <- function(one_mirna, query_table = c("predicted", "validated")) {
    query_table <- match.arg(query_table)
    res <- tryCatch(
      get_multimir(
        org = "hsa",
        mirna = one_mirna,
        table = query_table,
        predicted.cutoff = if (query_table == "predicted") PREDICTED_CUTOFF else NULL,
        predicted.cutoff.type = if (query_table == "predicted") PREDICTED_CUTOFF_TYPE else "p",
        predicted.site = if (query_table == "predicted") PREDICTED_SITE_MODE else "all",
        summary = USE_SUMMARY_RESULT,
        use.tibble = USE_TIBBLE_OUTPUT,
        limit = LIMIT_PER_QUERY,
        legacy.out = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(res)) return(NULL)
    tryCatch(as.data.frame(res@data), error = function(e) NULL)
  }

  pred_long_list <- list(); val_long_list <- list(); pred_summary_rows <- list(); val_summary_rows <- list()
  for (i in seq_len(nrow(locked_mirs))) {
    one_mir <- locked_mirs$miRNA[i]
    mir_tag <- sanitize_filename(one_mir)

    pred_raw <- run_multimir_query(one_mir, "predicted")
    if (is.null(pred_raw) || nrow(pred_raw) == 0) {
      pred_long_keep <- empty_predicted_df()
      safe_write_csv(pred_long_keep, file.path(PROJECT_ROOT, "01_Target_Prediction", "Raw_By_Database", paste0("multiMiR_predicted_", mir_tag, ".csv")))
      safe_write_csv(pred_long_keep, file.path(PROJECT_ROOT, "01_Target_Prediction", "Per_miRNA", paste0("PredictedTargets_", mir_tag, ".csv")))
      pred_long_list[[one_mir]] <- pred_long_keep
      pred_summary_rows[[one_mir]] <- data.frame(miRNA = one_mir, source_type = "predicted", n_rows_raw = 0,
                                                 n_unique_gene_symbols = 0, n_unique_databases = 0, stringsAsFactors = FALSE)
    } else {
      pred_raw$queried_miRNA <- rep(one_mir, nrow(pred_raw))
      pred_raw$source_type <- rep("predicted", nrow(pred_raw))
      safe_write_csv(pred_raw, file.path(PROJECT_ROOT, "01_Target_Prediction", "Raw_By_Database", paste0("multiMiR_predicted_", mir_tag, ".csv")))
      mir_col <- extract_mirna_column(pred_raw)
      gene_col <- extract_target_gene_column(pred_raw)
      db_col <- extract_database_column(pred_raw)
      score_col <- extract_score_column(pred_raw)
      pred_long <- pred_raw
      pred_long$miRNA <- if (!is.na(mir_col)) clean_text(pred_long[[mir_col]]) else one_mir
      pred_long$gene_symbol <- if (!is.na(gene_col)) clean_text(pred_long[[gene_col]]) else NA_character_
      pred_long$database <- if (!is.na(db_col)) clean_text(pred_long[[db_col]]) else "multiMiR_predicted"
      pred_long$score <- if (!is.na(score_col)) as.character(pred_long[[score_col]]) else NA_character_
      pred_long_keep <- pred_long %>% dplyr::mutate(miRNA = queried_miRNA, gene_symbol = gene_symbol,
                                                    database = database, score = score, source_type = "predicted")
      pred_long_list[[one_mir]] <- pred_long_keep
      safe_write_csv(pred_long_keep, file.path(PROJECT_ROOT, "01_Target_Prediction", "Per_miRNA", paste0("PredictedTargets_", mir_tag, ".csv")))
      pred_summary_rows[[one_mir]] <- data.frame(miRNA = one_mir, source_type = "predicted",
                                                 n_rows_raw = nrow(pred_raw),
                                                 n_unique_gene_symbols = length(unique(na.omit(pred_long_keep$gene_symbol))),
                                                 n_unique_databases = length(unique(na.omit(pred_long_keep$database))),
                                                 stringsAsFactors = FALSE)
    }

    val_raw <- run_multimir_query(one_mir, "validated")
    if (is.null(val_raw) || nrow(val_raw) == 0) {
      val_long_keep <- empty_validated_df()
      safe_write_csv(val_long_keep, file.path(PROJECT_ROOT, "02_Target_Validation", "Raw_Validated", paste0("multiMiR_validated_", mir_tag, ".csv")))
      safe_write_csv(val_long_keep, file.path(PROJECT_ROOT, "02_Target_Validation", "Tables", paste0("ValidatedTargets_", mir_tag, ".csv")))
      val_long_list[[one_mir]] <- val_long_keep
      val_summary_rows[[one_mir]] <- data.frame(miRNA = one_mir, source_type = "validated", n_rows_raw = 0,
                                                n_unique_gene_symbols = 0, n_unique_databases = 0, stringsAsFactors = FALSE)
    } else {
      val_raw$queried_miRNA <- rep(one_mir, nrow(val_raw))
      val_raw$source_type <- rep("validated", nrow(val_raw))
      safe_write_csv(val_raw, file.path(PROJECT_ROOT, "02_Target_Validation", "Raw_Validated", paste0("multiMiR_validated_", mir_tag, ".csv")))
      mir_col <- extract_mirna_column(val_raw)
      gene_col <- extract_target_gene_column(val_raw)
      db_col <- extract_database_column(val_raw)
      method_col <- extract_validation_method_column(val_raw)
      pmid_col <- extract_pubmed_column(val_raw)
      val_long <- val_raw
      val_long$miRNA <- if (!is.na(mir_col)) clean_text(val_long[[mir_col]]) else one_mir
      val_long$gene_symbol <- if (!is.na(gene_col)) clean_text(val_long[[gene_col]]) else NA_character_
      val_long$database <- if (!is.na(db_col)) clean_text(val_long[[db_col]]) else "multiMiR_validated"
      val_long$validation_method <- if (!is.na(method_col)) clean_text(val_long[[method_col]]) else NA_character_
      val_long$PMID <- if (!is.na(pmid_col)) clean_text(val_long[[pmid_col]]) else NA_character_
      val_long$evidence_strength <- NA_character_
      val_long_keep <- val_long %>% dplyr::mutate(miRNA = queried_miRNA, gene_symbol = gene_symbol,
                                                  database = database, validation_method = validation_method,
                                                  evidence_strength = evidence_strength, PMID = PMID,
                                                  source_type = "validated")
      val_long_list[[one_mir]] <- val_long_keep
      safe_write_csv(val_long_keep, file.path(PROJECT_ROOT, "02_Target_Validation", "Tables", paste0("ValidatedTargets_", mir_tag, ".csv")))
      val_summary_rows[[one_mir]] <- data.frame(miRNA = one_mir, source_type = "validated",
                                                n_rows_raw = nrow(val_raw),
                                                n_unique_gene_symbols = length(unique(na.omit(val_long_keep$gene_symbol))),
                                                n_unique_databases = length(unique(na.omit(val_long_keep$database))),
                                                stringsAsFactors = FALSE)
    }
  }

  pred_all <- dplyr::bind_rows(pred_long_list)
  val_all <- dplyr::bind_rows(val_long_list)
  pred_summary <- dplyr::bind_rows(pred_summary_rows)
  val_summary <- dplyr::bind_rows(val_summary_rows)
  safe_write_csv(pred_all, file.path(PROJECT_ROOT, "01_Target_Prediction", "Tables", "Table_All_Predicted_Targets_Long.csv"))
  safe_write_csv(val_all, file.path(PROJECT_ROOT, "02_Target_Validation", "Tables", "Table_All_Validated_Targets_Long.csv"))
  pred_compact <- if (nrow(pred_all) > 0) pred_all %>% dplyr::select(any_of(c("miRNA","gene_symbol","database","score","source_type"))) else empty_predicted_df()
  val_compact <- if (nrow(val_all) > 0) val_all %>% dplyr::select(any_of(c("miRNA","gene_symbol","database","validation_method","evidence_strength","PMID","source_type"))) else empty_validated_df()
  safe_write_csv(pred_compact, file.path(PROJECT_ROOT, "01_Target_Prediction", "Tables", "Table_All_Predicted_Targets_Compact.csv"))
  safe_write_csv(val_compact, file.path(PROJECT_ROOT, "02_Target_Validation", "Tables", "Table_All_Validated_Targets_Compact.csv"))
  safe_write_csv(dplyr::bind_rows(pred_summary, val_summary), file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Step16B_16C_Summary.csv"))
  sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16B_16C.txt")); print(sessionInfo()); sink()
}

###############################################################################
# STEP 16D
###############################################################################
if (RUN_16D) {
  log_message("Running Step 16D")
  pred_raw <- read.csv(file.path(PROJECT_ROOT, "01_Target_Prediction", "Tables", "Table_All_Predicted_Targets_Compact.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  val_raw <- read.csv(file.path(PROJECT_ROOT, "02_Target_Validation", "Tables", "Table_All_Validated_Targets_Compact.csv"), check.names = FALSE, stringsAsFactors = FALSE)

  standardize_predicted <- function(df) {
    if (nrow(df) == 0) return(data.frame(miRNA=character(), gene_symbol_original=character(), gene_symbol=character(), database=character(), score=character(), source_type=character(), stringsAsFactors = FALSE))
    needed <- c("miRNA","gene_symbol","database","score","source_type")
    for (nm in needed) if (!nm %in% colnames(df)) df[[nm]] <- NA_character_
    data.frame(miRNA = clean_text(df$miRNA), gene_symbol_original = clean_text(df$gene_symbol), gene_symbol = harmonize_gene_symbol(df$gene_symbol), database = clean_text(df$database), score = as.character(df$score), source_type = clean_text(df$source_type), stringsAsFactors = FALSE)
  }
  standardize_validated <- function(df) {
    if (nrow(df) == 0) return(data.frame(miRNA=character(), gene_symbol_original=character(), gene_symbol=character(), database=character(), validation_method=character(), evidence_strength=character(), PMID=character(), source_type=character(), stringsAsFactors = FALSE))
    needed <- c("miRNA","gene_symbol","database","validation_method","evidence_strength","PMID","source_type")
    for (nm in needed) if (!nm %in% colnames(df)) df[[nm]] <- NA_character_
    data.frame(miRNA = clean_text(df$miRNA), gene_symbol_original = clean_text(df$gene_symbol), gene_symbol = harmonize_gene_symbol(df$gene_symbol), database = clean_text(df$database), validation_method = clean_text(df$validation_method), evidence_strength = clean_text(df$evidence_strength), PMID = clean_text(df$PMID), source_type = clean_text(df$source_type), stringsAsFactors = FALSE)
  }
  pred_std <- standardize_predicted(pred_raw)
  val_std <- standardize_validated(val_raw)
  pred_drop_flag <- is_obviously_invalid_gene(pred_std$gene_symbol_original) | is.na(pred_std$miRNA) | pred_std$miRNA == ""
  val_drop_flag <- is_obviously_invalid_gene(val_std$gene_symbol_original) | is.na(val_std$miRNA) | val_std$miRNA == ""
  pred_dropped <- pred_std[pred_drop_flag, , drop = FALSE]
  val_dropped <- val_std[val_drop_flag, , drop = FALSE]
  pred_keep <- pred_std[!pred_drop_flag, , drop = FALSE]
  val_keep <- val_std[!val_drop_flag, , drop = FALSE]
  pred_dup_flag <- duplicated(pred_keep[, c("miRNA", "gene_symbol", "database")])
  val_key_cols <- c("miRNA", "gene_symbol", "database", "validation_method", "PMID")
  for (nm in val_key_cols) if (!nm %in% colnames(val_keep)) val_keep[[nm]] <- NA_character_
  val_dup_flag <- duplicated(val_keep[, val_key_cols])
  pred_harm <- pred_keep[!pred_dup_flag, , drop = FALSE]
  val_harm <- val_keep[!val_dup_flag, , drop = FALSE]
  safe_write_csv(pred_harm, file.path(PROJECT_ROOT, "03_Harmonized_Targets", "Table_All_Predicted_Targets_Harmonized.csv"))
  safe_write_csv(val_harm, file.path(PROJECT_ROOT, "03_Harmonized_Targets", "Table_All_Validated_Targets_Harmonized.csv"))
  safe_write_csv(pred_dropped, file.path(PROJECT_ROOT, "03_Harmonized_Targets", "Dropped_Rows_Predicted.csv"))
  safe_write_csv(val_dropped, file.path(PROJECT_ROOT, "03_Harmonized_Targets", "Dropped_Rows_Validated.csv"))
  safe_write_csv(pred_keep[pred_dup_flag, , drop = FALSE], file.path(PROJECT_ROOT, "03_Harmonized_Targets", "Duplicates_Removed_Predicted.csv"))
  safe_write_csv(val_keep[val_dup_flag, , drop = FALSE], file.path(PROJECT_ROOT, "03_Harmonized_Targets", "Duplicates_Removed_Validated.csv"))
  harm_log <- dplyr::bind_rows(
    if (nrow(pred_harm)>0) unique(data.frame(dataset="predicted", miRNA=pred_harm$miRNA, gene_symbol_original=pred_harm$gene_symbol_original, gene_symbol_harmonized=pred_harm$gene_symbol, changed=ifelse(pred_harm$gene_symbol_original != pred_harm$gene_symbol, "yes", "no"), stringsAsFactors = FALSE)) else data.frame(),
    if (nrow(val_harm)>0) unique(data.frame(dataset="validated", miRNA=val_harm$miRNA, gene_symbol_original=val_harm$gene_symbol_original, gene_symbol_harmonized=val_harm$gene_symbol, changed=ifelse(val_harm$gene_symbol_original != val_harm$gene_symbol, "yes", "no"), stringsAsFactors = FALSE)) else data.frame()
  )
  safe_write_csv(harm_log, file.path(PROJECT_ROOT, "03_Harmonized_Targets", "Table_GeneSymbol_Harmonization_Log.csv"))
  summary_tbl <- data.frame(dataset = c("predicted", "validated"),
                            n_input_rows = c(nrow(pred_raw), nrow(val_raw)),
                            n_after_standardization = c(nrow(pred_std), nrow(val_std)),
                            n_dropped_invalid_or_blank = c(nrow(pred_dropped), nrow(val_dropped)),
                            n_after_drop = c(nrow(pred_keep), nrow(val_keep)),
                            n_duplicates_removed = c(sum(pred_dup_flag), sum(val_dup_flag)),
                            n_final_rows = c(nrow(pred_harm), nrow(val_harm)),
                            n_unique_miRNAs = c(length(unique(pred_harm$miRNA)), length(unique(val_harm$miRNA))),
                            n_unique_genes = c(length(unique(pred_harm$gene_symbol)), length(unique(val_harm$gene_symbol))),
                            stringsAsFactors = FALSE)
  safe_write_csv(summary_tbl, file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Step16D_Harmonization_Summary.csv"))
  sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16D.txt")); print(sessionInfo()); sink()
}

###############################################################################
# STEP 16E
###############################################################################
if (RUN_16E) {
  log_message("Running Step 16E")
  pred <- read.csv(file.path(PROJECT_ROOT, "03_Harmonized_Targets", "Table_All_Predicted_Targets_Harmonized.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  pred$miRNA <- clean_text(pred$miRNA); pred$gene_symbol <- clean_text(pred$gene_symbol); pred$database <- clean_text(pred$database); pred$score <- clean_text(pred$score)
  pred <- pred[!is.na(pred$miRNA) & !is.na(pred$gene_symbol) & !is.na(pred$database), , drop = FALSE]
  db_usage <- pred %>% dplyr::group_by(database) %>% dplyr::summarise(n_rows = dplyr::n(), n_unique_miRNAs = dplyr::n_distinct(miRNA), n_unique_genes = dplyr::n_distinct(gene_symbol), .groups = "drop") %>% dplyr::arrange(dplyr::desc(n_rows), database)
  safe_write_csv(db_usage, file.path(PROJECT_ROOT, "04_Consensus_Predicted_Targets", "Table_Prediction_Database_Usage_Summary.csv"))
  pred_distinct <- pred %>% dplyr::distinct(miRNA, gene_symbol, database, .keep_all = TRUE)
  support_tbl <- pred_distinct %>% dplyr::group_by(miRNA, gene_symbol) %>% dplyr::summarise(
    n_prediction_databases = dplyr::n_distinct(database),
    supporting_databases = collapse_unique_sorted(database),
    consensus_class = consensus_class_from_n(dplyr::n_distinct(database)),
    .groups = "drop"
  ) %>% dplyr::arrange(miRNA, dplyr::desc(n_prediction_databases), gene_symbol)
  pair_row_counts <- pred %>% dplyr::group_by(miRNA, gene_symbol) %>% dplyr::summarise(n_rows_all_sources = dplyr::n(), .groups = "drop")
  support_tbl <- support_tbl %>% dplyr::left_join(pair_row_counts, by = c("miRNA","gene_symbol")) %>% dplyr::select(miRNA, gene_symbol, n_prediction_databases, supporting_databases, consensus_class, n_rows_all_sources)
  consensus_tbl <- support_tbl %>% dplyr::filter(n_prediction_databases >= MIN_DATABASE_SUPPORT)
  rejected_tbl <- support_tbl %>% dplyr::filter(n_prediction_databases < MIN_DATABASE_SUPPORT)
  per_mir_summary <- support_tbl %>% dplyr::group_by(miRNA) %>% dplyr::summarise(
    n_all_predicted_pairs = dplyr::n(),
    n_low_support_1db = sum(n_prediction_databases == 1, na.rm = TRUE),
    n_moderate_2db = sum(n_prediction_databases == 2, na.rm = TRUE),
    n_strong_3db = sum(n_prediction_databases == 3, na.rm = TRUE),
    n_very_strong_4plus = sum(n_prediction_databases >= 4, na.rm = TRUE),
    n_consensus_pairs_kept = sum(n_prediction_databases >= MIN_DATABASE_SUPPORT, na.rm = TRUE),
    .groups = "drop"
  ) %>% dplyr::arrange(miRNA)
  safe_write_csv(support_tbl, file.path(PROJECT_ROOT, "04_Consensus_Predicted_Targets", "Table_All_Predicted_Targets_DatabaseSupport.csv"))
  safe_write_csv(consensus_tbl, file.path(PROJECT_ROOT, "04_Consensus_Predicted_Targets", "Table_Consensus_Predicted_Targets.csv"))
  safe_write_csv(rejected_tbl, file.path(PROJECT_ROOT, "04_Consensus_Predicted_Targets", "Table_Rejected_NonConsensus_Predicted_Targets.csv"))
  safe_write_csv(per_mir_summary, file.path(PROJECT_ROOT, "04_Consensus_Predicted_Targets", "Table_Consensus_Predicted_Targets_By_miRNA.csv"))
  summary_tbl <- data.frame(
    metric = c("n_rows_input_harmonized_predicted","n_distinct_miRNA_gene_database_rows","n_unique_miRNA_gene_pairs_all","n_consensus_pairs_kept_support_ge_2","n_rejected_pairs_support_lt_2","n_unique_miRNAs_in_all_pairs","n_unique_miRNAs_in_consensus_pairs","n_unique_genes_in_all_pairs","n_unique_genes_in_consensus_pairs"),
    value = c(nrow(pred), nrow(pred_distinct), nrow(support_tbl), nrow(consensus_tbl), nrow(rejected_tbl), length(unique(support_tbl$miRNA)), if (nrow(consensus_tbl)>0) length(unique(consensus_tbl$miRNA)) else 0, length(unique(support_tbl$gene_symbol)), if (nrow(consensus_tbl)>0) length(unique(consensus_tbl$gene_symbol)) else 0),
    stringsAsFactors = FALSE
  )
  safe_write_csv(summary_tbl, file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Step16E_Consensus_Summary.csv"))
  sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16E.txt")); print(sessionInfo()); sink()
}

###############################################################################
# STEP 16F
###############################################################################
if (RUN_16F) {
  log_message("Running Step 16F")
  cons <- read.csv(file.path(PROJECT_ROOT, "04_Consensus_Predicted_Targets", "Table_Consensus_Predicted_Targets.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  val <- read.csv(file.path(PROJECT_ROOT, "03_Harmonized_Targets", "Table_All_Validated_Targets_Harmonized.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  cons$miRNA <- clean_text(cons$miRNA); cons$gene_symbol <- clean_text(cons$gene_symbol); cons$supporting_databases <- clean_text(cons$supporting_databases); cons$consensus_class <- clean_text(cons$consensus_class)
  val$miRNA <- clean_text(val$miRNA); val$gene_symbol <- clean_text(val$gene_symbol); val$database <- clean_text(val$database); val$validation_method <- clean_text(val$validation_method); val$PMID <- clean_text(val$PMID)
  cons <- cons[!is.na(cons$miRNA) & !is.na(cons$gene_symbol), , drop = FALSE]
  val <- val[!is.na(val$miRNA) & !is.na(val$gene_symbol), , drop = FALSE]
  val_collapsed <- val %>% dplyr::group_by(miRNA, gene_symbol) %>% dplyr::summarise(
    validated_databases = collapse_unique_sorted(database),
    validation_methods = collapse_unique_sorted(validation_method),
    PMIDs = collapse_unique_sorted(PMID),
    n_validated_rows = dplyr::n(), .groups = "drop")
  overlap_tbl <- cons %>% dplyr::inner_join(val_collapsed, by = c("miRNA","gene_symbol")) %>% dplyr::mutate(
    prediction_support_n = n_prediction_databases,
    predicted_databases = supporting_databases,
    evidence_class = ifelse(prediction_support_n >= 3, "Very high confidence", ifelse(prediction_support_n >= 2, "High confidence", "Other"))
  ) %>% dplyr::select(miRNA, gene_symbol, prediction_support_n, predicted_databases, consensus_class, n_rows_all_sources, validated_databases, validation_methods, PMIDs, evidence_class) %>% dplyr::arrange(miRNA, dplyr::desc(prediction_support_n), gene_symbol)
  safe_write_csv(overlap_tbl, file.path(PROJECT_ROOT, "05_Overlap_Targets", "Tables", "Table_All_Overlap_Targets.csv"))
  all_miRNAs <- sort(unique(cons$miRNA))
  overlap_summary_list <- list(); retained_list <- list(); excluded_list <- list()
  for (mir in all_miRNAs) {
    this_cons <- cons %>% dplyr::filter(miRNA == mir)
    this_val <- val_collapsed %>% dplyr::filter(miRNA == mir)
    this_ov <- overlap_tbl %>% dplyr::filter(miRNA == mir)
    safe_write_csv(this_ov, file.path(PROJECT_ROOT, "05_Overlap_Targets", "Per_miRNA", paste0("Table_Overlap_Targets_", sanitize_filename(mir), ".csv")))
    overlap_summary_list[[mir]] <- data.frame(
      miRNA = mir,
      n_consensus_predicted_targets = nrow(this_cons),
      n_validated_targets = nrow(this_val),
      n_overlap_targets = nrow(this_ov),
      retained_for_downstream_target_biology = ifelse(nrow(this_ov) > 0, "yes", "no"),
      exclusion_reason = ifelse(nrow(this_ov) > 0, NA_character_, "No overlap between consensus predicted and validated targets"),
      stringsAsFactors = FALSE
    )
    if (nrow(this_ov) > 0) {
      retained_list[[mir]] <- data.frame(miRNA = mir, n_overlap_targets = nrow(this_ov), downstream_status = "retained", rationale = "Has overlap-supported target set", stringsAsFactors = FALSE)
    } else {
      excluded_list[[mir]] <- data.frame(miRNA = mir, n_overlap_targets = 0, downstream_status = "excluded_from_overlap_core", rationale = "Insufficient overlap-supported target evidence", stringsAsFactors = FALSE)
    }
  }
  overlap_summary <- dplyr::bind_rows(overlap_summary_list)
  retained_tbl <- dplyr::bind_rows(retained_list)
  excluded_tbl <- dplyr::bind_rows(excluded_list)
  if (nrow(retained_tbl) == 0) retained_tbl <- data.frame(miRNA=character(), n_overlap_targets=numeric(), downstream_status=character(), rationale=character(), stringsAsFactors = FALSE)
  if (nrow(excluded_tbl) == 0) excluded_tbl <- data.frame(miRNA=character(), n_overlap_targets=numeric(), downstream_status=character(), rationale=character(), stringsAsFactors = FALSE)
  safe_write_csv(overlap_summary, file.path(PROJECT_ROOT, "05_Overlap_Targets", "Tables", "Table_Overlap_Targets_Summary_By_miRNA.csv"))
  safe_write_csv(retained_tbl, file.path(PROJECT_ROOT, "05_Overlap_Targets", "Tables", "Table_miRNAs_Retained_For_Downstream_Target_Biology.csv"))
  safe_write_csv(excluded_tbl, file.path(PROJECT_ROOT, "05_Overlap_Targets", "Tables", "Table_Excluded_miRNAs_From_Overlap_Downstream.csv"))
  summary_tbl <- data.frame(metric = c("n_consensus_predicted_rows_input", "n_validated_collapsed_pairs_input", "n_total_overlap_rows", "n_miRNAs_in_consensus_layer", "n_miRNAs_with_overlap_targets", "n_miRNAs_excluded_from_overlap_core", "n_unique_genes_in_overlap_targets"),
                            value = c(nrow(cons), nrow(val_collapsed), nrow(overlap_tbl), length(unique(cons$miRNA)), sum(overlap_summary$retained_for_downstream_target_biology == "yes", na.rm = TRUE), sum(overlap_summary$retained_for_downstream_target_biology == "no", na.rm = TRUE), if (nrow(overlap_tbl)>0) length(unique(overlap_tbl$gene_symbol)) else 0), stringsAsFactors = FALSE)
  safe_write_csv(summary_tbl, file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Step16F_Overlap_Summary.csv"))
  sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16F.txt")); print(sessionInfo()); sink()
}

###############################################################################
# STEP 16G + VENN
###############################################################################
if (RUN_16G) {
  log_message("Running Step 16G")
  ov <- read.csv(file.path(PROJECT_ROOT, "05_Overlap_Targets", "Tables", "Table_All_Overlap_Targets.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  retained <- read.csv(file.path(PROJECT_ROOT, "05_Overlap_Targets", "Tables", "Table_miRNAs_Retained_For_Downstream_Target_Biology.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  ov$miRNA <- clean_text(ov$miRNA); ov$gene_symbol <- clean_text(ov$gene_symbol); retained$miRNA <- clean_text(retained$miRNA)
  retained_mirs <- sort(unique(retained$miRNA[!is.na(retained$miRNA)]))
  ov <- ov %>% dplyr::filter(!is.na(miRNA), !is.na(gene_symbol), miRNA %in% retained_mirs) %>% dplyr::distinct(miRNA, gene_symbol, .keep_all = TRUE)
  membership_long <- ov %>% dplyr::select(miRNA, gene_symbol) %>% dplyr::mutate(value = 1L)
  membership_wide <- membership_long %>% tidyr::pivot_wider(names_from = miRNA, values_from = value, values_fill = 0)
  for (mir in retained_mirs) if (!mir %in% colnames(membership_wide)) membership_wide[[mir]] <- 0L
  membership_wide <- membership_wide[, c("gene_symbol", retained_mirs), drop = FALSE]
  membership_tbl <- data.frame(gene_symbol = membership_wide$gene_symbol, stringsAsFactors = FALSE)
  for (mir in retained_mirs) membership_tbl[[make_indicator_name(mir)]] <- as.integer(membership_wide[[mir]])
  membership_tbl$n_miRNAs_targeting <- rowSums(membership_wide[, retained_mirs, drop = FALSE])
  membership_tbl$miRNA_list <- apply(membership_wide[, retained_mirs, drop = FALSE], 1, function(v) {
    hits <- retained_mirs[which(v == 1)]
    if (length(hits) == 0) return(NA_character_)
    paste(hits, collapse = "; ")
  })
  membership_tbl$target_class <- ifelse(membership_tbl$n_miRNAs_targeting == 1, "unique", ifelse(membership_tbl$n_miRNAs_targeting == length(retained_mirs), "shared_by_all_retained_miRNAs", "shared_by_subset"))
  membership_tbl <- membership_tbl %>% dplyr::arrange(dplyr::desc(n_miRNAs_targeting), gene_symbol)
  unique_tbl <- membership_tbl %>% dplyr::filter(n_miRNAs_targeting == 1)
  shared_2plus_tbl <- membership_tbl %>% dplyr::filter(n_miRNAs_targeting >= 2)
  shared_all_tbl <- membership_tbl %>% dplyr::filter(n_miRNAs_targeting == length(retained_mirs))
  safe_write_csv(membership_tbl, file.path(PROJECT_ROOT, "06_Shared_Targets", "Table_Shared_Targets_Membership.csv"))
  safe_write_csv(unique_tbl, file.path(PROJECT_ROOT, "06_Shared_Targets", "Table_Unique_Targets.csv"))
  safe_write_csv(shared_2plus_tbl, file.path(PROJECT_ROOT, "06_Shared_Targets", "Table_Shared_Targets_2plus.csv"))
  safe_write_csv(shared_all_tbl, file.path(PROJECT_ROOT, "06_Shared_Targets", "Table_Shared_Targets_All_Retained_miRNAs.csv"))
  per_mir_summary_list <- list()
  for (mir in retained_mirs) {
    indicator_col <- make_indicator_name(mir)
    per_mir_summary_list[[mir]] <- data.frame(
      miRNA = mir,
      n_overlap_targets_total = sum(membership_tbl[[indicator_col]] == 1, na.rm = TRUE),
      n_unique_targets = sum(membership_tbl[[indicator_col]] == 1 & membership_tbl$n_miRNAs_targeting == 1, na.rm = TRUE),
      n_shared_targets_2plus = sum(membership_tbl[[indicator_col]] == 1 & membership_tbl$n_miRNAs_targeting >= 2, na.rm = TRUE),
      n_shared_targets_all_retained = sum(membership_tbl[[indicator_col]] == 1 & membership_tbl$n_miRNAs_targeting == length(retained_mirs), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  per_mir_summary <- dplyr::bind_rows(per_mir_summary_list) %>% dplyr::arrange(miRNA)
  safe_write_csv(per_mir_summary, file.path(PROJECT_ROOT, "06_Shared_Targets", "Table_Shared_Targets_By_miRNA.csv"))
  shared_targets_summary <- data.frame(category = c("unique","shared_by_2plus","shared_by_all_retained_miRNAs"), n_targets = c(nrow(unique_tbl), nrow(shared_2plus_tbl), nrow(shared_all_tbl)), stringsAsFactors = FALSE)
  safe_write_csv(shared_targets_summary, file.path(PROJECT_ROOT, "06_Shared_Targets", "Table_Shared_Targets_Summary.csv"))
  summary_tbl <- data.frame(metric = c("n_retained_miRNAs","n_total_gene_membership_rows","n_unique_targets","n_shared_targets_2plus","n_shared_targets_all_retained_miRNAs"),
                            value = c(length(retained_mirs), nrow(membership_tbl), nrow(unique_tbl), nrow(shared_2plus_tbl), nrow(shared_all_tbl)), stringsAsFactors = FALSE)
  safe_write_csv(summary_tbl, file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Step16G_SharedTargets_Summary.csv"))
  sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16G.txt")); print(sessionInfo()); sink()

  if (RUN_16G_VENN && length(retained_mirs) >= 2 && length(retained_mirs) <= 3) {
    log_message("Running Step 16G Venn")
    target_sets <- lapply(retained_mirs, function(mir) sort(unique(ov$gene_symbol[ov$miRNA == mir])))
    names(target_sets) <- retained_mirs
    set_summary <- data.frame(miRNA = retained_mirs, n_targets = vapply(target_sets, length, integer(1)), stringsAsFactors = FALSE)
    safe_write_csv(set_summary, file.path(PROJECT_ROOT, "06_Shared_Targets", "Plots", "Table_Venn_Set_Sizes.csv"))
    if (length(retained_mirs) == 2) {
      a <- retained_mirs[1]; b <- retained_mirs[2]
      venn_summary <- data.frame(metric = c(paste0("n_", a), paste0("n_", b), paste0("n_only_", a), paste0("n_only_", b), paste0("n_shared_", a, "_", b)),
                                 value = c(length(target_sets[[a]]), length(target_sets[[b]]), length(setdiff(target_sets[[a]], target_sets[[b]])), length(setdiff(target_sets[[b]], target_sets[[a]])), length(intersect(target_sets[[a]], target_sets[[b]]))), stringsAsFactors = FALSE)
      venn_grob <- venn.diagram(x = target_sets, category.names = retained_mirs, filename = NULL,
                                fill = c("skyblue", "tomato"), alpha = 0.5, cex = 1.4,
                                cat.cex = 1.2, cat.pos = c(-20, 20), cat.dist = c(0.05, 0.05), margin = 0.1)
    } else {
      a <- retained_mirs[1]; b <- retained_mirs[2]; c <- retained_mirs[3]
      venn_summary <- data.frame(metric = c(paste0("n_", a), paste0("n_", b), paste0("n_", c), paste0("n_shared_", a, "_", b), paste0("n_shared_", a, "_", c), paste0("n_shared_", b, "_", c), "n_shared_all"),
                                 value = c(length(target_sets[[a]]), length(target_sets[[b]]), length(target_sets[[c]]), length(intersect(target_sets[[a]], target_sets[[b]])), length(intersect(target_sets[[a]], target_sets[[c]])), length(intersect(target_sets[[b]], target_sets[[c]])), length(Reduce(intersect, target_sets))), stringsAsFactors = FALSE)
      venn_grob <- venn.diagram(x = target_sets, category.names = retained_mirs, filename = NULL,
                                fill = c("skyblue", "tomato", "palegreen3"), alpha = 0.5,
                                cex = 1.2, cat.cex = 1.1, margin = 0.1)
    }
    safe_write_csv(venn_summary, file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Step16G_Venn_Summary.csv"))
    png(file.path(PROJECT_ROOT, "06_Shared_Targets", "Plots", "Figure_Venn_Shared_Targets.png"), width = 2400, height = 2000, res = 300)
    grid.newpage(); grid.draw(venn_grob); dev.off()
    pdf(file.path(PROJECT_ROOT, "06_Shared_Targets", "Plots", "Figure_Venn_Shared_Targets.pdf"), width = 8, height = 7)
    grid.newpage(); grid.draw(venn_grob); dev.off()
    sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16G_Venn.txt")); print(sessionInfo()); sink()
  }
}

###############################################################################
# STEP 16H
###############################################################################
if (RUN_16H) {
  log_message("Running Step 16H")
  ov <- read.csv(file.path(PROJECT_ROOT, "05_Overlap_Targets", "Tables", "Table_All_Overlap_Targets.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  membership <- read.csv(file.path(PROJECT_ROOT, "06_Shared_Targets", "Table_Shared_Targets_Membership.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  retained <- read.csv(file.path(PROJECT_ROOT, "05_Overlap_Targets", "Tables", "Table_miRNAs_Retained_For_Downstream_Target_Biology.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  ov$miRNA <- clean_text(ov$miRNA); ov$gene_symbol <- clean_text(ov$gene_symbol); membership$gene_symbol <- clean_text(membership$gene_symbol); membership$miRNA_list <- clean_text(membership$miRNA_list); membership$target_class <- clean_text(membership$target_class); retained$miRNA <- clean_text(retained$miRNA)
  retained_mirs <- sort(unique(retained$miRNA[!is.na(retained$miRNA)]))
  ov <- ov %>% dplyr::filter(!is.na(miRNA), !is.na(gene_symbol), miRNA %in% retained_mirs) %>% dplyr::distinct(miRNA, gene_symbol, .keep_all = TRUE)
  membership <- membership %>% dplyr::filter(!is.na(gene_symbol)) %>% dplyr::distinct(gene_symbol, .keep_all = TRUE)
  edge_full <- ov %>% dplyr::left_join(membership %>% dplyr::select(gene_symbol, n_miRNAs_targeting, miRNA_list, target_class), by = "gene_symbol") %>% dplyr::mutate(source = miRNA, target = gene_symbol, interaction = "miRNA_targets_gene", support_n = prediction_support_n) %>% dplyr::select(source, target, interaction, support_n, evidence_class, predicted_databases, validated_databases, validation_methods, PMIDs, n_miRNAs_targeting, miRNA_list, target_class) %>% dplyr::arrange(source, dplyr::desc(support_n), target)
  safe_write_csv(edge_full, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Edge_Lists", "Table_miRNA_Target_Edges.csv"))
  safe_write_csv(edge_full, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Cytoscape_Files", "Cytoscape_EdgeList_Full.csv"))
  miRNA_nodes <- data.frame(node = retained_mirs, node_type = "miRNA", n_miRNAs_targeting = NA_real_, evidence_class = "locked_miRNA", target_class = NA_character_, stringsAsFactors = FALSE)
  gene_nodes <- edge_full %>% dplyr::group_by(target) %>% dplyr::summarise(n_miRNAs_targeting = dplyr::first(n_miRNAs_targeting), evidence_class = collapse_unique_sorted(evidence_class), target_class = dplyr::first(target_class), .groups = "drop") %>% dplyr::rename(node = target) %>% dplyr::mutate(node_type = "gene") %>% dplyr::select(node, node_type, n_miRNAs_targeting, evidence_class, target_class)
  degree_from_sources <- edge_full %>% dplyr::count(source, name = "degree_source") %>% dplyr::rename(node = source)
  degree_from_targets <- edge_full %>% dplyr::count(target, name = "degree_target") %>% dplyr::rename(node = target)
  node_full <- dplyr::bind_rows(miRNA_nodes, gene_nodes) %>% dplyr::left_join(degree_from_sources, by = "node") %>% dplyr::left_join(degree_from_targets, by = "node") %>% dplyr::mutate(degree_total = ifelse(is.na(degree_source), 0, degree_source) + ifelse(is.na(degree_target), 0, degree_target)) %>% dplyr::select(node, node_type, n_miRNAs_targeting, evidence_class, target_class, degree_total) %>% dplyr::arrange(node_type, dplyr::desc(degree_total), node)
  safe_write_csv(node_full, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Node_Attributes", "Table_miRNA_Target_Nodes.csv"))
  safe_write_csv(node_full, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Cytoscape_Files", "Cytoscape_NodeAttributes_Full.csv"))
  shared_genes <- membership %>% dplyr::filter(!is.na(n_miRNAs_targeting), n_miRNAs_targeting >= 2) %>% dplyr::pull(gene_symbol) %>% unique() %>% sort()
  edge_shared <- edge_full %>% dplyr::filter(target %in% shared_genes) %>% dplyr::arrange(source, dplyr::desc(support_n), target)
  if (nrow(edge_shared) == 0) {
    node_shared <- data.frame(node=character(), node_type=character(), n_miRNAs_targeting=numeric(), evidence_class=character(), target_class=character(), degree_total=numeric(), stringsAsFactors = FALSE)
    edge_shared <- data.frame(source=character(), target=character(), interaction=character(), support_n=numeric(), evidence_class=character(), predicted_databases=character(), validated_databases=character(), validation_methods=character(), PMIDs=character(), target_class=character(), stringsAsFactors = FALSE)
  } else {
    nodes_in_shared <- unique(c(edge_shared$source, edge_shared$target))
    node_shared <- node_full %>% dplyr::filter(node %in% nodes_in_shared) %>% dplyr::arrange(node_type, dplyr::desc(degree_total), node)
  }
  safe_write_csv(edge_shared, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Edge_Lists", "Table_miRNA_Target_Edges_SharedCore.csv"))
  safe_write_csv(node_shared, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Node_Attributes", "Table_miRNA_Target_Nodes_SharedCore.csv"))
  safe_write_csv(edge_shared, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Cytoscape_Files", "Cytoscape_EdgeList_SharedCore.csv"))
  safe_write_csv(node_shared, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Cytoscape_Files", "Cytoscape_NodeAttributes_SharedCore.csv"))
  edge_summary <- data.frame(network = c("full","shared_core"), n_edges = c(nrow(edge_full), nrow(edge_shared)), n_miRNA_nodes = c(length(unique(edge_full$source)), if (nrow(edge_shared)>0) length(unique(edge_shared$source)) else 0), n_gene_nodes = c(length(unique(edge_full$target)), if (nrow(edge_shared)>0) length(unique(edge_shared$target)) else 0), stringsAsFactors = FALSE)
  node_summary <- node_full %>% dplyr::group_by(node_type) %>% dplyr::summarise(n_nodes = dplyr::n(), median_degree = stats::median(degree_total, na.rm = TRUE), max_degree = max(degree_total, na.rm = TRUE), .groups = "drop")
  safe_write_csv(edge_summary, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Table_miRNA_Target_Network_Edge_Summary.csv"))
  safe_write_csv(node_summary, file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Table_miRNA_Target_Network_Node_Summary.csv"))
  summary_tbl <- data.frame(metric = c("n_retained_miRNAs","n_full_edges","n_full_gene_nodes","n_full_total_nodes","n_shared_core_edges","n_shared_core_gene_nodes","n_shared_core_total_nodes"),
                            value = c(length(retained_mirs), nrow(edge_full), sum(node_full$node_type == "gene"), nrow(node_full), nrow(edge_shared), if (nrow(node_shared)>0) sum(node_shared$node_type == "gene") else 0, nrow(node_shared)), stringsAsFactors = FALSE)
  safe_write_csv(summary_tbl, file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Step16H_miRNA_Target_Network_Summary.csv"))
  sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16H.txt")); print(sessionInfo()); sink()
}

###############################################################################
# STEP 16I
###############################################################################
if (RUN_16I) {
  log_message("Running Step 16I")
  node_full <- read.csv(file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Node_Attributes", "Table_miRNA_Target_Nodes.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  edge_full <- read.csv(file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Edge_Lists", "Table_miRNA_Target_Edges.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  node_shared <- read.csv(file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Node_Attributes", "Table_miRNA_Target_Nodes_SharedCore.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  edge_shared <- read.csv(file.path(PROJECT_ROOT, "07_miRNA_Target_Network", "Edge_Lists", "Table_miRNA_Target_Edges_SharedCore.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  node_full$node <- clean_text(node_full$node); node_full$node_type <- clean_text(node_full$node_type); node_shared$node <- clean_text(node_shared$node); node_shared$node_type <- clean_text(node_shared$node_type)
  edge_full$source <- clean_text(edge_full$source); edge_full$target <- clean_text(edge_full$target); edge_shared$source <- clean_text(edge_shared$source); edge_shared$target <- clean_text(edge_shared$target)
  gene_nodes_full <- node_full %>% dplyr::filter(node_type == "gene", !is.na(node)) %>% dplyr::rename(gene_symbol = node) %>% dplyr::distinct(gene_symbol, .keep_all = TRUE)
  gene_symbols_full <- sort(unique(gene_nodes_full$gene_symbol))
  gene_nodes_shared <- if (nrow(node_shared) > 0) node_shared %>% dplyr::filter(node_type == "gene", !is.na(node)) %>% dplyr::rename(gene_symbol = node) %>% dplyr::distinct(gene_symbol, .keep_all = TRUE) else data.frame(gene_symbol=character(), stringsAsFactors = FALSE)
  gene_symbols_shared <- if (nrow(gene_nodes_shared)>0) sort(unique(gene_nodes_shared$gene_symbol)) else character()
  safe_write_txt_vector(gene_symbols_full, file.path(PROJECT_ROOT, "08_PPI_Network", "Node_Attributes", "PPI_Input_Genes_Full.txt"))
  safe_write_txt_vector(gene_symbols_shared, file.path(PROJECT_ROOT, "08_PPI_Network", "Node_Attributes", "PPI_Input_Genes_SharedCore.txt"))
  safe_write_csv(data.frame(gene_symbol = gene_symbols_full, stringsAsFactors = FALSE), file.path(PROJECT_ROOT, "08_PPI_Network", "Node_Attributes", "Table_PPI_Input_Genes_Full.csv"))
  safe_write_csv(data.frame(gene_symbol = gene_symbols_shared, stringsAsFactors = FALSE), file.path(PROJECT_ROOT, "08_PPI_Network", "Node_Attributes", "Table_PPI_Input_Genes_SharedCore.csv"))
  build_gene_attr <- function(edge_df, gene_symbols) {
    if (length(gene_symbols) == 0 || nrow(edge_df) == 0) return(data.frame(gene_symbol=character(), n_miRNAs_targeting=numeric(), miRNA_list=character(), target_class=character(), evidence_class=character(), support_n_max=numeric(), support_n_mean=numeric(), n_edge_records=numeric(), predicted_databases_all=character(), validated_databases_all=character(), validation_methods_all=character(), PMIDs_all=character(), stringsAsFactors = FALSE))
    edge_df <- edge_df %>% dplyr::filter(!is.na(target), target %in% gene_symbols)
    if (nrow(edge_df) == 0) return(data.frame(gene_symbol=character(), n_miRNAs_targeting=numeric(), miRNA_list=character(), target_class=character(), evidence_class=character(), support_n_max=numeric(), support_n_mean=numeric(), n_edge_records=numeric(), predicted_databases_all=character(), validated_databases_all=character(), validation_methods_all=character(), PMIDs_all=character(), stringsAsFactors = FALSE))
    edge_df %>% dplyr::group_by(target) %>% dplyr::summarise(
      n_miRNAs_targeting = dplyr::n_distinct(source),
      miRNA_list = collapse_unique_sorted(source),
      target_class = ifelse(dplyr::n_distinct(source) == 1, "unique", "shared"),
      evidence_class = collapse_unique_sorted(evidence_class),
      support_n_max = suppressWarnings(max(as.numeric(support_n), na.rm = TRUE)),
      support_n_mean = suppressWarnings(mean(as.numeric(support_n), na.rm = TRUE)),
      n_edge_records = dplyr::n(),
      predicted_databases_all = collapse_unique_sorted(predicted_databases),
      validated_databases_all = collapse_unique_sorted(validated_databases),
      validation_methods_all = collapse_unique_sorted(validation_methods),
      PMIDs_all = collapse_unique_sorted(PMIDs),
      .groups = "drop"
    ) %>% dplyr::rename(gene_symbol = target) %>% dplyr::mutate(support_n_max = ifelse(is.infinite(support_n_max), NA_real_, support_n_max), support_n_mean = ifelse(is.nan(support_n_mean), NA_real_, support_n_mean)) %>% dplyr::arrange(dplyr::desc(n_miRNAs_targeting), dplyr::desc(support_n_max), gene_symbol)
  }
  gene_attr_full <- build_gene_attr(edge_full, gene_symbols_full)
  gene_attr_shared <- build_gene_attr(edge_shared, gene_symbols_shared)
  safe_write_csv(gene_attr_full, file.path(PROJECT_ROOT, "08_PPI_Network", "Node_Attributes", "Table_PPI_GeneAttributes_Full.csv"))
  safe_write_csv(gene_attr_shared, file.path(PROJECT_ROOT, "08_PPI_Network", "Node_Attributes", "Table_PPI_GeneAttributes_SharedCore.csv"))
  prehub_full <- gene_attr_full %>% dplyr::mutate(prehub_rank_primary = dplyr::min_rank(dplyr::desc(n_miRNAs_targeting)), prehub_rank_secondary = dplyr::min_rank(dplyr::desc(support_n_max))) %>% dplyr::arrange(dplyr::desc(n_miRNAs_targeting), dplyr::desc(support_n_max), gene_symbol)
  prehub_shared <- gene_attr_shared %>% dplyr::mutate(prehub_rank_primary = dplyr::min_rank(dplyr::desc(n_miRNAs_targeting)), prehub_rank_secondary = dplyr::min_rank(dplyr::desc(support_n_max))) %>% dplyr::arrange(dplyr::desc(n_miRNAs_targeting), dplyr::desc(support_n_max), gene_symbol)
  safe_write_csv(prehub_full, file.path(PROJECT_ROOT, "08_PPI_Network", "Hub_Genes", "Table_PreHub_Candidates_Full.csv"))
  safe_write_csv(prehub_shared, file.path(PROJECT_ROOT, "08_PPI_Network", "Hub_Genes", "Table_PreHub_Candidates_SharedCore.csv"))
  safe_write_csv(utils::head(prehub_full, 20), file.path(PROJECT_ROOT, "08_PPI_Network", "Hub_Genes", "Table_PreHub_Candidates_Top20_Full.csv"))
  safe_write_csv(utils::head(prehub_shared, 20), file.path(PROJECT_ROOT, "08_PPI_Network", "Hub_Genes", "Table_PreHub_Candidates_Top20_SharedCore.csv"))
  summary_tbl <- data.frame(metric = c("n_full_gene_symbols_for_ppi","n_shared_core_gene_symbols_for_ppi","n_full_prehub_candidates","n_shared_core_prehub_candidates","n_full_top20_written","n_shared_core_top20_written"), value = c(length(gene_symbols_full), length(gene_symbols_shared), nrow(prehub_full), nrow(prehub_shared), min(20, nrow(prehub_full)), min(20, nrow(prehub_shared))), stringsAsFactors = FALSE)
  safe_write_csv(summary_tbl, file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Step16I_PPI_Input_Summary.csv"))
  sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16I.txt")); print(sessionInfo()); sink()
}

###############################################################################
# STEP 16K CLEAN FINAL
###############################################################################
if (RUN_16K) {
  log_message("Running Step 16K CLEAN FINAL")
  overlap_genes <- extract_unique_genes(file.path(PROJECT_ROOT, "05_Overlap_Targets", "Tables", "Table_All_Overlap_Targets.csv"), gene_col = "gene_symbol")
  shared_genes <- extract_unique_genes(file.path(PROJECT_ROOT, "06_Shared_Targets", "Table_Shared_Targets_2plus.csv"), gene_col = "gene_symbol")
  final_hub_file <- file.path(PROJECT_ROOT, "08_PPI_Network", "Hub_Genes", "Table_Final_Hub_Gene_Set.csv")
  prehub_shared_file <- file.path(PROJECT_ROOT, "08_PPI_Network", "Hub_Genes", "Table_PreHub_Candidates_SharedCore.csv")
  if (file.exists(final_hub_file)) {
    hub_genes <- extract_unique_genes(final_hub_file, gene_col = "gene_symbol")
    hub_source_used <- "final_hub_gene_set"
  } else {
    hub_genes <- extract_unique_genes(prehub_shared_file, gene_col = "gene_symbol", top_n = HUB_TOP_N_FALLBACK)
    hub_source_used <- paste0("prehub_sharedcore_top", HUB_TOP_N_FALLBACK)
  }
  safe_write_csv(data.frame(source_gene_set = c(rep("OverlapTargets", length(overlap_genes)), rep("SharedTargets", length(shared_genes)), rep("HubGenes", length(hub_genes))), gene_symbol = c(overlap_genes, shared_genes, hub_genes), stringsAsFactors = FALSE), file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_Enrichment_Input_GeneSets.csv"))
  run_one_enrichment <- function(symbols, source_name, orgdb) {
    if (length(symbols) == 0) {
      empty <- standardize_enrichment_table(NULL, source_name, "GO_BP")
      return(list(bp_table = empty, mf_table = standardize_enrichment_table(NULL, source_name, "GO_MF"), cc_table = standardize_enrichment_table(NULL, source_name, "GO_CC"), kegg_table = standardize_enrichment_table(NULL, source_name, "KEGG"), summary = data.frame(source_gene_set = source_name, n_input_symbols = 0, n_mapped_entrez = 0, n_GO_BP_terms = 0, n_GO_MF_terms = 0, n_GO_CC_terms = 0, n_KEGG_terms = 0, stringsAsFactors = FALSE)))
    }
    conv <- tryCatch(suppressWarnings(clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)), error = function(e) NULL)
    if (is.null(conv) || nrow(conv) == 0) {
      return(list(bp_table = standardize_enrichment_table(NULL, source_name, "GO_BP"), mf_table = standardize_enrichment_table(NULL, source_name, "GO_MF"), cc_table = standardize_enrichment_table(NULL, source_name, "GO_CC"), kegg_table = standardize_enrichment_table(NULL, source_name, "KEGG"), summary = data.frame(source_gene_set = source_name, n_input_symbols = length(symbols), n_mapped_entrez = 0, n_GO_BP_terms = 0, n_GO_MF_terms = 0, n_GO_CC_terms = 0, n_KEGG_terms = 0, stringsAsFactors = FALSE)))
    }
    entrez_ids <- unique(conv$ENTREZID)
    ego_bp <- tryCatch(clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = orgdb, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = GO_PVALUE_CUTOFF, qvalueCutoff = GO_QVALUE_CUTOFF, readable = TRUE), error = function(e) NULL)
    ego_mf <- tryCatch(clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = orgdb, keyType = "ENTREZID", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = GO_PVALUE_CUTOFF, qvalueCutoff = GO_QVALUE_CUTOFF, readable = TRUE), error = function(e) NULL)
    ego_cc <- tryCatch(clusterProfiler::enrichGO(gene = entrez_ids, OrgDb = orgdb, keyType = "ENTREZID", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = GO_PVALUE_CUTOFF, qvalueCutoff = GO_QVALUE_CUTOFF, readable = TRUE), error = function(e) NULL)
    ekegg <- tryCatch(clusterProfiler::enrichKEGG(gene = entrez_ids, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = KEGG_PVALUE_CUTOFF, qvalueCutoff = KEGG_QVALUE_CUTOFF), error = function(e) NULL)
    bp_df <- standardize_enrichment_table(tryCatch(as.data.frame(ego_bp), error = function(e) NULL), source_name, "GO_BP")
    mf_df <- standardize_enrichment_table(tryCatch(as.data.frame(ego_mf), error = function(e) NULL), source_name, "GO_MF")
    cc_df <- standardize_enrichment_table(tryCatch(as.data.frame(ego_cc), error = function(e) NULL), source_name, "GO_CC")
    kegg_df <- standardize_enrichment_table(tryCatch(as.data.frame(ekegg), error = function(e) NULL), source_name, "KEGG")
    list(bp_table = bp_df, mf_table = mf_df, cc_table = cc_df, kegg_table = kegg_df,
         summary = data.frame(source_gene_set = source_name, n_input_symbols = length(symbols), n_mapped_entrez = length(entrez_ids), n_GO_BP_terms = nrow(bp_df), n_GO_MF_terms = nrow(mf_df), n_GO_CC_terms = nrow(cc_df), n_KEGG_terms = nrow(kegg_df), stringsAsFactors = FALSE))
  }
  res_overlap <- run_one_enrichment(overlap_genes, "OverlapTargets", org.Hs.eg.db)
  res_shared <- run_one_enrichment(shared_genes, "SharedTargets", org.Hs.eg.db)
  res_hubs <- run_one_enrichment(hub_genes, "HubGenes", org.Hs.eg.db)
  # save tables
  safe_write_csv(res_overlap$bp_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_GO_BP_OverlapTargets.csv"))
  safe_write_csv(res_overlap$mf_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_GO_MF_OverlapTargets.csv"))
  safe_write_csv(res_overlap$cc_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_GO_CC_OverlapTargets.csv"))
  safe_write_csv(res_overlap$kegg_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_KEGG_OverlapTargets.csv"))
  safe_write_csv(res_shared$bp_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_GO_BP_SharedTargets.csv"))
  safe_write_csv(res_shared$mf_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_GO_MF_SharedTargets.csv"))
  safe_write_csv(res_shared$cc_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_GO_CC_SharedTargets.csv"))
  safe_write_csv(res_shared$kegg_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_KEGG_SharedTargets.csv"))
  safe_write_csv(res_hubs$bp_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_GO_BP_HubGenes.csv"))
  safe_write_csv(res_hubs$mf_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_GO_MF_HubGenes.csv"))
  safe_write_csv(res_hubs$cc_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_GO_CC_HubGenes.csv"))
  safe_write_csv(res_hubs$kegg_table, file.path(PROJECT_ROOT, "09_Enrichment", "Tables", "Table_KEGG_HubGenes.csv"))
  # save colorful plots
  save_colorful_dotplot(res_overlap$bp_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_BP", "Dotplot_GO_BP_OverlapTargets"), "GO BP - Overlap Targets", TOP_N_PLOT)
  save_colorful_barplot(res_overlap$bp_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_BP", "Barplot_GO_BP_OverlapTargets"), "GO BP - Overlap Targets", TOP_N_PLOT)
  save_colorful_dotplot(res_shared$bp_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_BP", "Dotplot_GO_BP_SharedTargets"), "GO BP - Shared Targets", TOP_N_PLOT)
  save_colorful_barplot(res_shared$bp_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_BP", "Barplot_GO_BP_SharedTargets"), "GO BP - Shared Targets", TOP_N_PLOT)
  save_colorful_dotplot(res_hubs$bp_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_BP", "Dotplot_GO_BP_HubGenes"), "GO BP - Hub Genes", TOP_N_PLOT)
  save_colorful_barplot(res_hubs$bp_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_BP", "Barplot_GO_BP_HubGenes"), "GO BP - Hub Genes", TOP_N_PLOT)
  save_colorful_dotplot(res_overlap$mf_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_MF", "Dotplot_GO_MF_OverlapTargets"), "GO MF - Overlap Targets", TOP_N_PLOT)
  save_colorful_barplot(res_overlap$mf_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_MF", "Barplot_GO_MF_OverlapTargets"), "GO MF - Overlap Targets", TOP_N_PLOT)
  save_colorful_dotplot(res_shared$mf_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_MF", "Dotplot_GO_MF_SharedTargets"), "GO MF - Shared Targets", TOP_N_PLOT)
  save_colorful_barplot(res_shared$mf_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_MF", "Barplot_GO_MF_SharedTargets"), "GO MF - Shared Targets", TOP_N_PLOT)
  save_colorful_dotplot(res_hubs$mf_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_MF", "Dotplot_GO_MF_HubGenes"), "GO MF - Hub Genes", TOP_N_PLOT)
  save_colorful_barplot(res_hubs$mf_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_MF", "Barplot_GO_MF_HubGenes"), "GO MF - Hub Genes", TOP_N_PLOT)
  save_colorful_dotplot(res_overlap$cc_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_CC", "Dotplot_GO_CC_OverlapTargets"), "GO CC - Overlap Targets", TOP_N_PLOT)
  save_colorful_barplot(res_overlap$cc_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_CC", "Barplot_GO_CC_OverlapTargets"), "GO CC - Overlap Targets", TOP_N_PLOT)
  save_colorful_dotplot(res_shared$cc_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_CC", "Dotplot_GO_CC_SharedTargets"), "GO CC - Shared Targets", TOP_N_PLOT)
  save_colorful_barplot(res_shared$cc_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_CC", "Barplot_GO_CC_SharedTargets"), "GO CC - Shared Targets", TOP_N_PLOT)
  save_colorful_dotplot(res_hubs$cc_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_CC", "Dotplot_GO_CC_HubGenes"), "GO CC - Hub Genes", TOP_N_PLOT)
  save_colorful_barplot(res_hubs$cc_table, file.path(PROJECT_ROOT, "09_Enrichment", "GO_CC", "Barplot_GO_CC_HubGenes"), "GO CC - Hub Genes", TOP_N_PLOT)
  save_colorful_dotplot(res_overlap$kegg_table, file.path(PROJECT_ROOT, "09_Enrichment", "KEGG", "Dotplot_KEGG_OverlapTargets"), "KEGG - Overlap Targets", TOP_N_PLOT)
  save_colorful_barplot(res_overlap$kegg_table, file.path(PROJECT_ROOT, "09_Enrichment", "KEGG", "Barplot_KEGG_OverlapTargets"), "KEGG - Overlap Targets", TOP_N_PLOT)
  save_colorful_dotplot(res_shared$kegg_table, file.path(PROJECT_ROOT, "09_Enrichment", "KEGG", "Dotplot_KEGG_SharedTargets"), "KEGG - Shared Targets", TOP_N_PLOT)
  save_colorful_barplot(res_shared$kegg_table, file.path(PROJECT_ROOT, "09_Enrichment", "KEGG", "Barplot_KEGG_SharedTargets"), "KEGG - Shared Targets", TOP_N_PLOT)
  save_colorful_dotplot(res_hubs$kegg_table, file.path(PROJECT_ROOT, "09_Enrichment", "KEGG", "Dotplot_KEGG_HubGenes"), "KEGG - Hub Genes", TOP_N_PLOT)
  save_colorful_barplot(res_hubs$kegg_table, file.path(PROJECT_ROOT, "09_Enrichment", "KEGG", "Barplot_KEGG_HubGenes"), "KEGG - Hub Genes", TOP_N_PLOT)
  summary_tbl <- dplyr::bind_rows(res_overlap$summary, res_shared$summary, res_hubs$summary)
  summary_tbl$hub_gene_source_used <- ifelse(summary_tbl$source_gene_set == "HubGenes", hub_source_used, NA_character_)
  safe_write_csv(summary_tbl, file.path(PROJECT_ROOT, "11_Methods_and_Logs", "Step16K_CLEAN_FINAL_Summary.csv"))
  sink(file.path(PROJECT_ROOT, "11_Methods_and_Logs", "sessionInfo_Step16K_CLEAN_FINAL.txt")); print(sessionInfo()); sink()
}

log_message("STEP16 master pipeline finished.")
cat("\n==================== STEP16 MASTER PIPELINE DONE ====================\n")
cat("Project root:\n", PROJECT_ROOT, "\n")
cat("Main note: Step 16J is not included because it requires external STRING/Cytoscape exports.\n")
cat("===================================================================\n")
