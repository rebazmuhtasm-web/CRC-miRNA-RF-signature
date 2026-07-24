###############################################################################
# Reviewer-requested uncertainty, calibration, and permutation analyses
# FIXED SELF-CONTAINED VERSION 2
#
# This version fixes the previous error:
#   "The locked Random Forest model was not found"
#
# What is new in this fixed version:
#   - It does NOT require RBZ_step10_final_locked_modeling.rds to already exist.
#   - If the final RF object is missing, it reconstructs the original Train80/Test20
#     split from RBZ_step3_large_preprocessed.rds using the same seed and settings
#     used in the master script, then retrains the locked RF model using the locked
#     miRNAs from RBZ_step7_final_locked_signature.rds.
#   - It saves a reconstructed RF object so the reviewer-requested analyses can run.
#
# Recommended use:
#   1) Restart R/RStudio.
#   2) Run:
#      setwd("D:/RUN_miRDeep2_edgeR_QC_20260318_180642/mixed/Raw data/RF ONLY/RBZ_STEP1_INPUT_DISCOVERY_20260405_003306")
#      source("Reviewer_Requested_CI_Calibration_Permutation_FIXED_V2.R")
#
# Outputs are written to:
#   13_Reviewer_Requested_Uncertainty
###############################################################################

options(stringsAsFactors = FALSE)
set.seed(20260701)

###############################################################################
# 0) USER SETTINGS
###############################################################################

# Leave NULL if you run the script from inside the project folder.
# If auto-detection fails, paste your project path here:
# PROJECT_DIR_MANUAL <- "D:/RUN_miRDeep2_edgeR_QC_20260318_180642/mixed/Raw data/RF ONLY/RBZ_STEP1_INPUT_DISCOVERY_20260405_003306"
PROJECT_DIR_MANUAL <- NULL

# These must match the master workflow.
TRAIN_PROP <- 0.80
SPLIT_SEED <- 20260405
RF_NTREE <- 1000
DECISION_THRESHOLD <- 0.5

# Reviewer-requested uncertainty settings.
N_PERMUTATIONS <- 1000
BOOT_N_AUC_CI  <- 2000
N_CALIBRATION_BINS <- 5

# External validation will be attempted if the raw external count/meta files exist.
RUN_EXTERNAL_VALIDATION <- TRUE

###############################################################################
# 1) PACKAGE SETUP
###############################################################################

required_pkgs <- c("pROC", "dplyr", "ggplot2", "randomForest", "edgeR", "data.table", "caret")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Please install these R packages first: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(ggplot2)
  library(randomForest)
  library(edgeR)
  library(data.table)
  library(caret)
})

###############################################################################
# 2) PROJECT DIRECTORY AND OUTPUT FOLDERS
###############################################################################

if (!is.null(PROJECT_DIR_MANUAL) && nzchar(PROJECT_DIR_MANUAL)) {
  PROJECT_DIR <- normalizePath(PROJECT_DIR_MANUAL, winslash = "/", mustWork = TRUE)
} else {
  PROJECT_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT_DIR <- dirname(PROJECT_DIR)

DIR_04_OBJ <- file.path(PROJECT_DIR, "04_Objects")
DIR_05_LARGE <- file.path(PROJECT_DIR, "05_Large_Preprocessing")
DIR_06_SPLIT <- file.path(PROJECT_DIR, "06_Train_Test_Split")
DIR_10_MOD <- file.path(PROJECT_DIR, "10_Final_Locked_Modeling")
DIR_11_SUPP <- file.path(PROJECT_DIR, "11_Supportive_FullLarge_Locked_RF")
DIR_12_EXT <- file.path(PROJECT_DIR, "12_External_HTS_Validation")
DIR_13_REVIEWER <- file.path(PROJECT_DIR, "13_Reviewer_Requested_Uncertainty")

dir.create(DIR_04_OBJ, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_06_SPLIT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_10_MOD, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_11_SUPP, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_12_EXT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_13_REVIEWER, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(DIR_13_REVIEWER, "Run_Log_Reviewer_Requested_Uncertainty.txt")

log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

log_message("Reviewer-requested fixed V2 script started.")
log_message("Using PROJECT_DIR: ", PROJECT_DIR)
log_message("Using ROOT_DIR: ", ROOT_DIR)

###############################################################################
# 3) HELPER FUNCTIONS
###############################################################################

write_csv_safe <- function(df, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(as.data.frame(df), file = file, row.names = row.names, na = "")
  invisible(file)
}

sanitize_mir <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.", "-", x)
  x <- gsub("_", "-", x)
  x <- gsub("--+", "-", x)
  x <- sub("^hsa-mir", "hsa-miR", x, ignore.case = TRUE)
  x
}

sanitize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  gsub("[[:space:]]+", "", x)
}

normalize_group <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  ifelse(grepl("tumou?r|cancer|crc|case", x, ignore.case = TRUE), "Tumor",
         ifelse(grepl("normal|control|adjacent", x, ignore.case = TRUE), "Normal", x))
}

detect_sample_column <- function(meta, expr_ids) {
  expr_ids <- sanitize_id(expr_ids)
  scores <- sapply(seq_along(meta), function(i) sum(sanitize_id(meta[[i]]) %in% expr_ids, na.rm = TRUE))
  idx <- which.max(scores)
  if (length(idx) == 0 || scores[idx] == 0) return(colnames(meta)[1])
  colnames(meta)[idx]
}

detect_group_column <- function(meta) {
  cand <- c("Group", "group", "Class", "class", "Condition", "condition", "Type", "type", "status", "Status")
  hit <- cand[cand %in% colnames(meta)][1]
  if (!is.na(hit)) return(hit)
  for (nm in colnames(meta)) {
    vals <- unique(tolower(as.character(meta[[nm]])))
    if (any(grepl("tumou?r|cancer|crc|case", vals)) &&
        any(grepl("normal|control|adjacent", vals))) {
      return(nm)
    }
  }
  NA_character_
}

detect_pair_column <- function(meta) {
  cand <- c("pair_id", "PairID", "pair", "pairID", "pair.id", "Pair", "patient_id", "PatientID")
  hit <- cand[cand %in% colnames(meta)][1]
  if (!is.na(hit)) return(hit)
  NA_character_
}

read_count_dataset <- function(count_file, meta_file, dataset_name) {
  if (!file.exists(count_file)) stop("Count file not found for ", dataset_name, ": ", count_file)
  if (!file.exists(meta_file)) stop("Metadata file not found for ", dataset_name, ": ", meta_file)

  counts <- read.csv(count_file, stringsAsFactors = FALSE, check.names = FALSE)
  meta   <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)

  feature_col <- colnames(counts)[1]
  mirs <- sanitize_mir(counts[[feature_col]])
  cnt_df <- counts[, -1, drop = FALSE]
  cnt_mat <- as.matrix(data.frame(lapply(cnt_df, function(x) as.numeric(as.character(x))), check.names = FALSE))
  rownames(cnt_mat) <- mirs
  colnames(cnt_mat) <- sanitize_id(colnames(cnt_df))

  sample_col <- detect_sample_column(meta, colnames(cnt_mat))
  group_col  <- detect_group_column(meta)
  pair_col   <- detect_pair_column(meta)

  meta$SampleID <- sanitize_id(meta[[sample_col]])
  if (is.na(group_col)) stop("Could not detect group column in metadata for ", dataset_name)
  meta$Group <- normalize_group(meta[[group_col]])
  if (!is.na(pair_col)) meta$pair_id <- as.character(meta[[pair_col]]) else meta$pair_id <- NA_character_

  common_ids <- intersect(colnames(cnt_mat), meta$SampleID)
  if (length(common_ids) == 0) stop("No shared sample IDs for ", dataset_name)

  cnt_mat <- cnt_mat[, common_ids, drop = FALSE]
  meta <- meta[match(common_ids, meta$SampleID), , drop = FALSE]
  keep <- meta$Group %in% c("Normal", "Tumor")
  meta <- meta[keep, , drop = FALSE]
  cnt_mat <- cnt_mat[, meta$SampleID, drop = FALSE]

  list(dataset_name = dataset_name, counts = cnt_mat, info = meta)
}

calc_tmm_logcpm <- function(counts_mat) {
  dge <- edgeR::DGEList(counts = counts_mat)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR::cpm(dge, log = TRUE, prior.count = 1)
}

read_matrix_csv <- function(file) {
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  rn <- sanitize_mir(df[[1]])
  mat <- as.matrix(data.frame(lapply(df[, -1, drop = FALSE], function(x) as.numeric(as.character(x))), check.names = FALSE))
  rownames(mat) <- rn
  colnames(mat) <- sanitize_id(colnames(df)[-1])
  mat
}

extract_final_signature <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    if (!is.null(x$final_signature)) return(x$final_signature)
    if (!is.null(x$signature)) return(x$signature)
    if (!is.null(x$locked_signature)) return(x$locked_signature)
  }
  NULL
}

extract_large_after <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    if (!is.null(x$large_after)) return(x$large_after)
    if (!is.null(x$expr_native) && !is.null(x$info)) return(x)
  }
  NULL
}

extract_train_test <- function(x) {
  if (is.null(x)) return(list(train = NULL, test = NULL))
  if (is.list(x)) {
    tr <- NULL
    te <- NULL
    if (!is.null(x$train)) tr <- x$train
    if (!is.null(x$test)) te <- x$test
    if (!is.null(x$train_ds)) tr <- x$train_ds
    if (!is.null(x$test_ds)) te <- x$test_ds
    return(list(train = tr, test = te))
  }
  list(train = NULL, test = NULL)
}

calc_metrics_basic <- function(obs, prob, threshold = 0.5) {
  obs <- factor(obs, levels = c("Normal", "Tumor"))
  pred <- factor(ifelse(prob >= threshold, "Tumor", "Normal"), levels = c("Normal", "Tumor"))

  tp <- sum(obs == "Tumor" & pred == "Tumor", na.rm = TRUE)
  tn <- sum(obs == "Normal" & pred == "Normal", na.rm = TRUE)
  fp <- sum(obs == "Normal" & pred == "Tumor", na.rm = TRUE)
  fn <- sum(obs == "Tumor" & pred == "Normal", na.rm = TRUE)

  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  precision   <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  accuracy    <- mean(obs == pred, na.rm = TRUE)
  f1 <- ifelse(is.na(precision) || is.na(sensitivity) || (precision + sensitivity) == 0,
               NA_real_, 2 * precision * sensitivity / (precision + sensitivity))

  roc_obj <- tryCatch(
    pROC::roc(response = obs, predictor = prob, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  auc_val <- if (is.null(roc_obj)) NA_real_ else as.numeric(pROC::auc(roc_obj))

  data.frame(
    N = length(obs),
    Normal_n = sum(obs == "Normal"),
    Tumor_n = sum(obs == "Tumor"),
    AUC = auc_val,
    TP = tp, TN = tn, FP = fp, FN = fn,
    Accuracy = accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    Brier = mean((as.numeric(obs == "Tumor") - prob)^2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

binom_ci_row <- function(x, n) {
  if (is.na(x) || is.na(n) || n <= 0) return(c(NA_real_, NA_real_, NA_real_))
  bt <- binom.test(x, n)
  c(as.numeric(x / n), as.numeric(bt$conf.int[1]), as.numeric(bt$conf.int[2]))
}

auc_ci <- function(obs, prob, boot_n = 2000) {
  obs <- factor(obs, levels = c("Normal", "Tumor"))
  roc_obj <- tryCatch(
    pROC::roc(obs, prob, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(roc_obj)) return(c(AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_, Method = NA_character_))

  auc_val <- as.numeric(pROC::auc(roc_obj))

  ci_delong <- tryCatch(
    suppressWarnings(pROC::ci.auc(roc_obj, method = "delong")),
    error = function(e) NULL
  )
  if (!is.null(ci_delong) && all(is.finite(as.numeric(ci_delong)))) {
    return(c(AUC = auc_val, CI_low = as.numeric(ci_delong[1]), CI_high = as.numeric(ci_delong[3]), Method = "DeLong"))
  }

  ci_boot <- tryCatch(
    suppressWarnings(pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = boot_n)),
    error = function(e) NULL
  )
  if (!is.null(ci_boot) && all(is.finite(as.numeric(ci_boot)))) {
    return(c(AUC = auc_val, CI_low = as.numeric(ci_boot[1]), CI_high = as.numeric(ci_boot[3]), Method = "Bootstrap"))
  }

  c(AUC = auc_val, CI_low = NA_real_, CI_high = NA_real_, Method = "CI_not_estimable")
}

make_predict_table <- function(dataset_name, obs, prob, sample_id = NULL) {
  if (is.null(sample_id)) sample_id <- paste0(dataset_name, "_", seq_along(obs))
  pred <- ifelse(prob >= DECISION_THRESHOLD, "Tumor", "Normal")
  data.frame(
    Dataset = dataset_name,
    SampleID = sample_id,
    Observed = as.character(obs),
    Probability_Tumor = as.numeric(prob),
    Predicted = pred,
    stringsAsFactors = FALSE
  )
}

evaluate_prediction_table <- function(pred_df) {
  out <- list()

  for (ds in unique(pred_df$Dataset)) {
    sub <- pred_df[pred_df$Dataset == ds, , drop = FALSE]
    obs <- factor(sub$Observed, levels = c("Normal", "Tumor"))
    prob <- as.numeric(sub$Probability_Tumor)

    m <- calc_metrics_basic(obs, prob, DECISION_THRESHOLD)
    auc_row <- auc_ci(obs, prob, BOOT_N_AUC_CI)

    out[[ds]] <- list(metrics = m, auc_ci = auc_row)
  }

  out
}

calibration_bins <- function(pred_df, n_bins = 5) {
  all_bins <- list()
  brier_rows <- list()

  for (ds in unique(pred_df$Dataset)) {
    sub <- pred_df[pred_df$Dataset == ds, , drop = FALSE]
    sub$y <- as.numeric(sub$Observed == "Tumor")
    sub$p <- as.numeric(sub$Probability_Tumor)

    brier_rows[[ds]] <- data.frame(
      Dataset = ds,
      N = nrow(sub),
      Brier = mean((sub$y - sub$p)^2, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    if (nrow(sub) < 4 || length(unique(sub$p)) < 2) {
      sub$bin <- 1
    } else {
      br <- unique(stats::quantile(sub$p, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
      if (length(br) < 3) {
        sub$bin <- 1
      } else {
        sub$bin <- cut(sub$p, breaks = br, include.lowest = TRUE, labels = FALSE)
      }
    }

    tmp <- sub |>
      dplyr::group_by(bin) |>
      dplyr::summarise(
        Dataset = ds,
        N = dplyr::n(),
        Mean_Predicted_Probability = mean(p, na.rm = TRUE),
        Observed_Tumor_Fraction = mean(y, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::select(Dataset, bin, N, Mean_Predicted_Probability, Observed_Tumor_Fraction)

    all_bins[[ds]] <- as.data.frame(tmp)
  }

  list(
    bins = dplyr::bind_rows(all_bins),
    brier = dplyr::bind_rows(brier_rows)
  )
}

plot_calibration <- function(calib_df, ds, file) {
  sub <- calib_df[calib_df$Dataset == ds, , drop = FALSE]
  if (nrow(sub) == 0) return(invisible(NULL))
  p <- ggplot(sub, aes(x = Mean_Predicted_Probability, y = Observed_Tumor_Fraction)) +
    geom_point(size = 2) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = paste0("Calibration plot: ", ds),
      x = "Mean predicted tumour probability",
      y = "Observed tumour fraction"
    ) +
    theme_bw(base_size = 12)
  ggsave(file, plot = p, width = 6, height = 5, dpi = 300)
}

plot_roc <- function(pred_df, ds, file) {
  sub <- pred_df[pred_df$Dataset == ds, , drop = FALSE]
  obs <- factor(sub$Observed, levels = c("Normal", "Tumor"))
  prob <- as.numeric(sub$Probability_Tumor)
  roc_obj <- tryCatch(pROC::roc(obs, prob, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE),
                      error = function(e) NULL)
  if (is.null(roc_obj)) return(invisible(NULL))
  png(file, width = 1800, height = 1500, res = 300)
  plot(roc_obj, main = paste0("ROC: ", ds, " (AUC = ", sprintf("%.3f", as.numeric(pROC::auc(roc_obj))), ")"))
  abline(a = 0, b = 1, lty = 2)
  dev.off()
}

permutation_auc <- function(obs, prob, n_perm = 1000, dataset_name = "Dataset") {
  obs <- factor(obs, levels = c("Normal", "Tumor"))
  prob <- as.numeric(prob)

  roc_obs <- tryCatch(pROC::roc(obs, prob, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE),
                      error = function(e) NULL)
  observed_auc <- if (is.null(roc_obs)) NA_real_ else as.numeric(pROC::auc(roc_obs))

  perm_auc <- rep(NA_real_, n_perm)
  for (i in seq_len(n_perm)) {
    obs_perm <- sample(obs)
    roc_i <- tryCatch(pROC::roc(obs_perm, prob, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE),
                      error = function(e) NULL)
    perm_auc[i] <- if (is.null(roc_i)) NA_real_ else as.numeric(pROC::auc(roc_i))
  }

  p_emp <- ifelse(is.na(observed_auc), NA_real_,
                  (sum(perm_auc >= observed_auc, na.rm = TRUE) + 1) / (sum(!is.na(perm_auc)) + 1))

  data.frame(
    Dataset = dataset_name,
    Observed_AUC = observed_auc,
    N_Permutations = n_perm,
    Empirical_P = p_emp,
    Permuted_AUC_Mean = mean(perm_auc, na.rm = TRUE),
    Permuted_AUC_SD = sd(perm_auc, na.rm = TRUE),
    Permuted_AUC_2.5 = as.numeric(quantile(perm_auc, 0.025, na.rm = TRUE)),
    Permuted_AUC_97.5 = as.numeric(quantile(perm_auc, 0.975, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

plot_permutation <- function(obs_row, perm_values, ds, file) {
  png(file, width = 1800, height = 1500, res = 300)
  hist(perm_values, breaks = 30, main = paste0("Permutation AUC null distribution: ", ds),
       xlab = "Permuted-label AUC")
  abline(v = obs_row$Observed_AUC, lwd = 3)
  legend("topright", legend = paste0("Observed AUC = ", sprintf("%.3f", obs_row$Observed_AUC)), bty = "n")
  dev.off()
}

###############################################################################
# 4) LOAD AVAILABLE PROJECT OBJECTS
###############################################################################

load_if_exists <- function(file) {
  if (file.exists(file)) {
    log_message("Loading: ", file)
    return(readRDS(file))
  }
  NULL
}

# Large preprocessed object.
large_after <- NULL
large_rds_candidates <- c(
  file.path(DIR_04_OBJ, "RBZ_step3_large_preprocessed.rds"),
  file.path(DIR_04_OBJ, "RBZ_master_steps01_15_pipeline_object.rds")
)

for (f in large_rds_candidates) {
  z <- load_if_exists(f)
  large_after <- extract_large_after(z)
  if (!is.null(large_after)) break
}

if (is.null(large_after)) {
  large_mat <- file.path(DIR_05_LARGE, "Matrix_Large_TMM_logCPM.csv")
  large_meta <- file.path(DIR_05_LARGE, "sample_info_large_preprocessed.csv")
  large_counts <- file.path(DIR_05_LARGE, "Matrix_Large_Filtered_Counts.csv")
  if (file.exists(large_mat) && file.exists(large_meta)) {
    log_message("Loading large matrix/meta CSV fallback.")
    large_after <- list(
      expr_native = read_matrix_csv(large_mat),
      info = read.csv(large_meta, stringsAsFactors = FALSE, check.names = FALSE)
    )
    if (file.exists(large_counts)) large_after$counts <- read_matrix_csv(large_counts)
  }
}

if (is.null(large_after) || is.null(large_after$expr_native) || is.null(large_after$info)) {
  stop("Could not load large_after with expr_native and info. Check 04_Objects/RBZ_step3_large_preprocessed.rds.")
}

# Ensure row/sample names are clean.
rownames(large_after$expr_native) <- sanitize_mir(rownames(large_after$expr_native))
colnames(large_after$expr_native) <- sanitize_id(colnames(large_after$expr_native))
if (!is.null(large_after$counts)) {
  rownames(large_after$counts) <- sanitize_mir(rownames(large_after$counts))
  colnames(large_after$counts) <- sanitize_id(colnames(large_after$counts))
}

# Final locked signature.
final_signature <- NULL
sig_rds_candidates <- c(
  file.path(DIR_04_OBJ, "RBZ_step7_final_locked_signature.rds"),
  file.path(DIR_04_OBJ, "RBZ_master_steps01_15_pipeline_object.rds")
)

for (f in sig_rds_candidates) {
  z <- load_if_exists(f)
  final_signature <- extract_final_signature(z)
  if (!is.null(final_signature)) break
}

if (is.null(final_signature)) {
  sig_csv <- file.path(PROJECT_DIR, "09_Final_Locked_Signature", "Table_Final_Locked_Signature.csv")
  if (file.exists(sig_csv)) {
    log_message("Loading locked signature CSV fallback: ", sig_csv)
    final_signature <- read.csv(sig_csv, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

if (is.null(final_signature)) {
  warning("Could not load locked signature. Using hard-coded manuscript locked miRNAs.")
  locked_mirs <- c("hsa-miR-549a-5p", "hsa-miR-135b-5p", "hsa-miR-378a-5p")
} else {
  if (!"miRNA" %in% colnames(final_signature)) {
    colnames(final_signature)[1] <- "miRNA"
  }
  locked_mirs <- sanitize_mir(final_signature$miRNA)
  locked_mirs <- locked_mirs[!is.na(locked_mirs) & nzchar(locked_mirs)]
}

if (length(locked_mirs) == 0) {
  stop("Locked miRNA list is empty.")
}

missing_locked <- setdiff(locked_mirs, rownames(large_after$expr_native))
if (length(missing_locked) > 0) {
  stop("These locked miRNAs were not found in large_after$expr_native: ", paste(missing_locked, collapse = ", "))
}

log_message("Locked miRNAs: ", paste(locked_mirs, collapse = ", "))

###############################################################################
# 5) LOAD OR RECONSTRUCT TRAIN/TEST SPLIT
###############################################################################

train_ds <- NULL
test_ds <- NULL

split_rds <- file.path(DIR_04_OBJ, "RBZ_step4_train_test_split.rds")
z <- load_if_exists(split_rds)
tt <- extract_train_test(z)
train_ds <- tt$train
test_ds <- tt$test

if (is.null(train_ds) || is.null(test_ds)) {
  train_mat <- file.path(DIR_06_SPLIT, "Matrix_Train80_TMM_logCPM.csv")
  train_meta <- file.path(DIR_06_SPLIT, "sample_info_train80.csv")
  test_mat  <- file.path(DIR_06_SPLIT, "Matrix_Test20_TMM_logCPM.csv")
  test_meta <- file.path(DIR_06_SPLIT, "sample_info_test20.csv")

  if (file.exists(train_mat) && file.exists(train_meta) && file.exists(test_mat) && file.exists(test_meta)) {
    log_message("Loading train/test CSV fallback.")
    train_ds <- list(
      expr_native = read_matrix_csv(train_mat),
      info = read.csv(train_meta, stringsAsFactors = FALSE, check.names = FALSE),
      dataset_name = "Large_Train80"
    )
    test_ds <- list(
      expr_native = read_matrix_csv(test_mat),
      info = read.csv(test_meta, stringsAsFactors = FALSE, check.names = FALSE),
      dataset_name = "Large_Test20"
    )
  }
}

if (is.null(train_ds) || is.null(test_ds)) {
  log_message("Train/test object was not found. Reconstructing pair-aware 80:20 split from large_after using seed ", SPLIT_SEED, ".")

  meta_large <- large_after$info
  if (!"Group" %in% colnames(meta_large)) {
    group_col <- detect_group_column(meta_large)
    if (is.na(group_col)) stop("Could not detect Group column in large_after$info.")
    meta_large$Group <- normalize_group(meta_large[[group_col]])
  }

  if (!"pair_id" %in% colnames(meta_large)) {
    pair_col <- detect_pair_column(meta_large)
    if (!is.na(pair_col)) {
      meta_large$pair_id <- as.character(meta_large[[pair_col]])
    }
  }

  if (!"pair_id" %in% colnames(meta_large) || any(is.na(meta_large$pair_id))) {
    stop("Cannot reconstruct paired split because pair_id is missing in large_after$info.")
  }

  # Match sample order to matrix columns as much as possible.
  if ("SampleID" %in% colnames(meta_large)) {
    meta_large$SampleID <- sanitize_id(meta_large$SampleID)
    common <- intersect(colnames(large_after$expr_native), meta_large$SampleID)
    if (length(common) > 0) {
      large_after$expr_native <- large_after$expr_native[, common, drop = FALSE]
      if (!is.null(large_after$counts)) large_after$counts <- large_after$counts[, common, drop = FALSE]
      meta_large <- meta_large[match(common, meta_large$SampleID), , drop = FALSE]
    }
  }

  pair_df <- meta_large |>
    dplyr::distinct(pair_id, .keep_all = TRUE)

  set.seed(SPLIT_SEED)
  train_pairs_idx <- caret::createDataPartition(
    factor(pair_df$Group, levels = c("Normal", "Tumor")),
    p = TRAIN_PROP,
    list = FALSE
  )
  train_pair_ids <- pair_df$pair_id[train_pairs_idx]
  test_pair_ids <- setdiff(pair_df$pair_id, train_pair_ids)

  train_idx <- which(meta_large$pair_id %in% train_pair_ids)
  test_idx <- which(meta_large$pair_id %in% test_pair_ids)

  train_ds <- list(
    expr_native = large_after$expr_native[, train_idx, drop = FALSE],
    info = meta_large[train_idx, , drop = FALSE],
    dataset_name = "Large_Train80"
  )
  test_ds <- list(
    expr_native = large_after$expr_native[, test_idx, drop = FALSE],
    info = meta_large[test_idx, , drop = FALSE],
    dataset_name = "Large_Test20"
  )

  if (!is.null(large_after$counts)) {
    train_ds$counts <- large_after$counts[, train_idx, drop = FALSE]
    test_ds$counts <- large_after$counts[, test_idx, drop = FALSE]
  }

  write_csv_safe(train_ds$info, file.path(DIR_06_SPLIT, "sample_info_train80_RECONSTRUCTED.csv"))
  write_csv_safe(test_ds$info, file.path(DIR_06_SPLIT, "sample_info_test20_RECONSTRUCTED.csv"))
  write_csv_safe(train_ds$expr_native, file.path(DIR_06_SPLIT, "Matrix_Train80_TMM_logCPM_RECONSTRUCTED.csv"), row.names = TRUE)
  write_csv_safe(test_ds$expr_native, file.path(DIR_06_SPLIT, "Matrix_Test20_TMM_logCPM_RECONSTRUCTED.csv"), row.names = TRUE)

  saveRDS(list(train = train_ds, test = test_ds),
          file = file.path(DIR_04_OBJ, "RBZ_step4_train_test_split_RECONSTRUCTED_FOR_REVIEWER.rds"))
}

# Clean train/test names.
rownames(train_ds$expr_native) <- sanitize_mir(rownames(train_ds$expr_native))
rownames(test_ds$expr_native) <- sanitize_mir(rownames(test_ds$expr_native))
colnames(train_ds$expr_native) <- sanitize_id(colnames(train_ds$expr_native))
colnames(test_ds$expr_native) <- sanitize_id(colnames(test_ds$expr_native))

log_message("Train samples: ", ncol(train_ds$expr_native), "; Test samples: ", ncol(test_ds$expr_native))

###############################################################################
# 6) LOAD OR RECONSTRUCT LOCKED RANDOM FOREST MODEL
###############################################################################

rf_fit <- NULL

model_rds_candidates <- c(
  file.path(DIR_04_OBJ, "RBZ_step10_final_locked_modeling.rds"),
  file.path(DIR_04_OBJ, "RBZ_master_steps01_15_pipeline_object.rds"),
  file.path(DIR_04_OBJ, "RBZ_step10_final_locked_modeling_RECONSTRUCTED_FOR_REVIEWER.rds")
)

for (f in model_rds_candidates) {
  z <- load_if_exists(f)
  if (is.list(z)) {
    if (!is.null(z$rf_fit)) rf_fit <- z$rf_fit
    if (is.null(rf_fit) && !is.null(z$rf_model)) rf_fit <- z$rf_model
  }
  if (!is.null(rf_fit)) break
}

x_train <- t(train_ds$expr_native[locked_mirs, , drop = FALSE])
x_test  <- t(test_ds$expr_native[locked_mirs, , drop = FALSE])
y_train <- factor(train_ds$info$Group, levels = c("Normal", "Tumor"))
y_test  <- factor(test_ds$info$Group, levels = c("Normal", "Tumor"))

if (is.null(rf_fit)) {
  log_message("Locked RF model object was not found. Reconstructing RF model from Train80 using ntree = ", RF_NTREE, ".")
  set.seed(SPLIT_SEED)
  rf_fit <- randomForest::randomForest(
    x = x_train,
    y = y_train,
    ntree = RF_NTREE,
    importance = TRUE
  )

  saveRDS(
    list(
      rf_fit = rf_fit,
      locked_mirs = locked_mirs,
      train_ds = train_ds,
      test_ds = test_ds,
      reconstructed = TRUE,
      reconstruction_note = "RF reconstructed for reviewer-requested uncertainty analysis because RBZ_step10_final_locked_modeling.rds was missing.",
      split_seed = SPLIT_SEED,
      ntree = RF_NTREE,
      decision_threshold = DECISION_THRESHOLD
    ),
    file = file.path(DIR_04_OBJ, "RBZ_step10_final_locked_modeling_RECONSTRUCTED_FOR_REVIEWER.rds")
  )
} else {
  log_message("Locked RF model loaded successfully.")
}

###############################################################################
# 7) PREDICTIONS FOR DISCOVERY PARTITIONS AND FULL LARGE
###############################################################################

pred_tables <- list()

train_prob <- predict(rf_fit, newdata = x_train, type = "prob")[, "Tumor"]
test_prob  <- predict(rf_fit, newdata = x_test, type = "prob")[, "Tumor"]

pred_tables[["Train80"]] <- make_predict_table(
  "Train80", y_train, train_prob,
  sample_id = if ("SampleID" %in% colnames(train_ds$info)) train_ds$info$SampleID else colnames(train_ds$expr_native)
)

pred_tables[["Test20"]] <- make_predict_table(
  "Test20", y_test, test_prob,
  sample_id = if ("SampleID" %in% colnames(test_ds$info)) test_ds$info$SampleID else colnames(test_ds$expr_native)
)

# Full large supportive prediction.
x_full <- t(large_after$expr_native[locked_mirs, , drop = FALSE])
if (!"Group" %in% colnames(large_after$info)) {
  group_col <- detect_group_column(large_after$info)
  if (!is.na(group_col)) large_after$info$Group <- normalize_group(large_after$info[[group_col]])
}
y_full <- factor(large_after$info$Group, levels = c("Normal", "Tumor"))
full_prob <- predict(rf_fit, newdata = x_full, type = "prob")[, "Tumor"]

pred_tables[["FullLarge_supportive"]] <- make_predict_table(
  "FullLarge_supportive", y_full, full_prob,
  sample_id = if ("SampleID" %in% colnames(large_after$info)) large_after$info$SampleID else colnames(large_after$expr_native)
)

###############################################################################
# 8) EXTERNAL VALIDATION IF RAW FILES EXIST
###############################################################################

if (RUN_EXTERNAL_VALIDATION) {
  ext_files <- list(
    list(id = "GSE125961",
         count = file.path(ROOT_DIR, "GSE125961_count_matrix_raw.csv"),
         meta  = file.path(ROOT_DIR, "sample_info_internal_GSE125961.csv")),
    list(id = "GSE183437",
         count = file.path(ROOT_DIR, "GSE183437_count_matrix_raw.csv"),
         meta  = file.path(ROOT_DIR, "sample_info_internal_GSE183437.csv")),
    list(id = "GSE236991",
         count = file.path(ROOT_DIR, "GSE236991_count_matrix_raw.csv"),
         meta  = file.path(ROOT_DIR, "sample_info_internal_GSE236991.csv"))
  )

  external_success <- list()

  for (e in ext_files) {
    if (file.exists(e$count) && file.exists(e$meta)) {
      log_message("Running external validation: ", e$id)
      ext_ds <- tryCatch(read_count_dataset(e$count, e$meta, e$id), error = function(err) {
        log_message("External validation failed for ", e$id, ": ", conditionMessage(err))
        NULL
      })

      if (!is.null(ext_ds)) {
        missing_ext <- setdiff(locked_mirs, rownames(ext_ds$counts))
        if (length(missing_ext) > 0) {
          log_message("Skipping ", e$id, " because locked miRNAs are missing: ", paste(missing_ext, collapse = ", "))
        } else {
          ext_expr <- calc_tmm_logcpm(ext_ds$counts)
          rownames(ext_expr) <- sanitize_mir(rownames(ext_expr))
          x_ext <- t(ext_expr[locked_mirs, , drop = FALSE])
          y_ext <- factor(ext_ds$info$Group, levels = c("Normal", "Tumor"))
          ext_prob <- predict(rf_fit, newdata = x_ext, type = "prob")[, "Tumor"]

          pred_tables[[e$id]] <- make_predict_table(
            e$id, y_ext, ext_prob,
            sample_id = if ("SampleID" %in% colnames(ext_ds$info)) ext_ds$info$SampleID else colnames(ext_expr)
          )

          external_success[[e$id]] <- pred_tables[[e$id]]
        }
      }
    } else {
      log_message("External files not found for ", e$id, ". Expected: ", e$count, " and ", e$meta)
    }
  }

  if (length(external_success) > 0) {
    pooled_ext <- dplyr::bind_rows(external_success)
    pooled_ext$Dataset <- "Pooled_external_HTS"
    pred_tables[["Pooled_external_HTS"]] <- pooled_ext
  }
}

all_pred <- dplyr::bind_rows(pred_tables)
write_csv_safe(all_pred, file.path(DIR_13_REVIEWER, "All_Locked_RF_Predictions.csv"))

###############################################################################
# 9) AUC CI, CLASSIFICATION METRIC CI, CALIBRATION, PERMUTATION
###############################################################################

evals <- evaluate_prediction_table(all_pred)

auc_rows <- list()
metric_ci_rows <- list()

for (ds in names(evals)) {
  sub <- all_pred[all_pred$Dataset == ds, , drop = FALSE]
  obs <- factor(sub$Observed, levels = c("Normal", "Tumor"))
  prob <- as.numeric(sub$Probability_Tumor)

  m <- evals[[ds]]$metrics
  auc_ci_vec <- evals[[ds]]$auc_ci

  auc_rows[[ds]] <- data.frame(
    Dataset = ds,
    N = m$N,
    Normal_n = m$Normal_n,
    Tumor_n = m$Tumor_n,
    AUC = as.numeric(auc_ci_vec["AUC"]),
    AUC_CI_low = as.numeric(auc_ci_vec["CI_low"]),
    AUC_CI_high = as.numeric(auc_ci_vec["CI_high"]),
    AUC_CI_method = as.character(auc_ci_vec["Method"]),
    stringsAsFactors = FALSE
  )

  # Exact binomial CIs for classification metrics.
  acc <- binom_ci_row(m$TP + m$TN, m$N)
  sens <- binom_ci_row(m$TP, m$TP + m$FN)
  spec <- binom_ci_row(m$TN, m$TN + m$FP)
  prec <- binom_ci_row(m$TP, m$TP + m$FP)

  metric_ci_rows[[ds]] <- data.frame(
    Dataset = ds,
    N = m$N,
    TP = m$TP, TN = m$TN, FP = m$FP, FN = m$FN,
    Accuracy = acc[1], Accuracy_CI_low = acc[2], Accuracy_CI_high = acc[3],
    Sensitivity = sens[1], Sensitivity_CI_low = sens[2], Sensitivity_CI_high = sens[3],
    Specificity = spec[1], Specificity_CI_low = spec[2], Specificity_CI_high = spec[3],
    Precision = prec[1], Precision_CI_low = prec[2], Precision_CI_high = prec[3],
    F1 = m$F1,
    Brier = m$Brier,
    Threshold = DECISION_THRESHOLD,
    stringsAsFactors = FALSE
  )
}

auc_tbl <- dplyr::bind_rows(auc_rows)
metric_ci_tbl <- dplyr::bind_rows(metric_ci_rows)

write_csv_safe(auc_tbl, file.path(DIR_13_REVIEWER, "Table_AUC_95CI.csv"))
write_csv_safe(metric_ci_tbl, file.path(DIR_13_REVIEWER, "Table_Classification_Metrics_Exact95CI.csv"))

# Calibration
cal <- calibration_bins(all_pred, n_bins = N_CALIBRATION_BINS)
write_csv_safe(cal$bins, file.path(DIR_13_REVIEWER, "Table_Calibration_Bins.csv"))
write_csv_safe(cal$brier, file.path(DIR_13_REVIEWER, "Table_Brier_Scores.csv"))

for (ds in unique(all_pred$Dataset)) {
  safe_ds <- gsub("[^A-Za-z0-9_]+", "_", ds)
  plot_calibration(cal$bins, ds, file.path(DIR_13_REVIEWER, paste0("Calibration_", safe_ds, ".png")))
  plot_roc(all_pred, ds, file.path(DIR_13_REVIEWER, paste0("ROC_", safe_ds, ".png")))
}

# Permutation AUC nulls
perm_rows <- list()
perm_long_rows <- list()

for (ds in unique(all_pred$Dataset)) {
  log_message("Permutation AUC analysis: ", ds)

  sub <- all_pred[all_pred$Dataset == ds, , drop = FALSE]
  obs <- factor(sub$Observed, levels = c("Normal", "Tumor"))
  prob <- as.numeric(sub$Probability_Tumor)

  roc_obs <- tryCatch(pROC::roc(obs, prob, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE),
                      error = function(e) NULL)
  observed_auc <- if (is.null(roc_obs)) NA_real_ else as.numeric(pROC::auc(roc_obs))

  perm_auc <- rep(NA_real_, N_PERMUTATIONS)
  for (i in seq_len(N_PERMUTATIONS)) {
    obs_perm <- sample(obs)
    roc_i <- tryCatch(pROC::roc(obs_perm, prob, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE),
                      error = function(e) NULL)
    perm_auc[i] <- if (is.null(roc_i)) NA_real_ else as.numeric(pROC::auc(roc_i))
  }

  p_emp <- ifelse(is.na(observed_auc), NA_real_,
                  (sum(perm_auc >= observed_auc, na.rm = TRUE) + 1) / (sum(!is.na(perm_auc)) + 1))

  perm_rows[[ds]] <- data.frame(
    Dataset = ds,
    Observed_AUC = observed_auc,
    N_Permutations = N_PERMUTATIONS,
    Empirical_P = p_emp,
    Permuted_AUC_Mean = mean(perm_auc, na.rm = TRUE),
    Permuted_AUC_SD = sd(perm_auc, na.rm = TRUE),
    Permuted_AUC_2.5 = as.numeric(stats::quantile(perm_auc, 0.025, na.rm = TRUE)),
    Permuted_AUC_97.5 = as.numeric(stats::quantile(perm_auc, 0.975, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )

  perm_long_rows[[ds]] <- data.frame(
    Dataset = ds,
    Iteration = seq_along(perm_auc),
    Permuted_AUC = perm_auc,
    stringsAsFactors = FALSE
  )

  safe_ds <- gsub("[^A-Za-z0-9_]+", "_", ds)
  png(file.path(DIR_13_REVIEWER, paste0("Permutation_AUC_", safe_ds, ".png")),
      width = 1800, height = 1500, res = 300)
  hist(perm_auc, breaks = 30,
       main = paste0("Permutation AUC null distribution: ", ds),
       xlab = "Permuted-label AUC")
  abline(v = observed_auc, lwd = 3)
  legend("topright", legend = paste0("Observed AUC = ", sprintf("%.3f", observed_auc)), bty = "n")
  dev.off()
}

perm_tbl <- dplyr::bind_rows(perm_rows)
perm_long_tbl <- dplyr::bind_rows(perm_long_rows)

write_csv_safe(perm_tbl, file.path(DIR_13_REVIEWER, "Table_Permutation_AUC_Results.csv"))
write_csv_safe(perm_long_tbl, file.path(DIR_13_REVIEWER, "Table_Permutation_AUC_Long.csv"))

###############################################################################
# 10) REVIEWER-READY TEXT
###############################################################################

methods_text <- c(
  "Reviewer-requested uncertainty, calibration and permutation analyses",
  "",
  paste0("The locked three-miRNA Random Forest classifier was evaluated using a fixed tumour-probability decision threshold of ", DECISION_THRESHOLD, "."),
  paste0("The Random Forest model used ", RF_NTREE, " trees. When the saved final model object was unavailable, the model was reconstructed from the locked three-miRNA panel and the Train80 subset using the same split seed used in the original workflow (", SPLIT_SEED, ")."),
  "Receiver operating characteristic curves were summarized using AUC values and 95% confidence intervals estimated by DeLong's method where possible; bootstrap confidence intervals were used as fallback when required.",
  "Accuracy, sensitivity, specificity and precision were summarized with exact binomial 95% confidence intervals based on the corresponding confusion-matrix denominators.",
  paste0("Calibration was assessed by Brier score and grouped calibration tables using up to ", N_CALIBRATION_BINS, " probability bins."),
  paste0("Permutation-based null distributions were generated by randomly permuting class labels ", N_PERMUTATIONS, " times while keeping the locked model probabilities fixed; empirical P values were calculated as (number of permuted AUCs greater than or equal to the observed AUC + 1)/(number of valid permutations + 1)."),
  "",
  "Generated output folder:",
  DIR_13_REVIEWER
)

writeLines(methods_text, file.path(DIR_13_REVIEWER, "Reviewer_Methods_Text_To_Copy.txt"))

results_lines <- c(
  "Reviewer-requested analysis output summary",
  "",
  "AUC with 95% CI:",
  capture.output(print(auc_tbl)),
  "",
  "Classification metrics with exact 95% CI:",
  capture.output(print(metric_ci_tbl)),
  "",
  "Permutation AUC results:",
  capture.output(print(perm_tbl)),
  "",
  "Brier scores:",
  capture.output(print(cal$brier))
)

writeLines(results_lines, file.path(DIR_13_REVIEWER, "Reviewer_Results_Text_To_Copy.txt"))

###############################################################################
# 11) SESSION INFO AND DONE
###############################################################################

sink(file.path(DIR_13_REVIEWER, "sessionInfo_Reviewer_Requested_Uncertainty.txt"))
print(sessionInfo())
sink()

log_message("Completed successfully.")
message("Done. Reviewer-requested outputs are available in: ", DIR_13_REVIEWER)
