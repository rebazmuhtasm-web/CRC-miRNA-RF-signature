
###############################################################################
# RBZ MASTER CLEAN PIPELINE (RECONSTRUCTED FROM CHAT WORKFLOW)
# Covers RBZ Steps 01-15
#
# Main folders created:
# 01_Input_Snapshot
# 02_Logs
# 03_QC_Preprocessing
# 04_Objects
# 05_Large_Preprocessing
# 06_Train_Test_Split
# 07_TrainOnly_Filtering_and_DE
# 08_Train_Feature_Selection
# 09_Final_Locked_Signature
# 10_Final_Locked_Modeling
# 11_Supportive_FullLarge_Locked_RF
# 12_External_HTS_Validation
# 15_Internal_Supportive_Validation
#
# Notes:
# - This is a cleaned reconstructed master script based on the RBZ stepwise
#   workflow developed in chat.
# - It is designed to match the logic and output structure of the project,
#   but it is not claimed to be byte-identical to any original step files.
# - The script assumes count-like HTS matrices for Large and Internal data.
# - External HTS datasets are expected as count matrices plus sample_info files
#   in 01_Input_Snapshot.
###############################################################################

options(stringsAsFactors = FALSE)
set.seed(20260410)

suppressPackageStartupMessages({
  pkgs <- c(
    "data.table", "dplyr", "stringr", "ggplot2", "pROC",
    "randomForest", "edgeR", "limma", "glmnet", "Boruta", "caret"
  )
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) stop("Please install these packages first: ", paste(miss, collapse = ", "))
})

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(pROC)
library(randomForest)
library(edgeR)
library(limma)
library(glmnet)
library(Boruta)
library(caret)

###############################################################################
# USER SETTINGS
###############################################################################
ROOT_DIR <- "D:/RUN_miRDeep2_edgeR_QC_20260318_180642/mixed/Raw data/RF ONLY"
PROJECT_NAME <- paste0("RBZ_STEP1_INPUT_DISCOVERY_", format(Sys.time(), "%Y%m%d_%H%M%S"))
PROJECT_DIR <- file.path(ROOT_DIR, PROJECT_NAME)

# ---- INPUT FILES (edit only if needed) ----
INPUT_LARGE_COUNT  <- file.path(ROOT_DIR, "Large_COUNT_MATRIX_raw.csv.csv")
INPUT_LARGE_META   <- file.path(ROOT_DIR, "sample_info_large.csv.csv")
INPUT_INTERNAL_COUNT <- file.path(ROOT_DIR, "Internal_COUNT_MATRIX_raw.csv")
INPUT_INTERNAL_META  <- file.path(ROOT_DIR, "sample_info_internal.csv")

# External HTS validation datasets expected under ROOT_DIR or already copied
INPUT_EXT_FILES <- list(
  list(id = "GSE125961", count = file.path(ROOT_DIR, "GSE125961_count_matrix_raw.csv"),
       meta = file.path(ROOT_DIR, "sample_info_internal_GSE125961.csv")),
  list(id = "GSE183437", count = file.path(ROOT_DIR, "GSE183437_count_matrix_raw.csv"),
       meta = file.path(ROOT_DIR, "sample_info_internal_GSE183437.csv")),
  list(id = "GSE236991", count = file.path(ROOT_DIR, "GSE236991_count_matrix_raw.csv"),
       meta = file.path(ROOT_DIR, "sample_info_internal_GSE236991.csv"))
)

# ---- ANALYSIS PARAMETERS ----
TRAIN_PROP <- 0.80
N_FOLDS <- 5
N_REPEATS <- 10
FREQ_THRESHOLD_STRICT <- 0.70
MIN_FINAL_FEATURES <- 3
MAX_FINAL_FEATURES <- 10
POS_CLASS <- "Tumor"
NEG_CLASS <- "Normal"
CPM_MIN <- 1
CPM_MIN_SAMPLES <- 2
RF_NTREE <- 1000
PLOT_DPI <- 300

###############################################################################
# DIRECTORY STRUCTURE
###############################################################################
DIR_01_INPUT <- file.path(PROJECT_DIR, "01_Input_Snapshot")
DIR_02_LOGS  <- file.path(PROJECT_DIR, "02_Logs")
DIR_03_QC    <- file.path(PROJECT_DIR, "03_QC_Preprocessing")
DIR_04_OBJ   <- file.path(PROJECT_DIR, "04_Objects")
DIR_05_LARGE <- file.path(PROJECT_DIR, "05_Large_Preprocessing")
DIR_06_SPLIT <- file.path(PROJECT_DIR, "06_Train_Test_Split")
DIR_07_DE    <- file.path(PROJECT_DIR, "07_TrainOnly_Filtering_and_DE")
DIR_08_FS    <- file.path(PROJECT_DIR, "08_Train_Feature_Selection")
DIR_09_SIG   <- file.path(PROJECT_DIR, "09_Final_Locked_Signature")
DIR_10_MOD   <- file.path(PROJECT_DIR, "10_Final_Locked_Modeling")
DIR_11_SUPP  <- file.path(PROJECT_DIR, "11_Supportive_FullLarge_Locked_RF")
DIR_12_EXT   <- file.path(PROJECT_DIR, "12_External_HTS_Validation")
DIR_15_INT   <- file.path(PROJECT_DIR, "15_Internal_Supportive_Validation")

DIRS <- c(
  PROJECT_DIR, DIR_01_INPUT, DIR_02_LOGS, DIR_03_QC, DIR_04_OBJ, DIR_05_LARGE,
  DIR_06_SPLIT, DIR_07_DE, DIR_08_FS, DIR_09_SIG, DIR_10_MOD, DIR_11_SUPP,
  DIR_12_EXT, DIR_15_INT
)
for (d in DIRS) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(DIR_02_LOGS, "RBZ_master_pipeline_log.txt")
log_message <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

###############################################################################
# HELPERS
###############################################################################
write_csv6 <- function(df, file, row.names = FALSE) {
  df2 <- as.data.frame(df, stringsAsFactors = FALSE)
  num_idx <- vapply(df2, is.numeric, logical(1))
  if (any(num_idx)) {
    df2[num_idx] <- lapply(df2[num_idx], function(x) ifelse(is.na(x), NA_character_, sprintf("%.6f", x)))
  }
  utils::write.csv(df2, file = file, row.names = row.names, na = "")
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
  cand <- c("Group","group","Class","class","Condition","condition","Type","type","status","Status")
  hit <- cand[cand %in% colnames(meta)][1]
  if (!is.na(hit)) return(hit)
  for (nm in colnames(meta)) {
    vals <- unique(tolower(as.character(meta[[nm]])))
    if (any(grepl("tumou?r|cancer|crc|case", vals)) && any(grepl("normal|control|adjacent", vals))) return(nm)
  }
  NA_character_
}

detect_pair_column <- function(meta) {
  cand <- c("pair_id","PairID","pair","pairID","pair.id")
  hit <- cand[cand %in% colnames(meta)][1]
  if (!is.na(hit)) return(hit)
  NA_character_
}

read_count_dataset <- function(count_file, meta_file, dataset_name) {
  stopifnot(file.exists(count_file), file.exists(meta_file))
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
  keep <- meta$Group %in% c("Normal","Tumor")
  meta <- meta[keep, , drop = FALSE]
  cnt_mat <- cnt_mat[, meta$SampleID, drop = FALSE]

  list(dataset_name = dataset_name, counts = cnt_mat, info = meta)
}

calc_tmm_logcpm <- function(counts_mat) {
  dge <- DGEList(counts = counts_mat)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm(dge, log = TRUE, prior.count = 1)
}

zscore_by_row <- function(mat) {
  out <- t(scale(t(mat)))
  out[is.na(out)] <- 0
  out
}

calc_metrics <- function(obs, prob, positive = "Tumor", threshold = 0.5) {
  obs <- factor(obs, levels = c("Normal","Tumor"))
  pred <- ifelse(prob >= threshold, positive, "Normal")
  pred <- factor(pred, levels = c("Normal","Tumor"))
  tp <- sum(obs == positive & pred == positive, na.rm = TRUE)
  tn <- sum(obs != positive & pred != positive, na.rm = TRUE)
  fp <- sum(obs != positive & pred == positive, na.rm = TRUE)
  fn <- sum(obs == positive & pred != positive, na.rm = TRUE)
  sensitivity <- ifelse((tp+fn)==0, NA, tp/(tp+fn))
  specificity <- ifelse((tn+fp)==0, NA, tn/(tn+fp))
  precision   <- ifelse((tp+fp)==0, NA, tp/(tp+fp))
  accuracy <- mean(obs == pred, na.rm = TRUE)
  f1 <- ifelse(is.na(precision) || is.na(sensitivity) || (precision+sensitivity)==0, NA,
               2*precision*sensitivity/(precision+sensitivity))
  roc_obj <- tryCatch(
    pROC::roc(response = obs, predictor = prob, levels = c("Normal","Tumor"), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  auc_val <- if (is.null(roc_obj)) NA_real_ else as.numeric(pROC::auc(roc_obj))
  list(
    metrics = data.frame(
      N = length(obs), Normal_n = sum(obs=="Normal"), Tumor_n = sum(obs=="Tumor"),
      AUC = auc_val, TP = tp, TN = tn, FP = fp, FN = fn,
      Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity,
      Precision = precision, F1 = f1, stringsAsFactors = FALSE
    ),
    roc = roc_obj
  )
}

make_grouped_folds <- function(group, pair_id, k = 5, repeats = 10, seed = 1) {
  set.seed(seed)
  out <- list()
  if (all(!is.na(pair_id))) {
    pair_df <- data.frame(pair_id = pair_id, group = group, stringsAsFactors = FALSE) |>
      distinct(pair_id, .keep_all = TRUE)
    idx_pairs <- createMultiFolds(y = factor(pair_df$group, levels = c("Normal","Tumor")),
                                  k = k, times = repeats)
    for (nm in names(idx_pairs)) {
      tr_pairs <- pair_df$pair_id[idx_pairs[[nm]]]
      tr <- which(pair_id %in% tr_pairs)
      te <- setdiff(seq_along(group), tr)
      out[[nm]] <- list(train = tr, test = te)
    }
  } else {
    idx <- createMultiFolds(y = factor(group, levels = c("Normal","Tumor")), k = k, times = repeats)
    for (nm in names(idx)) {
      tr <- idx[[nm]]
      te <- setdiff(seq_along(group), tr)
      out[[nm]] <- list(train = tr, test = te)
    }
  }
  out
}

###############################################################################
# STEP 01: INPUT SNAPSHOT
###############################################################################
log_message("STEP 01: Input snapshot")

copy_if_exists <- function(src, dst_dir) {
  if (file.exists(src)) file.copy(src, file.path(dst_dir, basename(src)), overwrite = TRUE)
}

copy_if_exists(INPUT_LARGE_COUNT, DIR_01_INPUT)
copy_if_exists(INPUT_LARGE_META, DIR_01_INPUT)
copy_if_exists(INPUT_INTERNAL_COUNT, DIR_01_INPUT)
copy_if_exists(INPUT_INTERNAL_META, DIR_01_INPUT)
for (sp in INPUT_EXT_FILES) {
  copy_if_exists(sp$count, DIR_01_INPUT)
  copy_if_exists(sp$meta, DIR_01_INPUT)
}

input_manifest <- data.frame(
  file = c(INPUT_LARGE_COUNT, INPUT_LARGE_META, INPUT_INTERNAL_COUNT, INPUT_INTERNAL_META,
           unlist(lapply(INPUT_EXT_FILES, function(x) c(x$count, x$meta)))),
  exists = file.exists(c(INPUT_LARGE_COUNT, INPUT_LARGE_META, INPUT_INTERNAL_COUNT, INPUT_INTERNAL_META,
                         unlist(lapply(INPUT_EXT_FILES, function(x) c(x$count, x$meta))))),
  stringsAsFactors = FALSE
)
write_csv6(input_manifest, file.path(DIR_01_INPUT, "Table_input_manifest.csv"), row.names = FALSE)

###############################################################################
# STEP 03 + 05: QC + LARGE PREPROCESSING
###############################################################################
log_message("STEP 03/05: QC preprocessing and large preprocessing")

large_before <- read_count_dataset(INPUT_LARGE_COUNT, INPUT_LARGE_META, "Large")
large_after <- large_before
large_after$expr_native <- calc_tmm_logcpm(large_after$counts)

qc_large <- data.frame(
  dataset = "Large",
  n_samples = nrow(large_after$info),
  n_tumor = sum(large_after$info$Group == "Tumor"),
  n_normal = sum(large_after$info$Group == "Normal"),
  n_miRNAs = nrow(large_after$expr_native),
  scale = "TMM_logCPM",
  stringsAsFactors = FALSE
)
write_csv6(qc_large, file.path(DIR_03_QC, "Table_QC_Summary_Large.csv"))

saveRDS(list(large_before = large_before, large_after = large_after),
        file = file.path(DIR_04_OBJ, "RBZ_step3_large_preprocessed.rds"))

write_csv6(large_after$expr_native, file.path(DIR_05_LARGE, "Matrix_Large_TMM_logCPM.csv"), row.names = TRUE)
write_csv6(large_after$info, file.path(DIR_05_LARGE, "sample_info_large_preprocessed.csv"), row.names = FALSE)

###############################################################################
# STEP 06: TRAIN / TEST SPLIT (PAIRED)
###############################################################################
log_message("STEP 06: Paired train/test split")

meta_large <- large_after$info
pair_available <- all(!is.na(meta_large$pair_id))
if (!pair_available) stop("Large dataset metadata must contain pair_id for paired split.")

pair_df <- meta_large |>
  distinct(pair_id, .keep_all = TRUE)

set.seed(20260405)
train_pairs <- createDataPartition(factor(pair_df$Group, levels = c("Normal","Tumor")), p = TRAIN_PROP, list = FALSE)
train_pair_ids <- pair_df$pair_id[train_pairs]
test_pair_ids  <- setdiff(pair_df$pair_id, train_pair_ids)

train_idx <- which(meta_large$pair_id %in% train_pair_ids)
test_idx  <- which(meta_large$pair_id %in% test_pair_ids)

train_ds <- list(
  expr_native = large_after$expr_native[, train_idx, drop = FALSE],
  counts = large_after$counts[, train_idx, drop = FALSE],
  info = meta_large[train_idx, , drop = FALSE],
  dataset_name = "Large_Train80"
)
test_ds <- list(
  expr_native = large_after$expr_native[, test_idx, drop = FALSE],
  counts = large_after$counts[, test_idx, drop = FALSE],
  info = meta_large[test_idx, , drop = FALSE],
  dataset_name = "Large_Test20"
)

write_csv6(train_ds$info, file.path(DIR_06_SPLIT, "sample_info_train80.csv"))
write_csv6(test_ds$info,  file.path(DIR_06_SPLIT, "sample_info_test20.csv"))
write_csv6(train_ds$expr_native, file.path(DIR_06_SPLIT, "Matrix_Train80_TMM_logCPM.csv"), row.names = TRUE)
write_csv6(test_ds$expr_native,  file.path(DIR_06_SPLIT, "Matrix_Test20_TMM_logCPM.csv"), row.names = TRUE)

saveRDS(list(train = train_ds, test = test_ds),
        file = file.path(DIR_04_OBJ, "RBZ_step4_train_test_split.rds"))

###############################################################################
# STEP 07: TRAIN-ONLY FILTERING AND DE
###############################################################################
log_message("STEP 07: Train-only filtering and DE")

train_counts <- train_ds$counts
train_group <- factor(train_ds$info$Group, levels = c("Normal","Tumor"))

keep <- rowSums(edgeR::cpm(train_counts) >= CPM_MIN) >= CPM_MIN_SAMPLES
train_counts_filt <- train_counts[keep, , drop = FALSE]
dge <- DGEList(counts = train_counts_filt, group = train_group)
dge <- calcNormFactors(dge)
design <- model.matrix(~ train_group)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
de_tbl <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
de_tbl$miRNA <- rownames(de_tbl)
de_tbl <- de_tbl[, c("miRNA","logFC","AveExpr","t","P.Value","adj.P.Val","B")]

de_sig <- de_tbl |> filter(adj.P.Val <= 0.05 & abs(logFC) >= 1)

write_csv6(de_tbl, file.path(DIR_07_DE, "Table_Train80_All_DE.csv"))
write_csv6(de_sig, file.path(DIR_07_DE, "Table_Train80_Significant_DE.csv"))
write_csv6(v$E, file.path(DIR_07_DE, "Matrix_Train80_Filtered_voomExpression.csv"), row.names = TRUE)

saveRDS(list(
  keep_features = rownames(train_counts_filt),
  voom_expr = v$E,
  de_all = de_tbl,
  de_sig = de_sig
), file = file.path(DIR_04_OBJ, "RBZ_step5_train_only_filtering_and_de.rds"))

###############################################################################
# STEP 08: TRAIN FEATURE SELECTION (REPEATED CV)
###############################################################################
log_message("STEP 08: Train feature selection")

step5 <- readRDS(file.path(DIR_04_OBJ, "RBZ_step5_train_only_filtering_and_de.rds"))
train_expr <- step5$voom_expr
train_meta <- train_ds$info
folds <- make_grouped_folds(
  group = train_meta$Group,
  pair_id = train_meta$pair_id,
  k = N_FOLDS,
  repeats = N_REPEATS,
  seed = 20260405
)

selector_records <- list()

for (nm in names(folds)) {
  fld <- folds[[nm]]
  tr <- fld$train
  expr_tr <- train_expr[, tr, drop = FALSE]
  group_tr <- factor(train_meta$Group[tr], levels = c("Normal","Tumor"))

  # DE inside CV fold
  keep_fold <- rowSums(expr_tr > 0) >= CPM_MIN_SAMPLES
  expr_fold <- expr_tr[keep_fold, , drop = FALSE]
  design_fold <- model.matrix(~ group_tr)
  fit_fold <- eBayes(lmFit(expr_fold, design_fold))
  de_fold <- topTable(fit_fold, coef = 2, number = Inf, sort.by = "P")
  de_fold$miRNA <- rownames(de_fold)
  de_candidates <- de_fold |>
    filter(adj.P.Val <= 0.05 & abs(logFC) >= 1) |>
    pull(miRNA) |>
    unique()

  if (length(de_candidates) < 3) next

  x <- t(expr_fold[de_candidates, , drop = FALSE])
  y <- ifelse(group_tr == "Tumor", 1, 0)

  # Elastic Net
  cvfit <- tryCatch(cv.glmnet(x, y, family = "binomial", alpha = 0.5, nfolds = min(5, nrow(x))), error = function(e) NULL)
  enet_sel <- character(0)
  if (!is.null(cvfit)) {
    cf <- coef(cvfit, s = "lambda.min")
    rn <- rownames(cf)
    enet_sel <- rn[as.vector(cf != 0)]
    enet_sel <- setdiff(enet_sel, "(Intercept)")
  }

  # Boruta
  boruta_sel <- character(0)
  bor_res <- tryCatch(Boruta(x = as.data.frame(x), y = factor(group_tr), doTrace = 0), error = function(e) NULL)
  if (!is.null(bor_res)) {
    bor_final <- TentativeRoughFix(bor_res)
    boruta_sel <- getSelectedAttributes(bor_final, withTentative = FALSE)
  }

  # RF-RFE via caret
  rfe_sel <- character(0)
  ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = min(5, nrow(x)))
  sizes <- sort(unique(c(3, 5, 10, min(20, ncol(x)), ncol(x))))
  rfe_res <- tryCatch(rfe(x = x, y = factor(group_tr), sizes = sizes, rfeControl = ctrl), error = function(e) NULL)
  if (!is.null(rfe_res)) rfe_sel <- predictors(rfe_res)

  add_rec <- function(method, feats) {
    if (length(feats) == 0) return(NULL)
    data.frame(Fold = nm, Method = method, miRNA = feats, stringsAsFactors = FALSE)
  }

  selector_records[[length(selector_records) + 1]] <- add_rec("ENet", enet_sel)
  selector_records[[length(selector_records) + 1]] <- add_rec("Boruta", boruta_sel)
  selector_records[[length(selector_records) + 1]] <- add_rec("RFRFE", rfe_sel)
}

selector_records <- selector_records[!vapply(selector_records, is.null, logical(1))]
selector_tbl <- bind_rows(selector_records)

write_csv6(selector_tbl, file.path(DIR_08_FS, "Table_Train80_Selector_Records_Long.csv"))

# Summaries
total_runs <- length(unique(selector_tbl$Fold))
freq_tbl <- selector_tbl |>
  count(Method, miRNA, name = "Selection_Count") |>
  mutate(Selection_Rate = Selection_Count / total_runs)

wide_tbl <- freq_tbl |>
  select(Method, miRNA, Selection_Rate) |>
  tidyr::pivot_wider(names_from = Method, values_from = Selection_Rate, values_fill = 0)

wide_tbl$Selection_Count <- rowSums(wide_tbl[, setdiff(colnames(wide_tbl), "miRNA"), drop = FALSE] >= FREQ_THRESHOLD_STRICT)
wide_tbl$StrictIntersection_Rate <- apply(wide_tbl[, c("ENet","Boruta","RFRFE"), drop = FALSE], 1, min)
wide_tbl$Average_Method_Rate <- rowMeans(wide_tbl[, c("ENet","Boruta","RFRFE"), drop = FALSE])

write_csv6(freq_tbl, file.path(DIR_08_FS, "Table_Train80_Feature_Selection_Frequency_Long.csv"))
write_csv6(wide_tbl, file.path(DIR_08_FS, "Table_Train80_Feature_Selection_Frequency_Wide.csv"))

saveRDS(list(selector_tbl = selector_tbl, frequency_wide = wide_tbl),
        file = file.path(DIR_04_OBJ, "RBZ_step8_train_feature_selection.rds"))

###############################################################################
# STEP 09: FINAL LOCKED SIGNATURE
###############################################################################
log_message("STEP 09: Final locked signature")

step8 <- readRDS(file.path(DIR_04_OBJ, "RBZ_step8_train_feature_selection.rds"))
wide_tbl <- step8$frequency_wide
de_sig <- step5$de_sig

final_signature <- wide_tbl |>
  filter(ENet >= FREQ_THRESHOLD_STRICT & Boruta >= FREQ_THRESHOLD_STRICT & RFRFE >= FREQ_THRESHOLD_STRICT) |>
  arrange(desc(StrictIntersection_Rate), desc(Average_Method_Rate))

if (nrow(final_signature) < MIN_FINAL_FEATURES) {
  final_signature <- wide_tbl |>
    filter((ENet >= FREQ_THRESHOLD_STRICT) + (Boruta >= FREQ_THRESHOLD_STRICT) + (RFRFE >= FREQ_THRESHOLD_STRICT) >= 2) |>
    arrange(desc(Average_Method_Rate))
}

final_signature <- final_signature |>
  left_join(de_sig[, c("miRNA","logFC","adj.P.Val","P.Value","AveExpr","B")], by = "miRNA") |>
  mutate(abs_logFC = abs(logFC)) |>
  arrange(desc(StrictIntersection_Rate), desc(Average_Method_Rate), adj.P.Val, desc(abs_logFC)) |>
  head(3)

decision_tbl <- data.frame(
  Decision_Mode = ifelse(all(final_signature$StrictIntersection_Rate >= FREQ_THRESHOLD_STRICT), "Primary_StrictIntersection", "Fallback_AverageRate"),
  Strict_Rate_Threshold = FREQ_THRESHOLD_STRICT,
  Fallback_Average_Rate_Threshold = FREQ_THRESHOLD_STRICT,
  Min_Final_Features = MIN_FINAL_FEATURES,
  Max_Final_Features = MAX_FINAL_FEATURES,
  N_Train_DE_Candidates = nrow(de_sig),
  N_Strict_Stable_Set = sum(wide_tbl$StrictIntersection_Rate >= FREQ_THRESHOLD_STRICT),
  N_Fallback_Stable_Set = sum(((wide_tbl$ENet >= FREQ_THRESHOLD_STRICT) + (wide_tbl$Boruta >= FREQ_THRESHOLD_STRICT) + (wide_tbl$RFRFE >= FREQ_THRESHOLD_STRICT)) >= 2),
  N_Final_Locked_Signature = nrow(final_signature),
  stringsAsFactors = FALSE
)

write_csv6(decision_tbl, file.path(DIR_09_SIG, "Table_Final_Locked_Decision_Summary.csv"))
write_csv6(final_signature, file.path(DIR_09_SIG, "Table_Final_Locked_Signature.csv"))

saveRDS(list(final_signature = final_signature, decision_summary = decision_tbl),
        file = file.path(DIR_04_OBJ, "RBZ_step7_final_locked_signature.rds"))

###############################################################################
# STEP 10: FINAL LOCKED MODELING
###############################################################################
log_message("STEP 10: Final locked modeling")

locked_mirs <- final_signature$miRNA
x_train <- t(train_ds$expr_native[locked_mirs, , drop = FALSE])
x_test  <- t(test_ds$expr_native[locked_mirs, , drop = FALSE])
y_train <- factor(train_ds$info$Group, levels = c("Normal","Tumor"))
y_test  <- factor(test_ds$info$Group, levels = c("Normal","Tumor"))

rf_fit <- randomForest(x = x_train, y = y_train, ntree = RF_NTREE, importance = TRUE)

train_prob <- predict(rf_fit, newdata = x_train, type = "prob")[, "Tumor"]
test_prob  <- predict(rf_fit, newdata = x_test, type = "prob")[, "Tumor"]

train_perf <- calc_metrics(y_train, train_prob, threshold = 0.5)
test_perf  <- calc_metrics(y_test, test_prob, threshold = 0.5)

panel_tbl <- bind_rows(
  cbind(data.frame(Set = "Train80", Model = "Locked_RF_Panel"), train_perf$metrics),
  cbind(data.frame(Set = "Test20", Model = "Locked_RF_Panel"), test_perf$metrics)
)
write_csv6(panel_tbl, file.path(DIR_10_MOD, "Table_Final_Locked_Panel_Performance.csv"))

indiv_rows <- list()
for (mir in locked_mirs) {
  tr_res <- calc_metrics(y_train, as.numeric(x_train[, mir]), threshold = 0)
  te_res <- calc_metrics(y_test, as.numeric(x_test[, mir]), threshold = 0)
  if (!is.na(tr_res$metrics$AUC[1]) && tr_res$metrics$AUC[1] < 0.5) tr_res <- calc_metrics(y_train, -as.numeric(x_train[, mir]), threshold = 0)
  if (!is.na(te_res$metrics$AUC[1]) && te_res$metrics$AUC[1] < 0.5) te_res <- calc_metrics(y_test, -as.numeric(x_test[, mir]), threshold = 0)
  indiv_rows[[length(indiv_rows)+1]] <- cbind(data.frame(Set = "Train80", miRNA = mir, Model = "Individual_Descriptive"), tr_res$metrics)
  indiv_rows[[length(indiv_rows)+1]] <- cbind(data.frame(Set = "Test20", miRNA = mir, Model = "Individual_Descriptive"), te_res$metrics)
}
indiv_tbl <- bind_rows(indiv_rows)
write_csv6(indiv_tbl, file.path(DIR_10_MOD, "Table_Individual_FinalLocked_miRNA_Descriptive_Performance.csv"))

saveRDS(list(rf_fit = rf_fit, panel_summary = panel_tbl, individual_summary = indiv_tbl),
        file = file.path(DIR_04_OBJ, "RBZ_step10_final_locked_modeling.rds"))

###############################################################################
# STEP 11: SUPPORTIVE FULL LARGE LOCKED RF
###############################################################################
log_message("STEP 11: Supportive full large locked RF")

x_full <- t(large_after$expr_native[locked_mirs, , drop = FALSE])
y_full <- factor(large_after$info$Group, levels = c("Normal","Tumor"))

full_prob <- predict(rf_fit, newdata = x_full, type = "prob")[, "Tumor"]
full_perf <- calc_metrics(y_full, full_prob, threshold = 0.5)

full_panel_tbl <- cbind(data.frame(Dataset = "FullLarge", Model = "Locked_RF_Panel"), full_perf$metrics)
write_csv6(full_panel_tbl, file.path(DIR_11_SUPP, "Table_Supportive_FullLarge_Locked_Panel_Performance.csv"))

full_indiv_rows <- list()
for (mir in locked_mirs) {
  res <- calc_metrics(y_full, as.numeric(x_full[, mir]), threshold = 0)
  if (!is.na(res$metrics$AUC[1]) && res$metrics$AUC[1] < 0.5) res <- calc_metrics(y_full, -as.numeric(x_full[, mir]), threshold = 0)
  full_indiv_rows[[length(full_indiv_rows)+1]] <- cbind(data.frame(Dataset = "FullLarge", miRNA = mir, Model = "Individual_Descriptive"), res$metrics)
}
full_indiv_tbl <- bind_rows(full_indiv_rows)
write_csv6(full_indiv_tbl, file.path(DIR_11_SUPP, "Table_Supportive_FullLarge_Individual_miRNA_Performance.csv"))

# DE on full large for locked miRs
design_full <- model.matrix(~ y_full)
fit_full <- eBayes(lmFit(large_after$expr_native[locked_mirs, , drop = FALSE], design_full))
de_full <- topTable(fit_full, coef = 2, number = Inf, sort.by = "P")
de_full$miRNA <- rownames(de_full)
de_full$Mean_Normal <- rowMeans(large_after$expr_native[locked_mirs, y_full == "Normal", drop = FALSE])
de_full$Mean_Tumor  <- rowMeans(large_after$expr_native[locked_mirs, y_full == "Tumor", drop = FALSE])
de_full <- de_full[, c("miRNA","Mean_Normal","Mean_Tumor","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
write_csv6(de_full, file.path(DIR_11_SUPP, "Table_Supportive_FullLarge_FinalLocked_DE_Normal_vs_Tumor.csv"))

###############################################################################
# STEP 12: EXTERNAL HTS VALIDATION
###############################################################################
log_message("STEP 12: External HTS validation")

external_panel_rows <- list()
external_indiv_rows <- list()

for (sp in INPUT_EXT_FILES) {
  if (!file.exists(sp$count) || !file.exists(sp$meta)) next
  ds <- read_count_dataset(sp$count, sp$meta, sp$id)
  ds$expr_native <- calc_tmm_logcpm(ds$counts)

  common_mirs <- intersect(locked_mirs, rownames(ds$expr_native))
  if (length(common_mirs) < length(locked_mirs)) {
    log_message("Dataset ", sp$id, " is missing some locked miRNAs; skipping if < 3.")
    if (length(common_mirs) < 3) next
  }

  x_ext <- t(ds$expr_native[locked_mirs[locked_mirs %in% rownames(ds$expr_native)], , drop = FALSE])
  y_ext <- factor(ds$info$Group, levels = c("Normal","Tumor"))

  # If exactly the same features exist, use RF; otherwise use z-score mean fallback
  if (all(locked_mirs %in% rownames(ds$expr_native))) {
    prob <- predict(rf_fit, newdata = x_ext, type = "prob")[, "Tumor"]
    res_panel <- calc_metrics(y_ext, prob, threshold = 0.5)
    external_panel_rows[[length(external_panel_rows)+1]] <- cbind(
      data.frame(Dataset = sp$id, Validation_Type = "External_HTS", Model = "Locked_RF_Panel"), res_panel$metrics
    )
  } else {
    zmat <- zscore_by_row(ds$expr_native[common_mirs, , drop = FALSE])
    panel_score <- colMeans(zmat)
    res_panel <- calc_metrics(y_ext, panel_score, threshold = 0)
    external_panel_rows[[length(external_panel_rows)+1]] <- cbind(
      data.frame(Dataset = sp$id, Validation_Type = "External_HTS", Model = "Locked_ZscoreMean_Panel"), res_panel$metrics
    )
  }

  for (mir in intersect(locked_mirs, rownames(ds$expr_native))) {
    res <- calc_metrics(y_ext, as.numeric(ds$expr_native[mir, ]), threshold = 0)
    if (!is.na(res$metrics$AUC[1]) && res$metrics$AUC[1] < 0.5) res <- calc_metrics(y_ext, -as.numeric(ds$expr_native[mir, ]), threshold = 0)
    external_indiv_rows[[length(external_indiv_rows)+1]] <- cbind(
      data.frame(Dataset = sp$id, Validation_Type = "External_HTS", miRNA = mir, Model = "Individual_Descriptive"), res$metrics
    )
  }
}

if (length(external_panel_rows) > 0) {
  ext_panel_tbl <- bind_rows(external_panel_rows)
  write_csv6(ext_panel_tbl, file.path(DIR_12_EXT, "Table_External_HTS_Panel_Performance.csv"))
}
if (length(external_indiv_rows) > 0) {
  ext_indiv_tbl <- bind_rows(external_indiv_rows)
  write_csv6(ext_indiv_tbl, file.path(DIR_12_EXT, "Table_External_HTS_Individual_miRNA_Performance.csv"))
}

###############################################################################
# STEP 15: INTERNAL SUPPORTIVE VALIDATION
###############################################################################
log_message("STEP 15: Internal supportive validation")

internal_ds <- read_count_dataset(INPUT_INTERNAL_COUNT, INPUT_INTERNAL_META, "Internal")
internal_ds$expr_native <- calc_tmm_logcpm(internal_ds$counts)
common_mirs <- intersect(locked_mirs, rownames(internal_ds$expr_native))
if (length(common_mirs) < 1) stop("No locked miRNAs found in internal dataset.")

x_int <- t(internal_ds$expr_native[common_mirs, , drop = FALSE])
y_int <- factor(internal_ds$info$Group, levels = c("Normal","Tumor"))

# Supportive panel by z-score mean to avoid mismatch if not all features or scale differences
zmat_int <- zscore_by_row(internal_ds$expr_native[common_mirs, , drop = FALSE])
panel_score_int <- colMeans(zmat_int)
int_panel <- calc_metrics(y_int, panel_score_int, threshold = 0)

int_panel_tbl <- cbind(
  data.frame(Dataset = "Internal", Validation_Type = "Supportive_Internal", Model = "Locked_Panel_ZscoreMean"),
  int_panel$metrics
)
write_csv6(int_panel_tbl, file.path(DIR_15_INT, "Table_Internal_Supportive_Panel_Performance.csv"))

int_indiv_rows <- list()
for (mir in common_mirs) {
  res <- calc_metrics(y_int, as.numeric(internal_ds$expr_native[mir, ]), threshold = 0)
  if (!is.na(res$metrics$AUC[1]) && res$metrics$AUC[1] < 0.5) res <- calc_metrics(y_int, -as.numeric(internal_ds$expr_native[mir, ]), threshold = 0)
  int_indiv_rows[[length(int_indiv_rows)+1]] <- cbind(
    data.frame(Dataset = "Internal", Validation_Type = "Supportive_Internal", miRNA = mir, Model = "Individual_Descriptive"),
    res$metrics
  )
}
int_indiv_tbl <- bind_rows(int_indiv_rows)
write_csv6(int_indiv_tbl, file.path(DIR_15_INT, "Table_Internal_Supportive_Individual_miRNA_Performance.csv"))

design_int <- model.matrix(~ y_int)
fit_int <- eBayes(lmFit(internal_ds$expr_native[common_mirs, , drop = FALSE], design_int))
de_int <- topTable(fit_int, coef = 2, number = Inf, sort.by = "P")
de_int$miRNA <- rownames(de_int)
de_int$Mean_Normal <- rowMeans(internal_ds$expr_native[common_mirs, y_int == "Normal", drop = FALSE])
de_int$Mean_Tumor  <- rowMeans(internal_ds$expr_native[common_mirs, y_int == "Tumor", drop = FALSE])
de_int <- de_int[, c("miRNA","Mean_Normal","Mean_Tumor","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
write_csv6(de_int, file.path(DIR_15_INT, "Table_Internal_Supportive_FinalLocked_DE_Normal_vs_Tumor.csv"))

###############################################################################
# SAVE MASTER OBJECT + SESSION INFO
###############################################################################
log_message("Saving master outputs")

rbz_master <- list(
  project_dir = PROJECT_DIR,
  large_after = large_after,
  train_ds = train_ds,
  test_ds = test_ds,
  de_train = step5,
  feature_selection = step8,
  final_signature = final_signature,
  rf_model = rf_fit,
  panel_train_test = panel_tbl,
  panel_full_large = full_panel_tbl,
  internal_panel = int_panel_tbl
)

saveRDS(rbz_master, file = file.path(DIR_04_OBJ, "RBZ_master_steps01_15_pipeline_object.rds"))
capture.output(sessionInfo(), file = file.path(DIR_02_LOGS, "RBZ_master_steps01_15_sessionInfo.txt"))

cat("\n==================== DONE ====================\n")
cat("RBZ master reconstructed pipeline completed.\n")
cat("Project folder:\n", PROJECT_DIR, "\n")
cat("Locked miRNAs:\n", paste(locked_mirs, collapse = ", "), "\n")
