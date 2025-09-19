# TCGA-BRCA | PPARG KM + Cox (UV & MV) with 2 CSV tables

need <- c("data.table","dplyr","survival","survminer","ggplot2","broom","forcats")
for (p in need) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(survival); library(survminer); library(ggplot2)
  library(broom); library(forcats)
})

expr_path <- "/Users/ashah/Downloads/TCGA.BRCA.sampleMap_HiSeqV2"
clin_path <- "/Users/ashah/Downloads/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix"
stopifnot(file.exists(expr_path), file.exists(clin_path))

# Parameters
TARGET_GENE <- "PPARG"
CUTOFF      <- 6.903
ROUND_DP    <- 3
CSV_OUT_MV  <- sprintf("Table2_Cox_MV_%s_cutoff_%.3f.csv", TARGET_GENE, CUTOFF)

# Helper (Functions)
is_primary_barcode <- function(bc) substr(bc, 14, 15) == "01"
num <- function(x) suppressWarnings(as.numeric(x))
chr <- function(x) suppressWarnings(as.character(x))

map_meno <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[grepl("^pre", y)]  <- "pre"
  y[grepl("^peri", y)] <- "peri"
  y[grepl("^post", y)] <- "post"
  y[!(y %in% c("pre","peri","post"))] <- NA_character_
  y
}
map_posneg <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[grepl("pos", y)] <- "positive"
  y[grepl("neg", y)] <- "negative"
  y[!(y %in% c("positive","negative"))] <- NA_character_
  y
}

# Clinical Data
clin <- fread(clin_path, sep="\t", quote="", data.table=FALSE, check.names=FALSE, showProgress=FALSE)

# IDs
sample_col <- grep("^sample(?:ID)?$", names(clin), ignore.case=TRUE, value=TRUE)
if (!length(sample_col)) sample_col <- names(clin)[1]
names(clin)[match(sample_col, names(clin))] <- "sample"
clin$patient <- substr(clin$sample, 1, 12)

# Sample type (normalize)
st_candidates <- grep("^sample_type", names(clin), ignore.case=TRUE, value=TRUE)
stopifnot(length(st_candidates) >= 1)
clin <- clin %>% mutate(sample_type_norm = tolower(trimws(.data[[st_candidates[1]]])))

# OS
if (all(c("OS","OS.time") %in% names(clin))) {
  clin_os_patient <- clin %>%
    transmute(patient = substr(sample,1,12),
              event   = suppressWarnings(as.integer(OS)),
              time    = suppressWarnings(as.numeric(OS.time))) %>%
    filter(!is.na(event), !is.na(time)) %>%
    distinct(patient, .keep_all = TRUE)
} else {
  vs        <- tolower(chr(clin$vital_status))
  d_death   <- num(clin$days_to_death)
  d_alive   <- num(clin$days_to_last_known_alive)
  d_follow  <- num(clin$days_to_last_followup)
  d_contact <- if ("days_to_last_contact" %in% names(clin)) num(clin$days_to_last_contact) else NA_real_
  
  os_event <- ifelse(is.na(vs), NA_integer_,
                     ifelse(grepl("dead|deceased", vs), 1L,
                            ifelse(grepl("alive|living", vs), 0L, NA_integer_)))
  max_alive <- pmax(d_alive, d_follow, d_contact, na.rm=TRUE); max_alive[is.infinite(max_alive)] <- NA
  os_time  <- ifelse(os_event == 1L, ifelse(!is.na(d_death), d_death, max_alive), max_alive)
  
  clin_os_patient <- clin %>%
    transmute(patient = substr(sample,1,12), time=os_time, event=os_event) %>%
    filter(!is.na(time) & !is.na(event)) %>%
    distinct(patient, .keep_all = TRUE)
}
cat(sprintf("Patients with usable OS endpoint: %d\n", nrow(clin_os_patient)))

# Gene Expression
expr_dt <- fread(expr_path, sep="\t", quote="", data.table=FALSE, check.names=FALSE, showProgress=FALSE)
colnames(expr_dt)[1] <- "gene"
expr_dt$gene_clean <- toupper(sub("\\|.*$", "", expr_dt$gene))
gi <- which(expr_dt$gene_clean == toupper(TARGET_GENE)); stopifnot(length(gi) >= 1L)
row <- expr_dt[gi[1], , drop=FALSE]

value_cols <- setdiff(colnames(row), c("gene","gene_clean"))
col_idx <- seq_along(value_cols)
expr_vals  <- suppressWarnings(as.numeric(row[, value_cols, drop=TRUE])); names(expr_vals) <- value_cols

expr_df <- data.frame(
  sample   = names(expr_vals),
  expr_raw = expr_vals,
  patient  = substr(names(expr_vals), 1, 12),
  col_idx  = col_idx,
  stringsAsFactors = FALSE
) %>% filter(!is.na(expr_raw))

# Primary tumor only

ALLOWED <- c("primary tumor", "primary solid tumor")
expr_primary <- expr_df %>%
  left_join(clin %>% select(sample, sample_type_norm), by="sample") %>%
  mutate(barcode01 = is_primary_barcode(sample)) %>%
  filter(barcode01, !is.na(sample_type_norm), sample_type_norm %in% ALLOWED)

expr_primary_unique <- expr_primary %>%
  group_by(patient) %>%
  slice_min(order_by = col_idx, n=1, with_ties=FALSE) %>%
  ungroup()
cat(sprintf("Unique primary tumor patients (pre-OS join): %d\n", nrow(expr_primary_unique)))

# Drop single most-negative time
dat_km <- expr_primary_unique %>%
  inner_join(clin_os_patient, by="patient") %>%
  mutate(time = as.numeric(time), event = as.integer(event))

nonpos <- dat_km %>% filter(time <= 0) %>% arrange(time)
if (nrow(nonpos) > 0) {
  cat("Non-positive survival times (kept except the most negative one below):\n")
  print(nonpos %>% select(patient, time) %>% head(20))
}
drop_id <- if (nrow(nonpos) > 0) nonpos$patient[which.min(nonpos$time)] else NA_character_
if (!is.na(drop_id)) {
  dat_km <- dat_km %>% filter(patient != drop_id)
  cat(sprintf("Dropped patient with most-negative time to match Xena N: %s\n", drop_id))
}
cat(sprintf("Usable KM N (post-fix): %d\n", nrow(dat_km)))

# Group (used in MV Cox)
dat_km <- dat_km %>%
  mutate(expr_cmp = round(expr_raw, ROUND_DP),
         group = factor(ifelse(expr_cmp >= CUTOFF, paste0("≥ ", CUTOFF), paste0("< ", CUTOFF)),
                        levels = c(paste0("< ", CUTOFF), paste0("≥ ", CUTOFF))))


# Multivariable Cox 
table2_mv <- broom::tidy(cox_mv, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::filter(term != "menoperi") %>%   # remove "Perimenopausal (vs pre)"
  mutate(
    term = dplyr::case_when(
      term == paste0("group", levels(dat_cox$group)[2]) ~ paste0("PPARG ≥ ", CUTOFF, " (vs <)"),
      term == "age"                       ~ "Age (per year)",
      term == "menopost"                  ~ "Postmenopausal (vs pre)",
      term == "ERpositive"                ~ "ER positive (vs negative)",
      term == "PRpositive"                ~ "PR positive (vs negative)",
      term == "HER2positive"              ~ "HER2 positive (vs negative)",
      TRUE ~ term
    )
  ) %>%
  transmute(
    Term   = term,
    HR     = estimate,
    CI_low = conf.low,
    CI_high= conf.high,
    P      = p.value
  )

readr::write_csv(table2_mv, CSV_OUT_MV)
cat(sprintf("Saved Table 2 CSV: %s\n", normalizePath(CSV_OUT_MV)))
