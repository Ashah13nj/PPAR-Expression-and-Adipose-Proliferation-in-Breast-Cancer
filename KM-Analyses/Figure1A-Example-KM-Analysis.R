#Figure 1A - example analysis

need <- c("data.table","dplyr","survival","survminer","ggplot2")
for (p in need) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(survival); library(survminer); library(ggplot2)
})

# Paths (adjust for different system or path name)
expr_path <- "/Users/ashah/Downloads/TCGA.BRCA.sampleMap_HiSeqV2"
clin_path <- "/Users/ashah/Downloads/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix"
stopifnot(file.exists(expr_path), file.exists(clin_path))

# 
TARGET_GENE <- "PPARG"
CUTOFF      <- 6.903
ROUND_DP    <- 3
PDF_OUT     <- sprintf("KM_TCGA_BRCA_%s_cutoff_%.3f.pdf", TARGET_GENE, CUTOFF)

# Helpers
is_primary_barcode <- function(bc) substr(bc, 14, 15) == "01"
num <- function(x) suppressWarnings(as.numeric(x))
chr <- function(x) suppressWarnings(as.character(x))


# Clinical
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

# Expression of Gene 
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

# Primary Tumor Selection (Barcode 01)
ALLOWED <- c("primary tumor", "primary solid tumor")
expr_primary <- expr_df %>%
  left_join(clin %>% select(sample, sample_type_norm), by="sample") %>%
  mutate(barcode01 = is_primary_barcode(sample)) %>%
  filter(barcode01, !is.na(sample_type_norm), sample_type_norm %in% ALLOWED)

# Dedup to first column order
expr_primary_unique <- expr_primary %>%
  group_by(patient) %>%
  slice_min(order_by = col_idx, n=1, with_ties=FALSE) %>%
  ungroup()
cat(sprintf("Unique primary tumor patients (pre-OS join): %d\n", nrow(expr_primary_unique)))


# Merge with OS (KEEP time==0); drop the most-negative time only

dat_km <- expr_primary_unique %>%
  inner_join(clin_os_patient, by="patient") %>%
  mutate(time = as.numeric(time), event = as.integer(event))

# Identify non-positive times for visibility
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

# Group split and KM (with CI + at-risk, styled)
dat_km <- dat_km %>%
  mutate(expr_cmp = round(expr_raw, ROUND_DP),
         group = factor(ifelse(expr_cmp >= CUTOFF, paste0("≥ ", CUTOFF), paste0("< ", CUTOFF)),
                        levels = c(paste0("< ", CUTOFF), paste0("≥ ", CUTOFF))))

counts <- table(dat_km$group); print(counts)

fit <- survfit(Surv(time, event) ~ group, data=dat_km)
lr  <- survdiff(Surv(time, event) ~ group, data=dat_km, rho=0)
pval_logrank <- 1 - pchisq(lr$chisq, df=length(lr$n)-1)
pval_logrank_fmt <- sprintf("%.3f", pval_logrank)

# Cox PH 
cox <- coxph(Surv(time, event) ~ group, data = dat_km)
cs  <- summary(cox)
hr      <- unname(cs$coefficients[1, "exp(coef)"])
hr_p    <- unname(cs$coefficients[1, "Pr(>|z|)"])
ci_low  <- unname(cs$conf.int[1, "lower .95"])
ci_high <- unname(cs$conf.int[1, "upper .95"])
hr_str  <- sprintf("HR=%.2f (95%% CI %.2f–%.2f), P = %.3f", hr, ci_low, ci_high, hr_p)

# CI ribbons
palette2 <- c(` < 6.903`="#1f77b4", `≥ 6.903`="#d62728")  

palette2 <- setNames(c("#1f77b4","#d62728"), levels(dat_km$group))
legend_labs <- sprintf("%s (n=%d)", names(counts), as.integer(counts))

max_t <- max(dat_km$time, na.rm = TRUE)
ticks  <- pretty(c(0, max_t), n = 5)
step_by <- if (length(ticks) >= 2) median(diff(ticks)) else max_t/4
if (!is.finite(step_by) || step_by <= 0) step_by <- max_t/4

km <- ggsurvplot(
  fit, data = dat_km,
  conf.int = TRUE,                 # 95% ribbons
  risk.table = TRUE,               # at-risk table
  risk.table.col   = "black",      # numbers black
  risk.table.height= 0.16,         # compact
  risk.table.y.text= FALSE,        # hide left labels
  risk.table.title = NULL,         # no "Number at risk" header
  break.time.by    = step_by,
  xlim = c(0, max_t),
  xlab = NULL,                     
  ylab = "Survival probability",
  legend.title = NULL,
  legend.labs  = legend_labs,
  palette = unname(palette2),
  ggtheme = theme_classic(base_size = 12) 
)

# Panel styling 
km$plot <- km$plot +
  theme(
    axis.text  = element_text(face = "bold", colour = "black"),
    axis.title = element_text(face = "bold", colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_blank(),
    legend.position      = c(0.02, 0.05),
    legend.justification = c(0, 0),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.text          = element_text(face = "bold")
  ) +
  annotate("label", x = Inf, y = Inf, label = hr_str,
           hjust = 1.02, vjust = 1.8, size = 4.0, label.size = 0,
           fill = "white", fontface = "bold") +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("P = %s", pval_logrank_fmt),
           hjust = 1.02, vjust = 3.2, size = 4.0, label.size = 0,
           fill = "white", fontface = "bold")

# Risk table
km$table <- km$table +
  labs(x = "Days", y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.title.x = element_text(face = "bold", colour = "black"),
    axis.text    = element_text(colour = "black"),
    axis.ticks   = element_line(colour = "black"),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = palette2, guide = "none") +
  scale_color_manual(values = palette2, guide = "none")  

print(km)

# Save PDF

pdf(PDF_OUT, width=6.5, height=5.6)
print(km)
dev.off()
cat(sprintf("Saved PDF: %s\n", normalizePath(PDF_OUT)))

# Use the loaded Clinical Matrix to stratify by Tumor Type (Figure 2) or to perform analyses by menopause status (Figure 1B-C). Use variable names as listed in source dataset, eg. Her2_Final_Status_Nature2012 or Menopause_status. Use median to separate patient cohorts by gene expression, as shown above as 6.903 for overall cohort.
