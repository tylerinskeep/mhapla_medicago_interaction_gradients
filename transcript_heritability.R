### ESTIMATING SNP-BASED HERITABILITY FOR TRANSCRIPT ABUNDANCE ####
## Tyler Inskeep, Groen Lab, UC Riverside
## Adapted from infection-intensity heritability analysis
##
## Assumptions:
##   1. Genotypes come from variants/cohort.filtered.vcf.gz
##   2. The count matrix is in the working directory.
##   3. Count matrix format: transcripts/genes in rows, samples in columns,
##      with the first column containing transcript IDs.
##   4. Sample names in the count matrix match VCF sample names.
##
## Output:
##   transcript_heritability_results.csv
##
## IMPORTANT:
##   Raw counts are filtered and transformed to log2 CPM before modeling.
##   h2 here is SNP-based genomic heritability, not broad-sense heritability.

### SET WORKING DIRECTORY #######################################################

setwd("~/bigdata/haplaRIL")

### USER SETTINGS ##############################################################

# Change this to the actual count-matrix filename.
counts_file <- "counts.csv"

# Minimum expression filter:
# Keep transcripts with at least min_count reads in at least min_samples samples.
min_count <- 10
min_samples <- 10

# SNP QC before genomic relationship matrix construction.
max_snp_missing <- 0.10
min_maf <- 0.05

# Number of CPU cores used for transcript-wise models.
# Set to 1 for serial execution.
n_cores <- 4

# Optional covariate file.
# Set to NULL for an intercept-only fixed-effects model.
# If used, the first column must be named "sample".
covariate_file <- NULL

### LOAD PACKAGES ##############################################################

required_packages <- c(
  "vcfR",
  "rrBLUP",
  "sommer",
  "edgeR",
  "parallel"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(vcfR)
library(rrBLUP)
library(sommer)
library(edgeR)
library(parallel)

### IMPORT VCF AND BUILD NUMERIC GENOTYPE MATRIX ###############################

vcf_file <- "variants/cohort.filtered.vcf.gz"

if (!file.exists(vcf_file)) {
  stop("VCF file not found: ", vcf_file)
}

vcf <- read.vcfR(vcf_file, verbose = FALSE)

snp_meta <- as.data.frame(vcf@fix)
snp_ids <- paste0(snp_meta$CHROM, "_", snp_meta$POS)

gt_raw <- extract.gt(vcf, element = "GT")

# Robust conversion of phased and unphased diploid genotypes:
# 0/0 or 0|0 -> 0
# 0/1, 1/0, 0|1, 1|0 -> 1
# 1/1 or 1|1 -> 2
# everything else, including missing calls, -> NA
convert_gt <- function(x) {
  out <- rep(NA_real_, length(x))
  out[x %in% c("0/0", "0|0")] <- 0
  out[x %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  out[x %in% c("1/1", "1|1")] <- 2
  out
}

gt_numeric <- apply(gt_raw, 2, convert_gt)
gt_numeric <- t(gt_numeric)

colnames(gt_numeric) <- snp_ids
rownames(gt_numeric) <- colnames(gt_raw)

stopifnot(ncol(gt_numeric) == length(snp_ids))

### SNP QUALITY CONTROL ########################################################

snp_missing <- colMeans(is.na(gt_numeric))

# A.mat expects marker coding centered around -1, 0, 1 or equivalent dosage.
# For 0/1/2 dosage, allele frequency is mean dosage / 2.
allele_frequency <- colMeans(gt_numeric, na.rm = TRUE) / 2
snp_maf <- pmin(allele_frequency, 1 - allele_frequency)

keep_snp <- (
  is.finite(snp_maf) &
  snp_missing <= max_snp_missing &
  snp_maf >= min_maf
)

message(
  "Keeping ", sum(keep_snp), " of ", length(keep_snp),
  " SNPs after missingness and MAF filtering."
)

if (sum(keep_snp) < 100) {
  stop(
    "Fewer than 100 SNPs remain after QC. Inspect genotype missingness and ",
    "consider changing max_snp_missing or min_maf."
  )
}

gt_qc <- gt_numeric[, keep_snp, drop = FALSE]

### BUILD GENOMIC RELATIONSHIP MATRIX ##########################################

G <- rrBLUP::A.mat(
  gt_qc,
  impute.method = "mean"
)

rownames(G) <- rownames(gt_qc)
colnames(G) <- rownames(gt_qc)

### IMPORT COUNTS ##############################################################

if (!file.exists(counts_file)) {
  stop(
    "Count file not found: ", counts_file,
    "\nEdit counts_file in the USER SETTINGS section."
  )
}

counts_df <- read.csv(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(counts_df) < 3) {
  stop("The count matrix must contain one transcript-ID column and sample columns.")
}

transcript_ids <- counts_df[[1]]
counts <- as.matrix(counts_df[, -1, drop = FALSE])

storage.mode(counts) <- "numeric"
rownames(counts) <- make.unique(as.character(transcript_ids))

if (anyNA(counts)) {
  stop(
    "The count matrix contains missing or nonnumeric values. ",
    "Check the file format before fitting models."
  )
}

if (any(counts < 0)) {
  stop("Counts cannot contain negative values.")
}

### ALIGN RNA-SEQ SAMPLES TO GENOTYPE SAMPLES ##################################

common_samples <- intersect(colnames(counts), rownames(G))

message(
  "Found ", length(common_samples),
  " samples shared between the count matrix and genomic relationship matrix."
)

if (length(common_samples) < 20) {
  stop(
    "Too few matching samples. Check that count-matrix column names match ",
    "VCF sample names exactly."
  )
}

# Put both datasets in exactly the same order.
counts <- counts[, common_samples, drop = FALSE]
G <- G[common_samples, common_samples, drop = FALSE]

stopifnot(identical(colnames(counts), rownames(G)))
stopifnot(identical(rownames(G), colnames(G)))

### FILTER AND NORMALIZE TRANSCRIPT COUNTS #####################################

keep_transcript <- rowSums(counts >= min_count) >= min_samples

message(
  "Keeping ", sum(keep_transcript), " of ", nrow(counts),
  " transcripts after expression filtering."
)

counts_filtered <- counts[keep_transcript, , drop = FALSE]

dge <- edgeR::DGEList(counts = counts_filtered)
dge <- edgeR::calcNormFactors(dge, method = "TMM")

# log2 counts per million with a small prior count.
log_cpm <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 0.5
)

### OPTIONAL FIXED-EFFECT COVARIATES ###########################################

sample_data <- data.frame(
  sample = factor(common_samples, levels = common_samples),
  stringsAsFactors = FALSE
)

fixed_formula <- abundance ~ 1

if (!is.null(covariate_file)) {
  if (!file.exists(covariate_file)) {
    stop("Covariate file not found: ", covariate_file)
  }

  covariates <- read.csv(
    covariate_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"sample" %in% colnames(covariates)) {
    stop("The covariate file must contain a column named 'sample'.")
  }

  covariates <- covariates[
    match(common_samples, covariates$sample),
    ,
    drop = FALSE
  ]

  if (anyNA(covariates$sample)) {
    stop("Some RNA-seq samples are missing from the covariate file.")
  }

  covariate_names <- setdiff(colnames(covariates), "sample")
  sample_data <- cbind(sample_data, covariates[, covariate_names, drop = FALSE])

  fixed_formula <- as.formula(
    paste("abundance ~", paste(covariate_names, collapse = " + "))
  )
}

### FUNCTION TO FIT ONE TRANSCRIPT #############################################

fit_one_transcript <- function(i) {

  transcript_id <- rownames(log_cpm)[i]
  y <- as.numeric(log_cpm[i, ])

  dat <- sample_data
  dat$abundance <- y

  # Skip transcripts with essentially no variation after transformation.
  if (!is.finite(var(y)) || var(y) < 1e-10) {
    return(data.frame(
      transcript = transcript_id,
      n_samples = length(y),
      mean_logCPM = mean(y),
      variance_logCPM = var(y),
      genetic_variance = NA_real_,
      genetic_variance_se = NA_real_,
      residual_variance = NA_real_,
      residual_variance_se = NA_real_,
      h2 = NA_real_,
      converged = FALSE,
      boundary_estimate = NA,
      error = "No usable variance",
      stringsAsFactors = FALSE
    ))
  }

  fit <- tryCatch(
    sommer::mmer(
      fixed = fixed_formula,
      random = ~ sommer::vsr(sample, Gu = G),
      data = dat,
      verbose = FALSE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      transcript = transcript_id,
      n_samples = length(y),
      mean_logCPM = mean(y),
      variance_logCPM = var(y),
      genetic_variance = NA_real_,
      genetic_variance_se = NA_real_,
      residual_variance = NA_real_,
      residual_variance_se = NA_real_,
      h2 = NA_real_,
      converged = FALSE,
      boundary_estimate = NA,
      error = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }

  vc <- fit$varcomp

  if (is.null(vc) || nrow(vc) < 2) {
    return(data.frame(
      transcript = transcript_id,
      n_samples = length(y),
      mean_logCPM = mean(y),
      variance_logCPM = var(y),
      genetic_variance = NA_real_,
      genetic_variance_se = NA_real_,
      residual_variance = NA_real_,
      residual_variance_se = NA_real_,
      h2 = NA_real_,
      converged = FALSE,
      boundary_estimate = NA,
      error = "Variance components unavailable",
      stringsAsFactors = FALSE
    ))
  }

  vg <- as.numeric(vc[1, "VarComp"])
  ve <- as.numeric(vc[2, "VarComp"])

  vg_se <- if ("VarCompSE" %in% colnames(vc)) {
    as.numeric(vc[1, "VarCompSE"])
  } else {
    NA_real_
  }

  ve_se <- if ("VarCompSE" %in% colnames(vc)) {
    as.numeric(vc[2, "VarCompSE"])
  } else {
    NA_real_
  }

  h2 <- if (
    is.finite(vg) &&
    is.finite(ve) &&
    (vg + ve) > 0
  ) {
    vg / (vg + ve)
  } else {
    NA_real_
  }

  # sommer versions differ in how convergence information is stored.
  converged <- TRUE
  if (!is.null(fit$convergence)) {
    converged <- isTRUE(fit$convergence)
  }

  boundary_estimate <- is.finite(h2) && (h2 <= 1e-6 || h2 >= 1 - 1e-6)

  data.frame(
    transcript = transcript_id,
    n_samples = length(y),
    mean_logCPM = mean(y),
    variance_logCPM = var(y),
    genetic_variance = vg,
    genetic_variance_se = vg_se,
    residual_variance = ve,
    residual_variance_se = ve_se,
    h2 = h2,
    converged = converged,
    boundary_estimate = boundary_estimate,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
}

### FIT ALL TRANSCRIPTS ########################################################

transcript_indices <- seq_len(nrow(log_cpm))

if (.Platform$OS.type == "unix" && n_cores > 1) {
  result_list <- parallel::mclapply(
    transcript_indices,
    fit_one_transcript,
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
} else {
  result_list <- lapply(transcript_indices, fit_one_transcript)
}

results <- do.call(rbind, result_list)

### MULTIPLE-TESTING-ORIENTED DIAGNOSTIC #######################################

# A variance-component point estimate alone is not a formal significance test.
# This Z score is included only as a rough diagnostic because variance
# components lie on the boundary at zero.
results$vg_z <- with(
  results,
  genetic_variance / genetic_variance_se
)

results$vg_rough_p <- 2 * pnorm(
  -abs(results$vg_z)
)

results$vg_rough_fdr <- p.adjust(
  results$vg_rough_p,
  method = "BH"
)

# Rank transcripts by estimated h2.
results <- results[
  order(results$h2, decreasing = TRUE, na.last = TRUE),
]

write.csv(
  results,
  "transcript_heritability_results.csv",
  row.names = FALSE
)

message(
  "Finished. Results written to transcript_heritability_results.csv"
)

### QUICK SUMMARIES ############################################################

print(summary(results$h2))

cat(
  "\nSuccessfully fit:",
  sum(is.finite(results$h2)),
  "of",
  nrow(results),
  "transcripts\n"
)

cat(
  "Boundary estimates:",
  sum(results$boundary_estimate %in% TRUE, na.rm = TRUE),
  "\n"
)

head(
  results[
    ,
    c(
      "transcript",
      "h2",
      "genetic_variance",
      "genetic_variance_se",
      "residual_variance",
      "boundary_estimate",
      "error"
    )
  ],
  20
)
