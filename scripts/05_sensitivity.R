# =============================================================================
# 05_sensitivity.R -- GSE270472 reanalysis, KO-NSC vs IC1-NSC
#
# Three sensitivity analyses, all registered in METHODS_LOG §9 Extensions:
#   A. independentFiltering = FALSE  -- direct test of §5.13
#   B. leave-one-out dropping KO_4   -- robustness to the PC2 outlier (§9)
#   C. continuous passage covariate  -- exploratory (§5.10)
#
# Run from the repo root, after 03.
# Writes: data/processed/05_sensitivity.rds
#         results/tables/05_sensitivity_summary.csv
#         results/figures/05_*.png
#
# Cross-references like (PREP 7) point at notes/prep_python_to_r.R.
# =============================================================================


# =============================================================================
# PART 0 -- SETUP
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

dir_proc <- "data/processed"
dir_tab  <- "results/tables"
dir_fig  <- "results/figures"

LFC_CUT  <- 1.5
PADJ_CUT <- 0.05

# Registered prediction (METHODS_LOG §5.13): turning independent filtering
# OFF tests every gene, so BH correction is applied across more hypotheses
# and each surviving gene carries a heavier penalty. If §5.13 is right, the
# DEG count should FALL from 2,252 toward their 1,464.
BASELINE_DEG <- 2252
REF_DEG      <- 1464

set.seed(42)


# =============================================================================
# PART 1 -- LOAD 02 AND 03
# =============================================================================

d2 <- readRDS(file.path(dir_proc, "02_deseq2.rds"))
d3 <- readRDS(file.path(dir_proc, "03_concordance.rds"))

dds  <- d2$dds
meta <- d2$meta

# The Entrez -> Ensembl mapping settled in 03. One row per Ensembl ID and one
# per Entrez ID, so within this subset the mapping is 1:1 and every variant
# below is compared on identical footing.
ens_of      <- setNames(d3$mapped$ensembl, d3$mapped$entrez)
ref_deg_ens <- d3$ref_deg$ensembl

cat("Loaded. Baseline:", length(d3$recovered_ens), "recovered of",
    length(ref_deg_ens), sprintf("(%.1f%%)\n\n", d3$stats$rate_all))


# =============================================================================
# PART 2 -- A REUSABLE SCORER
# =============================================================================

# Every variant below produces a DESeqResults object and gets scored the same
# way. Writing this once means the variants cannot drift apart through
# copy-paste edits -- a real risk when the whole point is comparing them.

score_variant <- function(res, label) {

  # which() not bare indexing -- padj holds NAs (PREP 7.2)
  deg_entrez <- rownames(res)[which(abs(res$log2FoldChange) > LFC_CUT &
                                    res$padj < PADJ_CUT)]

  n_up   <- sum(res$log2FoldChange[which(abs(res$log2FoldChange) > LFC_CUT &
                                         res$padj < PADJ_CUT)] > 0)
  n_down <- length(deg_entrez) - n_up

  deg_ens <- unname(ens_of[deg_entrez])
  deg_ens <- deg_ens[!is.na(deg_ens)]

  recovered <- intersect(ref_deg_ens, deg_ens)

  data.frame(
    variant     = label,
    n_deg       = length(deg_entrez),
    n_up        = n_up,
    n_down      = n_down,
    n_deg_mapped = length(deg_ens),
    n_recovered = length(recovered),
    recovery    = round(100 * length(recovered) / length(ref_deg_ens), 1),
    n_na_padj   = sum(is.na(res$padj)),
    stringsAsFactors = FALSE
  )
}

# Baseline, rescored through the same function so the comparison is exact.
res_base <- d2$res_raw
summ <- score_variant(res_base, "baseline (filtering ON)")


# =============================================================================
# PART 3 -- VARIANT A: independentFiltering = FALSE  (tests §5.13)
# =============================================================================

# No need to rerun DESeq(). Filtering and Cook's are applied at the results()
# step, not during model fitting -- the Wald tests are already done.

res_nofilt <- results(
  dds,
  name                 = "condition_KO_vs_IC1",
  alpha                = PADJ_CUT,
  independentFiltering = FALSE,
  cooksCutoff          = FALSE
)

summ <- rbind(summ, score_variant(res_nofilt, "A: filtering OFF"))

cat("[A] independentFiltering = FALSE, cooksCutoff = FALSE\n")
cat("  padj = NA:", sum(is.na(res_nofilt$padj)),
    "(baseline had", sum(is.na(res_base$padj)), "-- should now be ~0,\n")
cat("      matching their table's zero NAs, §5.6)\n")
cat("  DEGs:", summ$n_deg[2], "vs baseline", summ$n_deg[1],
    "vs reference", REF_DEG, "\n\n")

delta <- summ$n_deg[2] - summ$n_deg[1]

cat("  §5.13 PREDICTION: DEG count should FALL toward 1,464.\n")
if (delta < 0) {
  cat(sprintf("  -> Fell by %d (%.1f%% of the baseline-to-reference gap of %d).\n",
              -delta, 100 * -delta / (BASELINE_DEG - REF_DEG),
              BASELINE_DEG - REF_DEG))
  cat("     Direction supports §5.13. Judge the MAGNITUDE, not just the sign:\n")
  cat("     closing most of the gap is strong support; closing a tenth of it\n")
  cat("     means filtering is one contributor among several.\n\n")
} else {
  cat(sprintf("  -> ROSE by %d. §5.13 is NOT supported. Rewrite it as a\n", delta))
  cat("     rejected explanation and state what is left unexplained.\n\n")
}


# =============================================================================
# PART 4 -- VARIANT B: leave-one-out, drop KO_4
# =============================================================================

# KO_4 sits ~28 units from the other KO samples on PC2 (§9). If the DEG list
# barely moves without it, that is a clean robustness statement. If it moves
# a lot, KO_4 is driving results and that matters more than anything else here.

dds_loo <- dds[, colnames(dds) != meta$gsm[meta$label == "KO_4"]]
dds_loo$condition <- droplevels(dds_loo$condition)

# droplevels() removes factor levels no longer present. Not strictly needed
# here (both conditions survive) but forgetting it after subsetting is a
# common source of "level has no samples" errors.

cat("[B] Leave-one-out: dropping KO_4\n")
cat("  design now:", sum(dds_loo$condition == "IC1"), "IC1 vs",
    sum(dds_loo$condition == "KO"), "KO\n")

dds_loo <- DESeq(dds_loo, quiet = TRUE)
res_loo <- results(dds_loo, name = "condition_KO_vs_IC1", alpha = PADJ_CUT)

summ <- rbind(summ, score_variant(res_loo, "B: drop KO_4"))

# Jaccard index: |intersection| / |union|. 1 means identical sets, 0 disjoint.
deg_base <- rownames(res_base)[which(abs(res_base$log2FoldChange) > LFC_CUT &
                                     res_base$padj < PADJ_CUT)]
deg_loo  <- rownames(res_loo)[which(abs(res_loo$log2FoldChange) > LFC_CUT &
                                    res_loo$padj < PADJ_CUT)]

jac <- length(intersect(deg_base, deg_loo)) / length(union(deg_base, deg_loo))

cat("  DEGs:", length(deg_loo), "vs baseline", length(deg_base), "\n")
cat(sprintf("  Jaccard overlap with baseline DEG list: %.3f\n", jac))
cat("  Above ~0.8 is reassuring; below ~0.6 means KO_4 is load-bearing\n")
cat("  and the writeup must say so.\n\n")


# =============================================================================
# PART 5 -- VARIANT C: continuous passage  (exploratory, §5.10)
# =============================================================================

# §5.10: passage as a FACTOR spends 3 df and leaves 1 residual df from 6
# samples. Continuous passage spends 1 and leaves 3. The PCA does not support
# a shared additive passage effect across lines, so this is exploratory --
# reported, not adopted as a corrected model.

dds_pass <- dds
dds_pass$passage_num <- as.numeric(as.character(dds_pass$passage))

# as.numeric() on a FACTOR returns the level CODES (1,2,3,4), not the labels.
# as.character() first, then as.numeric(). This silently corrupts analyses
# and is one of R's sharpest edges (PREP 3).

cat("[C] Continuous passage. Values:",
    paste(dds_pass$passage_num, collapse = ", "), "\n")
stopifnot(all(dds_pass$passage_num %in% 4:7))

design(dds_pass) <- ~ passage_num + condition
dds_pass <- DESeq(dds_pass, quiet = TRUE)

cat("  resultsNames:", paste(resultsNames(dds_pass), collapse = ", "), "\n")

res_pass <- results(dds_pass, name = "condition_KO_vs_IC1", alpha = PADJ_CUT)
summ <- rbind(summ, score_variant(res_pass, "C: + continuous passage"))

deg_pass <- rownames(res_pass)[which(abs(res_pass$log2FoldChange) > LFC_CUT &
                                     res_pass$padj < PADJ_CUT)]
jac_pass <- length(intersect(deg_base, deg_pass)) / length(union(deg_base, deg_pass))

cat("  DEGs:", length(deg_pass), "vs baseline", length(deg_base), "\n")
cat(sprintf("  Jaccard overlap with baseline: %.3f\n", jac_pass))
cat("  A high overlap means passage absorbs little -- consistent with the\n")
cat("  PCA showing no passage gradient (§5.10).\n\n")


# =============================================================================
# PART 6 -- SUMMARY
# =============================================================================

cat("=============================================================\n")
print(summ, row.names = FALSE)
cat("=============================================================\n")
cat("reference (Table S3):", REF_DEG, "DEGs\n\n")

write.csv(summ, file.path(dir_tab, "05_sensitivity_summary.csv"), row.names = FALSE)


# =============================================================================
# PART 7 -- FIGURES
# =============================================================================

png(file.path(dir_fig, "05_deg_counts.png"), width = 1400, height = 1000, res = 150)
par(mar = c(8, 5, 4, 2))
bp <- barplot(summ$n_deg, names.arg = summ$variant, las = 2,
              col = "steelblue", border = NA, ylim = c(0, max(summ$n_deg) * 1.2),
              ylab = "DEGs at |log2FC| > 1.5, padj < 0.05",
              main = "DEG count across sensitivity variants")
text(bp, summ$n_deg + max(summ$n_deg) * 0.04, labels = summ$n_deg, cex = 0.9)
abline(h = REF_DEG, lty = 2, col = "firebrick", lwd = 2)
text(x = max(bp), y = REF_DEG, labels = "Table S3 (1,464)",
     pos = 3, col = "firebrick", cex = 0.85)
par(mar = c(5, 4, 4, 2) + 0.1)
dev.off()

# par(mar = ...) widens the bottom margin so rotated labels fit. Reset it
# afterwards or every later plot inherits the change.

png(file.path(dir_fig, "05_recovery.png"), width = 1400, height = 1000, res = 150)
par(mar = c(8, 5, 4, 2))
bp <- barplot(summ$recovery, names.arg = summ$variant, las = 2,
              col = "darkseagreen4", border = NA, ylim = c(0, 100),
              ylab = "Recovery of reference DEGs (%)",
              main = "Recovery across sensitivity variants")
text(bp, summ$recovery + 4, labels = sprintf("%.1f%%", summ$recovery), cex = 0.9)
par(mar = c(5, 4, 4, 2) + 0.1)
dev.off()

cat("Figures written to", dir_fig, "\n\n")


# =============================================================================
# PART 8 -- SAVE
# =============================================================================

saveRDS(
  list(
    summary    = summ,
    res_nofilt = res_nofilt,
    res_loo    = res_loo,
    res_pass   = res_pass,
    jaccard    = list(loo = jac, passage = jac_pass),
    provenance = list(run_at = Sys.time(), session = utils::sessionInfo())
  ),
  file = file.path(dir_proc, "05_sensitivity.rds")
)

cat("Saved:", file.path(dir_proc, "05_sensitivity.rds"), "\n")
cat("Next: scripts/04_enrichment.R\n")


# =============================================================================
# NUMBERS FOR METHODS_LOG §9 EXTENSIONS
# =============================================================================
#   A: DEG count with filtering OFF, and how much of the 788-gene gap between
#      baseline (2,252) and reference (1,464) it closes. This is the §5.13
#      verdict -- update §5.13 with the outcome either way.
#   B: DEG count without KO_4, and the Jaccard overlap.
#   C: DEG count with continuous passage, and the Jaccard overlap.


# =============================================================================
# SELF-CHECK
# =============================================================================
#
# 1. Why does variant A not need DESeq() rerun, while B and C do?
# 2. as.numeric(as.character(x)) on a factor -- what breaks if you drop the
#    as.character() step?
# 3. Variant A closes only a small part of the gap. Does that falsify §5.13?
# 4. Why score every variant through one function instead of writing the
#    comparison inline three times?
#
# Answers: 1) Independent filtering and Cook's are applied when results() is
#             called; the Wald tests were already computed. B and C change
#             the samples or the design, so the model must be refit.
#          2) as.numeric() on a factor returns the internal level CODES
#             (1,2,3,4), not the labels (4,5,6,7). Silent corruption.
#          3) No -- it weakens it. §5.13 claims filtering is the leading
#             explanation, not the only one. Report the fraction of the gap
#             closed and name what is left over.
#          4) Copy-paste variants drift. One scorer guarantees the four rows
#             are computed identically, which is the entire point.
# =============================================================================
