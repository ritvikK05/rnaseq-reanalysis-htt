# =============================================================================
# 02_deseq2.R -- GSE270472 reanalysis, KO-NSC vs IC1-NSC
#
# Differential expression. Reads data/processed/01_loaded.rds, builds the
# DESeqDataSet, runs the model, applies shrinkage, and writes the results
# table that 03_concordance.R compares against Supplementary Table 3.
#
# Run from the repo root.
# Writes: data/processed/02_deseq2.rds
#         results/tables/02_results_full.csv
#         results/figures/02_*.png
#
# Cross-references like (PREP 7) point at notes/prep_python_to_r.R.
# =============================================================================


# =============================================================================
# PART 0 -- SETUP
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

# suppressPackageStartupMessages() hides the load banners. Bioconductor
# packages print a lot on attach and it buries your actual output.

dir_proc <- "data/processed"
dir_tab  <- "results/tables"
dir_fig  <- "results/figures"

for (d in c(dir_tab, dir_fig)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Thresholds -- match the paper (METHODS_LOG §4)
LFC_CUT  <- 1.5
PADJ_CUT <- 0.05

# Expectations registered before running (METHODS_LOG §7)
HTT_ENTREZ        <- "3064"
HTT_EXPECTED_L2FC <- 1.34      # Table S3, ENSG00000197386

set.seed(42)


# =============================================================================
# PART 1 -- LOAD 01 OUTPUT
# =============================================================================

loaded <- readRDS(file.path(dir_proc, "01_loaded.rds"))

counts    <- loaded$counts
meta      <- loaded$meta
annot     <- loaded$annot
symbol_of <- loaded$symbol_of

cat("Loaded from 01:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Contrast:", loaded$provenance$contrast, "\n\n")

# Re-assert the invariant. 01 checked it, but 02 must not TRUST 01 -- an RDS
# from an older run would load silently.
stopifnot(all(colnames(counts) == rownames(meta)))
stopifnot(levels(meta$condition)[1] == "IC1")


# =============================================================================
# PART 2 -- BUILD THE DESeqDataSet
# =============================================================================

# A DESeqDataSet is an S4 object (PREP 9) that bundles counts, metadata and
# the design formula together. You do not reach into it with $ -- you use
# accessor functions: counts(dds), colData(dds), sizeFactors(dds).

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition
)

# `~ condition` is a FORMULA object, R's way of writing a model. The left side
# is empty because the response is the count matrix itself. Everything named
# on the right must be a column of colData.

cat("DESeqDataSet built. Design:", deparse(design(dds)), "\n")
cat("Genes before prefiltering:", nrow(dds), "\n")


# --- 2.1 Prefilter -----------------------------------------------------------
# METHODS_LOG §6: rowSums(counts) >= 10. This is NOT a statistical filter --
# DESeq2 does its own independent filtering later. This just drops genes with
# essentially no information so the dispersion fit is not dominated by zeros.

keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

cat("Genes after prefiltering:", nrow(dds),
    sprintf("(dropped %d)\n", sum(!keep)))
cat("  -> METHODS_LOG §9: 'Genes tested after pre-filtering'\n")
cat("  their table has 31,264 rows\n\n")

# `dds[keep, ]` subsets rows of an S4 object the same way as a matrix. The
# method is defined by DESeq2; the syntax is the syntax you already know.


# =============================================================================
# PART 3 -- RUN THE MODEL
# =============================================================================

# DESeq() is three steps in one call:
#   estimateSizeFactors()   -- median-of-ratios normalisation
#   estimateDispersions()   -- gene-wise, then fitted trend, then shrunk
#   nbinomWaldTest()        -- the actual test

dds <- DESeq(dds)

cat("\n")


# --- 3.1 Size factors --------------------------------------------------------
# Sanity: these should track library size but NOT equal it. If they were
# proportional to library size you would just be doing CPM.

sf <- sizeFactors(dds)
lib <- colSums(counts(dds))

print(data.frame(
  sample      = meta$label,
  condition   = as.character(meta$condition),
  lib_M       = round(lib / 1e6, 1),
  size_factor = round(sf, 3),
  ratio       = round(sf / (lib / mean(lib)), 3)
), row.names = FALSE)

cat("\n  size factors span", round(min(sf), 2), "to", round(max(sf), 2), "\n")
cat("  the 'ratio' column is size factor vs pure library-size scaling.\n")
cat("  Departures from 1 are composition effects -- exactly what\n")
cat("  median-of-ratios exists to correct and CPM does not.\n\n")


# --- 3.2 Verify the coefficient name -----------------------------------------
# THE check. If this is not condition_KO_vs_IC1, every sign below is flipped.

cat("resultsNames(dds):\n")
print(resultsNames(dds))

stopifnot("condition_KO_vs_IC1" %in% resultsNames(dds))
cat("  -> KO vs IC1, IC1 as reference. OK\n\n")


# --- GOTCHA: n=2 in the control group ----------------------------------------
# IC1 has two replicates (METHODS_LOG §5.1). Consequences:
#   - dispersion estimates for IC1 lean heavily on the fitted trend rather
#     than on the gene's own data
#   - DESeq2 does not replace Cook's-distance outliers below 7 samples
#     (minReplicatesForReplace = 7), so no outlier replacement happens here
# This is a real power limitation, not an error. Expect to miss borderline
# genes -- which is §5.1's prediction.


# =============================================================================
# PART 4 -- EXTRACT RESULTS
# =============================================================================

# alpha must match the padj cutoff you will actually use. DESeq2's independent
# filtering optimises the filter threshold to maximise discoveries AT alpha --
# leave it at the 0.1 default while testing at 0.05 and you lose genes for no
# reason. Easy to miss; costs real discoveries.

res_raw <- results(dds, name = "condition_KO_vs_IC1", alpha = PADJ_CUT)

cat("Unshrunken results:\n")
summary(res_raw)


# --- 4.1 Shrinkage -----------------------------------------------------------
# apeglm shrinks log2 fold changes toward zero in proportion to their
# uncertainty. A gene at 3 counts vs 0 counts has a huge point estimate and
# almost no evidence; shrinkage pulls it back. p-values are NOT changed --
# passing res = res_raw carries them through untouched.

res_shr <- lfcShrink(
  dds,
  coef = "condition_KO_vs_IC1",
  type = "apeglm",
  res  = res_raw
)

cat("\nShrinkage applied (apeglm).\n")
cat("  max |log2FC| unshrunken:", round(max(abs(res_raw$log2FoldChange), na.rm = TRUE), 2), "\n")
cat("  max |log2FC| shrunken:  ", round(max(abs(res_shr$log2FoldChange), na.rm = TRUE), 2), "\n\n")


# --- 4.2 WHICH log2FC gets thresholded? --------------------------------------
#
# DECISION, not yet in METHODS_LOG §6 -- add it.
#
# The paper's Table S3 contains UNSHRUNKEN estimates (§5.7). Applying the
# |log2FC| > 1.5 cutoff to shrunken values while they applied it to
# unshrunken values would confound two different things: it would measure
# shrinkage, not concordance.
#
# So:
#   res_raw  -> used for the matched-threshold DEG set (apples to apples)
#   res_shr  -> used for ranking, plots, and effect-size reporting
#
# Both are saved. 03_concordance.R uses the unshrunken set as primary and
# reports the shrunken set as a sensitivity analysis.


# =============================================================================
# PART 5 -- CALL DEGs
# =============================================================================

# --- GOTCHA: padj contains NAs (PREP 7.2, METHODS_LOG §6) --------------------
# A bare boolean index on a vector containing NA returns NA ROWS -- silently,
# padded into your result. which() drops NAs instead of propagating them.
#
#   res[res$padj < 0.05, ]          WRONG -- NA rows appear
#   res[which(res$padj < 0.05), ]   RIGHT

call_degs <- function(res, lfc = LFC_CUT, padj = PADJ_CUT) {
  idx <- which(abs(res$log2FoldChange) > lfc & res$padj < padj)
  res[idx, ]
}

deg_raw <- call_degs(res_raw)
deg_shr <- call_degs(res_shr)

n_na <- sum(is.na(res_raw$padj))

cat("NA padj (independent filtering + zero counts):", n_na,
    sprintf("(%.1f%% of tested genes)\n", 100 * n_na / nrow(res_raw)))
cat("  -> METHODS_LOG §5.6: their run had ZERO NAs. These genes are\n")
cat("     unrecoverable by construction in the concordance step.\n\n")

report_degs <- function(deg, label) {
  up   <- sum(deg$log2FoldChange > 0)
  down <- sum(deg$log2FoldChange < 0)
  cat(sprintf("%-28s %5d DEGs  (%d up / %d down)\n", label, nrow(deg), up, down))
}

cat("At |log2FC| >", LFC_CUT, "and padj <", PADJ_CUT, ":\n")
report_degs(deg_raw, "  unshrunken (matched):")
report_degs(deg_shr, "  shrunken (apeglm):")
cat("  paper, recomputed from S3:  1464 DEGs  (1290 up / 174 down)\n\n")

cat("  -> METHODS_LOG §9: 'My DEG count at matched thresholds'\n")
cat("     Use the UNSHRUNKEN row. The gap between the two rows is the\n")
cat("     magnitude of §5.7 and is itself a result worth reporting.\n\n")


# =============================================================================
# PART 6 -- CONTROL GENES
# =============================================================================

# --- 6.1 PRIMARY CONTROL: HTT ------------------------------------------------

cat("[6.1] PRIMARY CONTROL -- HTT (Entrez", HTT_ENTREZ, ")\n")

if (!(HTT_ENTREZ %in% rownames(res_raw))) {
  cat("  *** HTT was removed by prefiltering. Something is very wrong. ***\n")
} else {
  cat(sprintf("  unshrunken log2FC: %+.2f   (Table S3: %+.2f)\n",
              res_raw[HTT_ENTREZ, "log2FoldChange"], HTT_EXPECTED_L2FC))
  cat(sprintf("  shrunken log2FC:   %+.2f\n",
              res_shr[HTT_ENTREZ, "log2FoldChange"]))
  cat(sprintf("  padj:              %.3g\n", res_raw[HTT_ENTREZ, "padj"]))
  cat(sprintf("  baseMean:          %.0f\n", res_raw[HTT_ENTREZ, "baseMean"]))

  if (res_raw[HTT_ENTREZ, "log2FoldChange"] > 0) {
    cat("  -> positive, matching Table S3 and the CPM check in 01. OK\n\n")
  } else {
    cat("  -> *** FAIL: sign disagrees with Table S3. STOP. ***\n\n")
  }
}


# --- 6.2 SECONDARY CONTROLS --------------------------------------------------

sec_up   <- c("TWIST1", "SIX1", "TBX1", "TBX15", "MSX2", "MEOX2", "FOXD1")
sec_down <- c("PAX6")

sec <- data.frame(
  symbol   = c(sec_up, sec_down),
  expected = c(rep("up", length(sec_up)), rep("down", length(sec_down)))
)
sec$entrez <- annot$GeneID[match(sec$symbol, annot$Symbol)]

found <- !is.na(sec$entrez) & sec$entrez %in% rownames(res_raw)

sec$baseMean <- sec$l2fc_raw <- sec$l2fc_shr <- sec$padj <- NA
sec$baseMean[found] <- round(res_raw[sec$entrez[found], "baseMean"], 0)
sec$l2fc_raw[found] <- round(res_raw[sec$entrez[found], "log2FoldChange"], 2)
sec$l2fc_shr[found] <- round(res_shr[sec$entrez[found], "log2FoldChange"], 2)
sec$padj[found]     <- signif(res_raw[sec$entrez[found], "padj"], 3)

sec$ok <- ifelse(is.na(sec$l2fc_raw), NA,
                 ifelse(sec$expected == "up", sec$l2fc_raw > 0, sec$l2fc_raw < 0))

cat("[6.2] SECONDARY CONTROLS\n")
print(sec[, c("symbol", "baseMean", "l2fc_raw", "l2fc_shr", "padj", "expected", "ok")],
      row.names = FALSE)
cat("\n  correct direction:", sum(sec$ok, na.rm = TRUE), "of", sum(!is.na(sec$ok)), "\n")
cat("  PAX6 in the paper: -2.38\n")
cat("  Compare l2fc_raw against l2fc_shr for the low-baseMean genes --\n")
cat("  that gap is §5.12's prediction becoming visible.\n\n")


# =============================================================================
# PART 7 -- DIAGNOSTIC PLOTS
# =============================================================================

# Base graphics only. png() opens a device, you draw, dev.off() closes and
# writes the file. Forget dev.off() and the file stays empty and locked.

# --- 7.1 Dispersion fit ------------------------------------------------------
png(file.path(dir_fig, "02_dispersion.png"), width = 1400, height = 1100, res = 150)
plotDispEsts(dds, main = "Dispersion estimates (KO vs IC1, n=4 vs 2)")
dev.off()

# Read it: black = per-gene estimates, red = fitted trend, blue = final
# shrunken values. Black points should scatter around red. With n=2 controls
# expect heavy shrinkage toward the trend -- most blue sitting on red.

# --- 7.2 MA plots, before and after shrinkage --------------------------------
png(file.path(dir_fig, "02_ma_comparison.png"), width = 2000, height = 900, res = 150)
par(mfrow = c(1, 2))
DESeq2::plotMA(res_raw, ylim = c(-12, 12), main = "Unshrunken")
DESeq2::plotMA(res_shr, ylim = c(-12, 12), main = "apeglm shrunken")
par(mfrow = c(1, 1))
dev.off()

# par(mfrow = c(1, 2)) splits the device into 1 row x 2 columns. Reset it
# afterwards or every later plot inherits the split.
# The left panel's low-baseMean fan is what shrinkage removes.

# --- 7.3 PCA -----------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)

# blind = TRUE ignores the design -- correct for QC, because you want to see
# whether the groups separate WITHOUT having told the transformation about them.

pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))

png(file.path(dir_fig, "02_pca.png"), width = 1200, height = 1100, res = 150)
plot(pca$PC1, pca$PC2,
     col  = ifelse(pca$condition == "KO", "firebrick", "steelblue"),
     pch  = 19, cex = 2,
     xlab = paste0("PC1: ", pct[1], "% variance"),
     ylab = paste0("PC2: ", pct[2], "% variance"),
     main = "VST PCA")
text(pca$PC1, pca$PC2, labels = meta$label, pos = 3, cex = 0.8)
legend("topright", legend = c("IC1", "KO"),
       col = c("steelblue", "firebrick"), pch = 19, bty = "n")
dev.off()

cat("PC1 variance explained:", pct[1], "%\n")
cat("  If PC1 does not separate KO from IC1, stop and reconsider before 03.\n")
cat("  Also check whether passage orders along PC2 -- that is §5.10.\n\n")


# =============================================================================
# PART 8 -- ASSEMBLE AND SAVE
# =============================================================================

# Build one tidy results table with both fold changes and gene symbols.
# as.data.frame() on a DESeqResults object gives a plain data.frame.

out <- as.data.frame(res_raw)
out$entrez     <- rownames(out)
out$symbol     <- symbol_of[out$entrez]
out$ensembl    <- annot$EnsemblGeneID[match(out$entrez, annot$GeneID)]
out$gene_type  <- annot$GeneType[match(out$entrez, annot$GeneID)]
out$l2fc_shrunk <- res_shr$log2FoldChange[match(out$entrez, rownames(res_shr))]

out <- out[, c("entrez", "symbol", "ensembl", "gene_type", "baseMean",
               "log2FoldChange", "l2fc_shrunk", "lfcSE", "pvalue", "padj")]
names(out)[names(out) == "log2FoldChange"] <- "l2fc_raw"

out <- out[order(out$padj, na.last = TRUE), ]

write.csv(out, file.path(dir_tab, "02_results_full.csv"), row.names = FALSE)

saveRDS(
  list(
    dds        = dds,
    res_raw    = res_raw,
    res_shr    = res_shr,
    results_df = out,
    deg_raw    = rownames(deg_raw),
    deg_shr    = rownames(deg_shr),
    meta       = meta,
    annot      = annot,
    symbol_of  = symbol_of,
    params     = list(lfc_cut = LFC_CUT, padj_cut = PADJ_CUT,
                      design = "~ condition", n_prefiltered = nrow(dds)),
    provenance = list(run_at  = Sys.time(),
                      session = utils::sessionInfo())
  ),
  file = file.path(dir_proc, "02_deseq2.rds")
)

cat("Saved:", file.path(dir_proc, "02_deseq2.rds"), "\n")
cat("Saved:", file.path(dir_tab, "02_results_full.csv"), "\n")
cat("Figures in:", dir_fig, "\n")
cat("Next: scripts/03_concordance.R\n")


# =============================================================================
# NUMBERS TO WRITE INTO METHODS_LOG §9 BEFORE MOVING ON
# =============================================================================
#   - Genes tested after prefiltering  (PART 2.1)
#   - DEG count, unshrunken, up/down   (PART 5)
#   - DEG count, shrunken, up/down     (PART 5, the §5.7 magnitude)
#   - NA padj count                    (PART 5, the §5.6 ceiling)
#   - HTT log2FC and padj              (PART 6.1)
#   - PC1 variance explained           (PART 7.3)
# And add the PART 4.2 thresholding decision to §6.


# =============================================================================
# SELF-CHECK
# =============================================================================
#
# 1. Why does results() take alpha = 0.05 rather than being left at default?
# 2. lfcShrink changed log2FoldChange but not padj. Why is that correct?
# 3. Why call DEGs on res_raw rather than res_shr when comparing to the paper?
# 4. res$padj has NAs. What does res[res$padj < 0.05, ] return, and why is
#    which() the fix?
# 5. vst(dds, blind = TRUE) -- what would blind = FALSE change, and why is
#    TRUE right for a QC plot?
#
# Answers: 1) alpha sets the target for independent filtering. Filtering
#             optimised for 0.1 while you test at 0.05 loses discoveries.
#          2) Shrinkage is a change to the effect-size ESTIMATE, not to the
#             evidence against the null. The Wald test already happened.
#          3) Table S3 is unshrunken. Thresholding shrunken values against
#             their unshrunken ones would measure shrinkage, not concordance.
#          4) Rows of NA, silently padded into the result. which() returns
#             only positions where the condition is TRUE, dropping NAs.
#          5) blind = FALSE uses the design when estimating dispersion for the
#             transformation. For QC you want the groups NOT baked in, so that
#             separation in the plot is evidence rather than assumption.
# =============================================================================
