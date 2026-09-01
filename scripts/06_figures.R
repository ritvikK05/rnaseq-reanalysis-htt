# =============================================================================
# 06_figures.R -- GSE270472 reanalysis, KO-NSC vs IC1-NSC
#
# Publication figures for the README. The diagnostic plots from 02-05 stay
# where they are; these two are the ones that carry the argument.
#
#   06_volcano.png          -- DE result with the validated control genes named
#   06_concordance_panel.png -- where the 1,464 reference DEGs went, and why
#
# Run from the repo root, after 05.
# =============================================================================


# =============================================================================
# PART 0 -- SETUP
# =============================================================================

dir_proc <- "data/processed"
dir_fig  <- "results/figures"

LFC_CUT  <- 1.5
PADJ_CUT <- 0.05

d2 <- readRDS(file.path(dir_proc, "02_deseq2.rds"))
d3 <- readRDS(file.path(dir_proc, "03_concordance.rds"))

res       <- as.data.frame(d2$res_raw)
symbol_of <- d2$symbol_of

res$symbol <- symbol_of[rownames(res)]

set.seed(42)


# =============================================================================
# PART 1 -- VOLCANO
# =============================================================================

# Genes with padj = NA were never tested to completion; drop rather than plot
# them at an arbitrary height.
v <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]

# -log10(0) is Inf. Cap at the largest finite value so the axis stays usable.
v$neglog <- -log10(v$padj)
finite_max <- max(v$neglog[is.finite(v$neglog)])
v$neglog[!is.finite(v$neglog)] <- finite_max

v$sig <- with(v, ifelse(padj < PADJ_CUT & log2FoldChange >  LFC_CUT, "up",
                 ifelse(padj < PADJ_CUT & log2FoldChange < -LFC_CUT, "down",
                        "ns")))

cols <- c(up = "firebrick", down = "steelblue", ns = "grey75")

# Controls from METHODS_LOG §7, all RT-qPCR validated by the authors
labels_want <- c("HTT", "PAX6", "TWIST1", "SIX1", "TBX1",
                 "TBX15", "MSX2", "MEOX2", "FOXD1")
lab <- v[v$symbol %in% labels_want, ]

png(file.path(dir_fig, "06_volcano.png"), width = 1500, height = 1300, res = 150)

plot(v$log2FoldChange, v$neglog,
     col  = cols[v$sig], pch = 19, cex = 0.35,
     xlim = c(-12, 12),
     xlab = expression(log[2] * " fold change (KO vs IC1, unshrunken)"),
     ylab = expression(-log[10] * " adjusted p"),
     main = "Differential expression, HTT knockout vs isogenic control")

abline(v = c(-LFC_CUT, LFC_CUT), lty = 2, col = "grey40")
abline(h = -log10(PADJ_CUT),     lty = 2, col = "grey40")

points(lab$log2FoldChange, lab$neglog, pch = 21, bg = "gold", cex = 1.4)
text(lab$log2FoldChange, lab$neglog, labels = lab$symbol,
     pos = 4, cex = 0.75, font = 3)

legend("topleft", bty = "n", cex = 0.85,
       legend = c(sprintf("up (%d)",   sum(v$sig == "up")),
                  sprintf("down (%d)", sum(v$sig == "down")),
                  "not significant",
                  "author-validated controls"),
       col = c("firebrick", "steelblue", "grey75", "black"),
       pt.bg = c(NA, NA, NA, "gold"),
       pch = c(19, 19, 19, 21))

mtext(sprintf("thresholds |log2FC| > %.1f, padj < %.2f (matched to the paper)",
              LFC_CUT, PADJ_CUT),
      side = 3, line = 0.2, cex = 0.75, col = "grey30")

dev.off()

cat("Volcano written. Controls found and labelled:", nrow(lab), "of",
    length(labels_want), "\n")


# =============================================================================
# PART 2 -- CONCORDANCE PANEL
# =============================================================================

st <- d3$stats

# Panel A: what happened to the 1,464 reference DEGs, plus the extras
# Panel B: the 466 misses decomposed by cause

png(file.path(dir_fig, "06_concordance_panel.png"), width = 2000, height = 950, res = 150)
par(mfrow = c(1, 2), mar = c(9, 5, 4, 2))

# --- Panel A -----------------------------------------------------------------
a_vals <- c(st$n_recovered, st$n_missed, st$n_extra)
a_labs <- c("recovered\n(both)", "missed\n(theirs only)", "extra\n(mine only)")
a_cols <- c("firebrick", "darkorange", "steelblue")

bp <- barplot(a_vals, names.arg = a_labs, col = a_cols, border = NA, las = 1,
              ylim = c(0, max(a_vals) * 1.25),
              ylab = "genes",
              main = sprintf("Concordance with Table S3\n%.1f%% of 1,464 recovered",
                             st$rate_all))
text(bp, a_vals + max(a_vals) * 0.05, labels = a_vals, cex = 0.95)

# --- Panel B -----------------------------------------------------------------
# Values from 03 PART 5. Hard-coded because the residual category is a
# subtraction rather than a stored statistic; verify against the console
# output if any upstream parameter changes.
b_vals <- c(244, 10, 72, 32, 108)
b_labs <- c("never tested here\n(unmappable / filtered)",
            "padj = NA",
            "padj 0.05-0.15\n(near miss)",
            "significant but\n|log2FC| < 1.5",
            "genuine\ndisagreement")
b_cols <- c("grey45", "grey65", "darkorange", "orange", "firebrick")

stopifnot(sum(b_vals) == st$n_missed)   # fails loudly if the decomposition drifts

bp2 <- barplot(b_vals, names.arg = b_labs, col = b_cols, border = NA, las = 2,
               cex.names = 0.7, ylim = c(0, max(b_vals) * 1.25),
               ylab = "genes",
               main = sprintf("Why %d reference DEGs were missed", st$n_missed))
text(bp2, b_vals + max(b_vals) * 0.05, labels = b_vals, cex = 0.9)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()

cat("Concordance panel written.\n")


# =============================================================================
# PART 3 -- REPORT THE NUMBERS THE README QUOTES
# =============================================================================

cat("\n=============================================================\n")
cat("NUMBERS QUOTED IN THE README -- verify these match\n")
cat("=============================================================\n")
cat("  genes tested (prefiltered):   ", nrow(res), "\n")
cat("  padj = NA:                    ", sum(is.na(res$padj)), "\n")
cat("  DEGs up / down:               ", sum(v$sig == "up"), "/", sum(v$sig == "down"), "\n")
cat("  reference DEGs:               ", st$n_ref_deg, "\n")
cat("  recovered / missed / extra:   ", st$n_recovered, "/", st$n_missed, "/", st$n_extra, "\n")
cat(sprintf("  recovery all / mappable:       %.1f%% / %.1f%%\n",
            st$rate_all, st$rate_mappable))
cat(sprintf("  sign agreement:                %.1f%%\n", st$pct_sign_agree))
cat(sprintf("  Spearman all / DEGs:           %.3f / %.3f\n", st$rho_all, st$rho_deg))
cat("  HTT log2FC:                   ",
    round(res[rownames(res) == "3064", "log2FoldChange"], 2), "(Table S3: +1.34)\n")
cat("=============================================================\n")
