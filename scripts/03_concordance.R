# =============================================================================
# 03_concordance.R -- GSE270472 reanalysis, KO-NSC vs IC1-NSC
#
# Compare this reanalysis against Supplementary Table 3. Maps Entrez to
# Ensembl, computes the reference DEG set, and quantifies recovery, sign
# agreement, and effect-size correlation. Tests two registered predictions:
#   §5.12 -- recovery lowest among low-baseMean genes
#   §9    -- the excess downregulated genes are compositionally distinct
#
# Run from the repo root.
# Writes: data/processed/03_concordance.rds
#         results/tables/03_*.csv
#         results/figures/03_*.png
#
# Cross-references like (PREP 7) point at notes/prep_python_to_r.R.
# =============================================================================


# =============================================================================
# PART 0 -- SETUP
# =============================================================================

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("readxl not installed. Run: install.packages('readxl')")
}

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# requireNamespace() tests availability WITHOUT attaching the package. Use it
# for a dependency you will call with :: rather than attach wholesale.

dir_proc <- "data/processed"
dir_tab  <- "results/tables"
dir_fig  <- "results/figures"
dir_raw  <- "data/raw"

for (d in c(dir_tab, dir_fig)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

path_s3 <- file.path(dir_raw, "13578_2025_1443_MOESM5_ESM.xlsx")

LFC_CUT  <- 1.5
PADJ_CUT <- 0.05

# Registered expectation (METHODS_LOG §4)
REF_DEG_EXPECTED <- 1464

set.seed(42)


# =============================================================================
# PART 1 -- LOAD 02 OUTPUT
# =============================================================================

d2 <- readRDS(file.path(dir_proc, "02_deseq2.rds"))

mine  <- d2$results_df        # entrez, symbol, ensembl, gene_type, baseMean,
                              # l2fc_raw, l2fc_shrunk, lfcSE, pvalue, padj
annot <- d2$annot

cat("Loaded from 02:", nrow(mine), "tested genes\n")
cat("My DEGs (unshrunken, matched thresholds):", length(d2$deg_raw), "\n\n")

stopifnot(nrow(mine) > 0, "ensembl" %in% names(mine))


# =============================================================================
# PART 2 -- LOAD THE TRUTH TABLE
# =============================================================================

# METHODS_LOG §3: sheet "KO-NSC vs IC1-NSC" needs skip = 1. The HD sheet needs
# skip = 2 -- the offsets differ between sheets. Do not assume.

ref <- readxl::read_excel(path_s3, sheet = "KO-NSC vs IC1-NSC", skip = 1)
ref <- as.data.frame(ref)

cat("Table S3 loaded:", nrow(ref), "rows (expect 31,264)\n")
cat("Columns:", paste(names(ref), collapse = " | "), "\n")

# --- GOTCHA: column names with hyphens ---------------------------------------
# "P-value" is not a syntactic R name. df$`P-value` needs backticks and breaks
# the moment anything renames it. Normalise once, immediately.

names(ref)[names(ref) == "gene_ID"]        <- "ensembl"
names(ref)[names(ref) == "Gene_name"]      <- "symbol_ref"
names(ref)[names(ref) == "log2FoldChange"] <- "l2fc_ref"
names(ref)[names(ref) == "P_adjusted"]     <- "padj_ref"

stopifnot(all(c("ensembl", "l2fc_ref", "padj_ref") %in% names(ref)))

ref$l2fc_ref <- as.numeric(ref$l2fc_ref)
ref$padj_ref <- as.numeric(ref$padj_ref)

cat("NA padj in their table:", sum(is.na(ref$padj_ref)),
    " (expect 0 -- §5.6, they disabled independent filtering)\n")
cat("Duplicate Ensembl IDs in their table:", sum(duplicated(ref$ensembl)), "\n\n")


# --- 2.1 Reference DEG set ---------------------------------------------------

ref_deg_idx <- which(abs(ref$l2fc_ref) > LFC_CUT & ref$padj_ref < PADJ_CUT)
ref_deg     <- ref[ref_deg_idx, ]

cat("Reference DEGs recomputed:", nrow(ref_deg),
    sprintf("(%d up / %d down)\n", sum(ref_deg$l2fc_ref > 0),
            sum(ref_deg$l2fc_ref < 0)))
cat("  registered expectation:", REF_DEG_EXPECTED, "(1290 up / 174 down)\n")

if (nrow(ref_deg) != REF_DEG_EXPECTED) {
  cat("  *** MISMATCH -- check the skip= offset before continuing ***\n\n")
} else {
  cat("  -> matches §4. OK\n\n")
}


# =============================================================================
# PART 3 -- ENTREZ TO ENSEMBL MAPPING (§5.4)
# =============================================================================

# 02 already carried the GEO annotation's EnsemblGeneID through. Start there,
# then try to recover the gaps from org.Hs.eg.db.

n_have <- sum(!is.na(mine$ensembl))
cat("[3.1] From GEO annotation:", n_have, "of", nrow(mine),
    sprintf("(%.1f%%) have an Ensembl ID\n", 100 * n_have / nrow(mine)))


# --- 3.2 Recover missing IDs from org.Hs.eg.db --------------------------------

missing_entrez <- mine$entrez[is.na(mine$ensembl)]

# AnnotationDbi::select() -- namespace it explicitly. dplyr and others export
# a select() and whichever attached last wins (METHODS_LOG §6).
# It returns a LONG data.frame with one row per mapping, so an Entrez ID with
# three Ensembl IDs produces three rows. Never assume nrow(out) == length(keys).

recovered <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys     = missing_entrez,
    keytype  = "ENTREZID",
    columns  = "ENSEMBL"
  )
)
recovered <- recovered[!is.na(recovered$ENSEMBL), ]

# Collapse 1:many by keeping the first mapping per Entrez ID. Arbitrary, but
# these are genes the primary annotation could not map at all -- mostly ncRNA
# and pseudogenes (§9), so the choice has little leverage on a
# protein-coding-dominated DEG list. Recorded rather than hidden.
recovered <- recovered[!duplicated(recovered$ENTREZID), ]

idx <- match(mine$entrez, recovered$ENTREZID)
mine$ensembl[is.na(mine$ensembl) & !is.na(idx)] <-
  recovered$ENSEMBL[idx[is.na(mine$ensembl) & !is.na(idx)]]

n_after <- sum(!is.na(mine$ensembl))
cat("[3.2] Recovered via org.Hs.eg.db:", n_after - n_have, "additional genes\n")
cat("      Total mapped:", n_after,
    sprintf("(%.1f%%)\n", 100 * n_after / nrow(mine)))
cat("      -> METHODS_LOG §9: 'Additional genes recovered'\n\n")


# --- 3.3 Collapse duplicate Ensembl IDs on my side ----------------------------

# Several Entrez IDs can map to one Ensembl ID. For a gene-level comparison
# each Ensembl ID must appear once. Keep the most significant (lowest padj),
# breaking ties on baseMean -- the entry with the most evidence behind it.

mapped <- mine[!is.na(mine$ensembl), ]
mapped <- mapped[order(mapped$padj, -mapped$baseMean, na.last = TRUE), ]

n_dup <- sum(duplicated(mapped$ensembl))
mapped <- mapped[!duplicated(mapped$ensembl), ]

cat("[3.3] Duplicate Ensembl IDs collapsed:", n_dup,
    "(kept lowest padj, tie-break on baseMean)\n")
cat("      Unique Ensembl IDs on my side:", nrow(mapped), "\n")
cat("      -> METHODS_LOG §9: 'Duplicate resolution method'\n\n")


# --- 3.4 What is unrecoverable by construction --------------------------------

ref_deg_mappable <- ref_deg$ensembl %in% mapped$ensembl
n_unmappable     <- sum(!ref_deg_mappable)

cat("[3.4] Reference DEGs absent from my tested gene universe:", n_unmappable,
    sprintf("(%.1f%% of %d)\n", 100 * n_unmappable / nrow(ref_deg), nrow(ref_deg)))
cat("      These cannot be recovered regardless of statistics -- they were\n")
cat("      never tested here. Causes: §5.2 annotation source, §5.4 mapping\n")
cat("      loss, and prefiltering.\n")
cat("      -> METHODS_LOG §9: 'Reference DEGs unmappable by construction'\n\n")


# =============================================================================
# PART 4 -- CONCORDANCE
# =============================================================================

my_deg_ens  <- mapped$ensembl[which(abs(mapped$l2fc_raw) > LFC_CUT &
                                    mapped$padj < PADJ_CUT)]
ref_deg_ens <- ref_deg$ensembl

recovered_ens <- intersect(ref_deg_ens, my_deg_ens)
missed_ens    <- setdiff(ref_deg_ens, my_deg_ens)
extra_ens     <- setdiff(my_deg_ens, ref_deg_ens)

# intersect / setdiff are R's set operations (PREP 8). setdiff is directional:
# setdiff(A, B) is "in A but not B".

cat("[4.1] HEADLINE CONCORDANCE\n")
cat("  Reference DEGs (theirs):     ", length(ref_deg_ens), "\n")
cat("  My DEGs (mapped, unshrunken):", length(my_deg_ens), "\n")
cat("  Recovered (shared):          ", length(recovered_ens), "\n")
cat("  Missed (theirs, not mine):   ", length(missed_ens), "\n")
cat("  Extra (mine, not theirs):    ", length(extra_ens), "\n\n")

rate_all      <- 100 * length(recovered_ens) / length(ref_deg_ens)
rate_mappable <- 100 * length(recovered_ens) / sum(ref_deg_mappable)

cat(sprintf("  Recovery rate, all reference DEGs:      %.1f%%\n", rate_all))
cat(sprintf("  Recovery rate, mappable reference DEGs: %.1f%%\n", rate_mappable))
cat("  Report BOTH. The first is the honest headline; the second isolates\n")
cat("  statistical disagreement from unmappability (§5.4).\n\n")


# --- 4.2 Sign agreement -------------------------------------------------------

shared <- merge(
  mapped[, c("ensembl", "symbol", "gene_type", "baseMean",
             "l2fc_raw", "l2fc_shrunk", "padj")],
  ref[, c("ensembl", "l2fc_ref", "padj_ref")],
  by = "ensembl"
)

# merge() is an inner join by default -- pandas' pd.merge(how="inner").

cat("[4.2] Genes tested in both analyses:", nrow(shared), "\n")

shared_deg  <- shared[shared$ensembl %in% recovered_ens, ]
sign_agree  <- sign(shared_deg$l2fc_raw) == sign(shared_deg$l2fc_ref)
pct_sign    <- 100 * mean(sign_agree, na.rm = TRUE)

cat(sprintf("  Sign agreement among recovered DEGs: %.1f%%\n", pct_sign))
if (pct_sign < 95) {
  cat("  *** Below 95% -- check reference levels before interpreting. ***\n")
} else {
  cat("  -> as expected. OK\n")
}


# --- 4.3 Effect-size correlation ----------------------------------------------

rho_all <- cor(shared$l2fc_raw, shared$l2fc_ref,
               method = "spearman", use = "complete.obs")
rho_deg <- cor(shared_deg$l2fc_raw, shared_deg$l2fc_ref,
               method = "spearman", use = "complete.obs")

cat(sprintf("\n[4.3] Spearman rho, all shared genes:     %.3f\n", rho_all))
cat(sprintf("      Spearman rho, recovered DEGs only: %.3f\n", rho_deg))
cat("      Spearman (rank-based) not Pearson: their table contains\n")
cat("      near-zero-denominator artifacts like FKBPL at log2FC 22 (§5.7)\n")
cat("      that would dominate a Pearson correlation.\n\n")


# =============================================================================
# PART 5 -- WHY WERE GENES MISSED?
# =============================================================================

missed <- shared[shared$ensembl %in% missed_ens, ]

# Genes in missed_ens but NOT in `shared` were never tested here at all.
n_missed_untested <- length(missed_ens) - nrow(missed)

cat("[5.1] Missed genes:", length(missed_ens), "total\n")
cat("  never tested here (unmappable / prefiltered):", n_missed_untested, "\n")
cat("  tested but not called:", nrow(missed), "\n\n")

n_na    <- sum(is.na(missed$padj))
n_close <- sum(missed$padj >= PADJ_CUT & missed$padj < 0.15, na.rm = TRUE)
n_lfc   <- sum(missed$padj < PADJ_CUT & abs(missed$l2fc_raw) <= LFC_CUT, na.rm = TRUE)

cat("  Of those tested but not called:\n")
cat(sprintf("    padj = NA (independent filtering, §5.6): %4d (%.1f%%)\n",
            n_na, 100 * n_na / nrow(missed)))
cat(sprintf("    padj between 0.05 and 0.15 (near miss):  %4d (%.1f%%)\n",
            n_close, 100 * n_close / nrow(missed)))
cat(sprintf("    significant but |log2FC| below cutoff:   %4d (%.1f%%)\n",
            n_lfc, 100 * n_lfc / nrow(missed)))
cat("  -> METHODS_LOG §9, Concordance block\n\n")


# =============================================================================
# PART 6 -- PREDICTION TEST: §5.12 baseMean stratification
# =============================================================================

# Registered before running 02: recovery should be LOWEST among low-baseMean
# genes, because their DEG list is dominated by genes switching on from near
# zero and those are exactly the ones with the least evidence here.

ref_in_shared <- shared[shared$ensembl %in% ref_deg_ens, ]
ref_in_shared$recovered <- ref_in_shared$ensembl %in% recovered_ens

qs <- quantile(ref_in_shared$baseMean, probs = seq(0, 1, 0.25), na.rm = TRUE)
ref_in_shared$bm_quartile <- cut(ref_in_shared$baseMean, breaks = qs,
                                 include.lowest = TRUE,
                                 labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"))

# cut() bins a continuous variable into a factor -- pandas' pd.cut().

recov_by_q <- tapply(ref_in_shared$recovered, ref_in_shared$bm_quartile,
                     function(x) 100 * mean(x))

# tapply() = split-apply-combine on one grouping factor: df.groupby(g).agg(f)

cat("[6] §5.12 PREDICTION TEST -- recovery by baseMean quartile\n")
print(round(recov_by_q, 1))
cat("\n  Prediction: recovery rises left to right.\n")
if (!any(is.na(recov_by_q)) && recov_by_q[1] < recov_by_q[4]) {
  cat("  -> Q1 below Q4. Direction supports §5.12.\n\n")
} else {
  cat("  -> Direction does NOT support §5.12. Report the null result.\n\n")
}


# =============================================================================
# PART 7 -- PREDICTION TEST: composition of the excess downregulated genes
# =============================================================================

# From §9: my up count is close to theirs (+15%) but my down count is 4.4x
# theirs. If the excess is compositionally distinct -- ncRNA-heavy, say --
# that points at §5.5 / §5.2 rather than at a uniform power difference.

extra <- mapped[mapped$ensembl %in% extra_ens, ]
extra$direction <- ifelse(extra$l2fc_raw > 0, "up", "down")

comp <- table(extra$gene_type, extra$direction)
comp <- comp[order(-rowSums(comp)), , drop = FALSE]

cat("[7] Gene-type composition of DEGs called here but not in Table S3\n")
print(head(comp, 10))

pct_nc_down <- 100 * sum(extra$direction == "down" &
                         grepl("ncRNA|pseudo", extra$gene_type, ignore.case = TRUE)) /
               max(1, sum(extra$direction == "down"))
pct_nc_up   <- 100 * sum(extra$direction == "up" &
                         grepl("ncRNA|pseudo", extra$gene_type, ignore.case = TRUE)) /
               max(1, sum(extra$direction == "up"))

cat(sprintf("\n  ncRNA + pseudogene share of extra DOWN: %.1f%%\n", pct_nc_down))
cat(sprintf("  ncRNA + pseudogene share of extra UP:   %.1f%%\n", pct_nc_up))
cat("  A large gap implicates §5.5 (rRNA depletion changes what dominates\n")
cat("  the library, and therefore size factors) or §5.2 (annotation source).\n")
cat("  Comparable shares mean the asymmetry needs a different explanation.\n\n")


# =============================================================================
# PART 8 -- FIGURES
# =============================================================================

# --- 8.1 log2FC scatter, mine vs theirs --------------------------------------

png(file.path(dir_fig, "03_lfc_scatter.png"), width = 1300, height = 1200, res = 150)

cat_col <- ifelse(shared$ensembl %in% recovered_ens, "firebrick",
           ifelse(shared$ensembl %in% missed_ens,    "darkorange",
           ifelse(shared$ensembl %in% extra_ens,     "steelblue",
                                                     "grey80")))

plot(shared$l2fc_ref, shared$l2fc_raw,
     col  = cat_col, pch = 19, cex = 0.4,
     xlim = c(-10, 10), ylim = c(-10, 10),
     xlab = "log2FC, Table S3 (unshrunken)",
     ylab = "log2FC, this reanalysis (unshrunken)",
     main = sprintf("Effect-size concordance (Spearman rho = %.3f)", rho_all))
abline(0, 1, lty = 2, col = "grey40")
abline(h = 0, v = 0, lty = 3, col = "grey60")
legend("topleft",
       legend = c("recovered", "missed", "extra", "neither"),
       col    = c("firebrick", "darkorange", "steelblue", "grey80"),
       pch = 19, bty = "n", cex = 0.9)
dev.off()

# --- 8.2 Recovery by baseMean quartile ---------------------------------------

png(file.path(dir_fig, "03_recovery_by_basemean.png"), width = 1200, height = 1000, res = 150)
bp <- barplot(recov_by_q, ylim = c(0, 100),
              col = "steelblue", border = NA,
              ylab = "Recovery rate (%)",
              xlab = "baseMean quartile (this reanalysis)",
              main = "Recovery of reference DEGs by expression level")
text(bp, recov_by_q + 4, labels = sprintf("%.0f%%", recov_by_q), cex = 0.9)
dev.off()

cat("Figures written to", dir_fig, "\n\n")


# =============================================================================
# PART 9 -- SAVE
# =============================================================================

write.csv(shared, file.path(dir_tab, "03_shared_genes.csv"), row.names = FALSE)
write.csv(missed, file.path(dir_tab, "03_missed_genes.csv"), row.names = FALSE)
write.csv(extra,  file.path(dir_tab, "03_extra_genes.csv"),  row.names = FALSE)

saveRDS(
  list(
    mapped        = mapped,
    ref           = ref,
    ref_deg       = ref_deg,
    shared        = shared,
    recovered_ens = recovered_ens,
    missed_ens    = missed_ens,
    extra_ens     = extra_ens,
    stats = list(
      n_ref_deg      = length(ref_deg_ens),
      n_my_deg       = length(my_deg_ens),
      n_recovered    = length(recovered_ens),
      n_missed       = length(missed_ens),
      n_extra        = length(extra_ens),
      n_unmappable   = n_unmappable,
      rate_all       = rate_all,
      rate_mappable  = rate_mappable,
      pct_sign_agree = pct_sign,
      rho_all        = rho_all,
      rho_deg        = rho_deg,
      recov_by_q     = recov_by_q
    ),
    provenance = list(run_at = Sys.time(), session = utils::sessionInfo())
  ),
  file = file.path(dir_proc, "03_concordance.rds")
)

cat("Saved:", file.path(dir_proc, "03_concordance.rds"), "\n")
cat("Next: scripts/04_enrichment.R\n")


# =============================================================================
# NUMBERS FOR METHODS_LOG §9
# =============================================================================
#   Mapping block:    additional genes recovered (3.2), duplicate resolution
#                     method (3.3), reference DEGs unmappable (3.4)
#   Concordance:      recovered / missed / extra (4.1), both recovery rates,
#                     sign agreement (4.2), Spearman rho (4.3),
#                     missed-gene breakdown (5.1),
#                     recovery by baseMean quartile (6),
#                     gene-type composition of extras (7)


# =============================================================================
# SELF-CHECK
# =============================================================================
#
# 1. Why report two recovery rates rather than one?
# 2. AnnotationDbi::select() returned more rows than keys. Why, and what
#    would break if you assumed otherwise?
# 3. Why Spearman rather than Pearson for the log2FC correlation?
# 4. A missed gene has padj = NA here. Is that a disagreement about biology?
# 5. §5.12 predicted recovery rises with baseMean. If Q1 recovery came back
#    HIGHER than Q4, what would you write in the log?
#
# Answers: 1) The all-DEG rate is the honest headline; the mappable-only rate
#             isolates statistical disagreement from genes that were never
#             testable here (§5.4). Reporting only the second flatters.
#          2) 1:many Entrez-to-Ensembl mappings produce one row per pair.
#             Assuming 1:1 silently misaligns every downstream match().
#          3) Their table contains near-zero-denominator artifacts at log2FC
#             around 22 (§5.7). Pearson would be dominated by a few such
#             points; rank correlation is not.
#          4) No. It is a gene independent filtering never tested here but
#             which they tested, because they disabled filtering (§5.6).
#             Unrecoverable by construction, not a biological disagreement.
#          5) The null result, plainly, plus what it rules out. A registered
#             prediction that fails is still evidence -- it means the
#             asymmetry needs an explanation other than baseline expression.
# =============================================================================
