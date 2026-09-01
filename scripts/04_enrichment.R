# =============================================================================
# 04_enrichment.R -- GSE270472 reanalysis, KO-NSC vs IC1-NSC
#
# GO enrichment with clusterProfiler. Four gene sets:
#   A. upregulated DEGs
#   B. downregulated DEGs
#   C. "extra" DEGs -- called here, not in Table S3. Coherent or noise?
#   D. "recovered" DEGs -- the shared set, as a positive control against the
#      biology the paper reports
#
# Run from the repo root, after 03.
# Writes: data/processed/04_enrichment.rds
#         results/tables/04_go_*.csv
#         results/figures/04_go_*.png
#
# Cross-references like (PREP 7) point at notes/prep_python_to_r.R.
# =============================================================================


# =============================================================================
# PART 0 -- SETUP
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

dir_proc <- "data/processed"
dir_tab  <- "results/tables"
dir_fig  <- "results/figures"

LFC_CUT  <- 1.5
PADJ_CUT <- 0.05
GO_PADJ  <- 0.05
GO_ONT   <- "BP"     # Biological Process. CC and MF exist; BP is what the
                     # paper's neural-crest / mesoderm story lives in.

set.seed(42)


# =============================================================================
# PART 1 -- LOAD AND DEFINE GENE SETS
# =============================================================================

d2 <- readRDS(file.path(dir_proc, "02_deseq2.rds"))
d3 <- readRDS(file.path(dir_proc, "03_concordance.rds"))

res    <- d2$res_raw
mapped <- d3$mapped

# --- 1.1 The universe --------------------------------------------------------
# THE most consequential parameter in an enrichment analysis, and the one most
# often left at a bad default. The universe is the set of genes that COULD
# have been called. Using "all genes in org.Hs.eg.db" instead of "genes I
# actually tested" inflates every p-value, because genes that were never
# testable here cannot be enriched or depleted.
#
# METHODS_LOG §5.11: their PANTHER run used a different universe again, so GO
# results differ for this reason alone before any biology enters.

universe <- rownames(res)[!is.na(res$padj)]

cat("Universe (tested genes with a padj):", length(universe), "\n")
cat("  of", nrow(res), "prefiltered genes\n\n")


# --- 1.2 The four gene sets --------------------------------------------------

deg_idx  <- which(abs(res$log2FoldChange) > LFC_CUT & res$padj < PADJ_CUT)
deg_up   <- rownames(res)[deg_idx[res$log2FoldChange[deg_idx] > 0]]
deg_down <- rownames(res)[deg_idx[res$log2FoldChange[deg_idx] < 0]]

# 03 worked in Ensembl space; enrichGO here works in Entrez. Translate back
# through `mapped`, which holds both.
entrez_of_ens <- setNames(mapped$entrez, mapped$ensembl)

deg_extra     <- unname(entrez_of_ens[d3$extra_ens])
deg_recovered <- unname(entrez_of_ens[d3$recovered_ens])
deg_extra     <- deg_extra[!is.na(deg_extra)]
deg_recovered <- deg_recovered[!is.na(deg_recovered)]

sets <- list(
  up        = deg_up,
  down      = deg_down,
  extra     = deg_extra,
  recovered = deg_recovered
)

cat("Gene sets:\n")
for (nm in names(sets)) cat(sprintf("  %-10s %5d genes\n", nm, length(sets[[nm]])))
cat("\n")


# =============================================================================
# PART 2 -- RUN ENRICHMENT
# =============================================================================

# enrichGO does a hypergeometric (Fisher) test per GO term: given the universe
# and the size of your gene set, is this term over-represented? Then corrects
# across terms with BH.

run_go <- function(genes, label) {
  cat("Running enrichGO:", label, "...\n")

  ego <- clusterProfiler::enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = GO_ONT,
    pAdjustMethod = "BH",
    pvalueCutoff  = GO_PADJ,
    qvalueCutoff  = 0.2,
    readable      = TRUE          # convert Entrez to symbols in the output
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("  no significant terms\n")
    return(NULL)
  }

  cat("  ", nrow(as.data.frame(ego)), "significant terms\n")
  ego
}

ego_raw <- lapply(sets, function(g) run_go(g, "set"))
names(ego_raw) <- names(sets)

# lapply() applies a function over a list and returns a list -- the vectorised
# alternative to a for loop that accumulates into a growing object (PREP 1.5).

cat("\n")


# --- 2.1 Simplify redundant terms --------------------------------------------
# GO is a hierarchy, so a real signal shows up as a dozen near-identical terms
# at different depths. simplify() collapses terms with semantic similarity
# above the cutoff, keeping the most significant representative.

ego <- lapply(ego_raw, function(e) {
  if (is.null(e)) return(NULL)
  clusterProfiler::simplify(e, cutoff = 0.7, by = "p.adjust", select_fun = min)
})

cat("After simplify(cutoff = 0.7):\n")
for (nm in names(ego)) {
  n <- if (is.null(ego[[nm]])) 0 else nrow(as.data.frame(ego[[nm]]))
  cat(sprintf("  %-10s %4d terms\n", nm, n))
}
cat("\n")


# =============================================================================
# PART 3 -- INSPECT
# =============================================================================

show_top <- function(e, label, n = 12) {
  cat("=============================================================\n")
  cat(label, "\n")
  cat("=============================================================\n")
  if (is.null(e)) { cat("  no significant terms\n\n"); return(invisible(NULL)) }
  df <- as.data.frame(e)
  df <- df[order(df$p.adjust), ]
  print(df[seq_len(min(n, nrow(df))), c("Description", "GeneRatio", "p.adjust")],
        row.names = FALSE)
  cat("\n")
}

show_top(ego$up,        "A. UPREGULATED in KO")
show_top(ego$down,      "B. DOWNREGULATED in KO")
show_top(ego$recovered, "D. RECOVERED (shared with Table S3) -- positive control")
show_top(ego$extra,     "C. EXTRA (called here, not in Table S3)")

cat("READ THESE AGAINST THE PAPER (METHODS_LOG §7):\n")
cat("  Expect UP to hit neural crest / mesoderm / skeletal patterning --\n")
cat("  TWIST1, SIX1, TBX1, TBX15, MSX2, MEOX2, FOXD1 are all in that story.\n")
cat("  Expect DOWN to hit neurogenesis / CNS development -- PAX6 is the\n")
cat("  canonical NSC identity marker and is strongly reduced.\n")
cat("  Set D should look like the paper. If it does not, the concordance\n")
cat("  numbers in 03 are agreeing on gene lists but not on biology.\n\n")

cat("THE QUESTION FOR SET C:\n")
cat("  If the 786 extra DEGs enrich for the SAME themes as the recovered\n")
cat("  set, they are additional signal in the same biology -- consistent\n")
cat("  with a power difference, not with artefact.\n")
cat("  If they enrich for something unrelated, or for nothing at all, that\n")
cat("  points at technical noise and belongs in §5.13.\n\n")


# =============================================================================
# PART 4 -- THEME OVERLAP BETWEEN RECOVERED AND EXTRA
# =============================================================================

# Quantify the set C question rather than eyeballing it.

terms_of <- function(e) {
  if (is.null(e)) return(character(0))
  as.data.frame(e)$ID
}

t_rec   <- terms_of(ego$recovered)
t_extra <- terms_of(ego$extra)

shared_terms <- intersect(t_rec, t_extra)
jac_terms    <- if (length(union(t_rec, t_extra)) == 0) NA else
                length(shared_terms) / length(union(t_rec, t_extra))

cat("[4] GO term overlap, recovered vs extra\n")
cat("  recovered terms:", length(t_rec), "\n")
cat("  extra terms:    ", length(t_extra), "\n")
cat("  shared:         ", length(shared_terms), "\n")
cat(sprintf("  Jaccard:         %.3f\n", jac_terms))
cat("  High overlap = same biology, more of it. Low = different signal.\n\n")

if (length(shared_terms) > 0) {
  rec_df <- as.data.frame(ego$recovered)
  cat("  Terms found in BOTH (top 10 by recovered-set significance):\n")
  sub <- rec_df[rec_df$ID %in% shared_terms, ]
  sub <- sub[order(sub$p.adjust), ]
  print(head(sub[, c("Description", "p.adjust")], 10), row.names = FALSE)
  cat("\n")
}


# =============================================================================
# PART 5 -- FIGURES
# =============================================================================

# Base graphics. clusterProfiler's own barplot()/dotplot() return ggplot
# objects; drawing these by hand keeps the repo consistent with the rest of
# the scripts and gives control over label truncation.

plot_go <- function(e, label, file, n = 15, col = "steelblue") {
  if (is.null(e)) return(invisible(NULL))
  df <- as.data.frame(e)
  df <- df[order(df$p.adjust), ]
  df <- df[seq_len(min(n, nrow(df))), ]
  df <- df[order(df$p.adjust, decreasing = TRUE), ]   # barplot draws bottom-up

  # Truncate long GO descriptions so labels stay legible
  labs <- ifelse(nchar(df$Description) > 45,
                 paste0(substr(df$Description, 1, 42), "..."),
                 df$Description)

  png(file.path(dir_fig, file), width = 1700, height = 1100, res = 150)
  par(mar = c(5, 22, 4, 2))                # wide left margin for term names
  barplot(-log10(df$p.adjust), names.arg = labs, horiz = TRUE, las = 1,
          col = col, border = NA, cex.names = 0.75,
          xlab = "-log10(adjusted p)", main = label)
  par(mar = c(5, 4, 4, 2) + 0.1)
  dev.off()
}

plot_go(ego$up,        "GO BP: upregulated in KO",   "04_go_up.png",        col = "firebrick")
plot_go(ego$down,      "GO BP: downregulated in KO", "04_go_down.png",      col = "steelblue")
plot_go(ego$recovered, "GO BP: recovered DEGs",      "04_go_recovered.png", col = "darkseagreen4")
plot_go(ego$extra,     "GO BP: extra DEGs",          "04_go_extra.png",     col = "darkorange3")

cat("Figures written to", dir_fig, "\n\n")


# =============================================================================
# PART 6 -- SAVE
# =============================================================================

for (nm in names(ego)) {
  if (is.null(ego[[nm]])) next
  write.csv(as.data.frame(ego[[nm]]),
            file.path(dir_tab, paste0("04_go_", nm, ".csv")),
            row.names = FALSE)
}

saveRDS(
  list(
    ego          = ego,
    ego_raw      = ego_raw,
    sets         = sets,
    universe     = universe,
    shared_terms = shared_terms,
    jac_terms    = jac_terms,
    params       = list(ont = GO_ONT, padj = GO_PADJ, simplify_cutoff = 0.7),
    provenance   = list(run_at = Sys.time(), session = utils::sessionInfo())
  ),
  file = file.path(dir_proc, "04_enrichment.rds")
)

cat("Saved:", file.path(dir_proc, "04_enrichment.rds"), "\n")
cat("Next: figures for the writeup, then README.\n")


# =============================================================================
# NUMBERS FOR METHODS_LOG
# =============================================================================
#   §9: term counts per set; whether up/down themes match the paper's
#       neural-crest and neurogenesis story; GO term Jaccard between
#       recovered and extra (PART 4) and what it implies for §5.13.
#   §5.11: this used the tested-gene universe (n above), BH, GO:BP,
#       simplify(0.7). Their PANTHER run used a different universe and
#       algorithm — record both so the comparison is interpretable.


# =============================================================================
# SELF-CHECK
# =============================================================================
#
# 1. Why is `universe` the tested genes rather than every gene in org.Hs.eg.db?
# 2. What does simplify(cutoff = 0.7) do, and why is the raw output kept too?
# 3. Set C enriches for the same terms as set D. What does that support?
# 4. Set C returns no significant terms at all. What would that support?
# 5. Your GO results differ from the paper's. Name two reasons that have
#    nothing to do with your DEG list being different.
#
# Answers: 1) Enrichment asks whether a term is over-represented among genes
#             that COULD have been called. Genes never tested here cannot be,
#             so including them inflates every p-value.
#          2) Collapses semantically redundant GO terms, keeping the most
#             significant representative. ego_raw is kept because the
#             simplification threshold is a judgement call someone may want
#             to vary.
#          3) That the extra DEGs are more of the same biology, consistent
#             with a power or sensitivity difference rather than artefact.
#          4) That they are probably technical -- scattered genes with no
#             coherent function, which is what noise looks like.
#          5) Different universe, and different tool/algorithm: PANTHER v17
#             Fisher vs clusterProfiler BH (§5.11). Either alone changes
#             which terms clear significance.
# =============================================================================
