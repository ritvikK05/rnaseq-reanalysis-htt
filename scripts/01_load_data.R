# =============================================================================
# 01_load_data.R -- GSE270472 reanalysis, KO-NSC vs IC1-NSC
#
# Import the NCBI raw count matrix + GEO annotation, build the sample table,
# subset to the KO vs IC1 contrast, and run pre-flight sanity checks.
#
# Base R only. DESeq2 gets loaded in 02_deseq2.R.
# Run from the repo root -- check getwd() if the paths fail.
# Writes: data/processed/01_loaded.rds
#
# Cross-references below like (PREP 3) point at notes/prep_python_to_r.R.
# =============================================================================


# =============================================================================
# PART 0 -- SETUP AND PATHS
# =============================================================================

# Every hard-coded string lives here so the script is auditable at a glance.

path_counts <- "data/raw/GSE270472_raw_counts_GRCh38.p13_NCBI.tsv.gz"
path_annot  <- "data/raw/Human.GRCh38.p13.annot.tsv.gz"
dir_out     <- "data/processed"

# IC1 is written FIRST because the first factor level is the reference (PREP 3).
groups_keep <- c("IC1", "KO")

HTT_ENTREZ <- "3064"                    # primary control, METHODS_LOG §7
HTT_EXPECTED_L2FC <- 1.34       # Table S3, ENSG00000197386 -- KO is HIGHER

# Secondary controls -- authors validated these by RT-qPCR
genes_up_expected   <- c("TWIST1", "SIX1", "TBX1", "TBX15", "MSX2", "MEOX2", "FOXD1")
genes_down_expected <- c("PAX6")

set.seed(42)


# =============================================================================
# PART 1 -- LOAD THE COUNT MATRIX
# =============================================================================

# read.delim() reads .gz directly. No decompression step.

counts_raw <- read.delim(
  path_counts,
  header       = TRUE,
  check.names  = FALSE,   # do NOT let R "repair" the GSM column names
  quote        = "",      # see GOTCHA below
  comment.char = ""       # a "#" in a field would otherwise truncate the line
)

# --- GOTCHA: quote = "" is not optional on GEO TSVs --------------------------
# Gene descriptions contain apostrophes ("5' end of..."). R's default quoting
# treats that apostrophe as an opening string delimiter and swallows every row
# until it finds a closing one. You get a truncated table and NO error.
# If a row count ever comes back wrong, check this argument first.

cat("Counts file:", nrow(counts_raw), "rows x", ncol(counts_raw), "cols\n")
str(counts_raw[, 1:4])                  # your df.info() -- use it constantly

# --- 1.1 data.frame -> matrix ------------------------------------------------
# A data.frame is a list of columns that can each have a different type (pandas
# DataFrame). A matrix is one homogeneous block. DESeq2 wants a MATRIX of
# integers with gene IDs as rownames -- not a data.frame with a GeneID column.
# Convert now; the alternative is a cryptic error two scripts from here.

stopifnot(colnames(counts_raw)[1] == "GeneID")

gene_ids   <- as.character(counts_raw$GeneID)
counts_mat <- as.matrix(counts_raw[, -1, drop = FALSE])
rownames(counts_mat) <- gene_ids

# `[, -1]` = "drop the first column", NOT "the last column" (PREP 1.3).
# drop = FALSE stops R collapsing a one-column result into a bare vector.

storage.mode(counts_mat) <- "integer"   # DESeq2 requires integers (PREP 4)

cat("Count matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")
cat("Samples:", paste(colnames(counts_mat), collapse = ", "), "\n\n")


# =============================================================================
# PART 2 -- LOAD THE GENE ANNOTATION
# =============================================================================

annot <- read.delim(
  path_annot,
  header       = TRUE,
  check.names  = FALSE,
  quote        = "",
  comment.char = ""
)

annot$GeneID <- as.character(annot$GeneID)   # must match the counts rownames type

# --- GOTCHA: blanks import as "" not NA (METHODS_LOG §3) ---------------------
# An empty string is a real value to R. It passes !is.na() happily. Normalize
# it now or every downstream missingness check quietly lies to you.
annot$EnsemblGeneID[annot$EnsemblGeneID == ""] <- NA

cat("Annotation:", nrow(annot), "genes\n")
cat("Missing Ensembl ID:", sum(is.na(annot$EnsemblGeneID)),
    sprintf("(%.1f%%)  -- expect 12,299 / 31.2%%\n\n",
            100 * mean(is.na(annot$EnsemblGeneID))))

# A NAMED VECTOR is base R's dict: index it with a character key.
#   symbol_of[["3064"]]  ->  "HTT"
symbol_of <- setNames(annot$Symbol, annot$GeneID)


# =============================================================================
# PART 3 -- SAMPLE METADATA
# =============================================================================

# Transcribed from METHODS_LOG §2. GSM8343537 is listed here ON PURPOSE so the
# "expected absence" check below is explicit rather than an accident.

meta_all <- data.frame(
  gsm     = c("GSM8343537", "GSM8343538", "GSM8343539",
              "GSM8343540", "GSM8343541", "GSM8343542", "GSM8343543",
              "GSM8343544", "GSM8343545", "GSM8343546", "GSM8343547"),
  label   = c("IC1_1", "IC1_2", "IC1_3",
              "HD_1", "HD_2", "HD_3", "HD_4",
              "KO_1", "KO_2", "KO_3", "KO_4"),
  group   = c("IC1", "IC1", "IC1",
              "HD", "HD", "HD", "HD",
              "KO", "KO", "KO", "KO"),
  passage = c(4, 5, 6,
              4, 5, 6, 7,
              4, 5, 6, 7)
)

missing_gsm <- setdiff(meta_all$gsm, colnames(counts_mat))   # PREP 8
cat("In GEO but absent from the NCBI counts:",
    if (length(missing_gsm)) paste(missing_gsm, collapse = ", ") else "none", "\n")

stopifnot(identical(missing_gsm, "GSM8343537"))
cat("  -> matches the documented absence of IC1_1. OK\n\n")

# stopifnot() is R's assert. Use it for anything that makes everything
# downstream meaningless if it turns out false.


# =============================================================================
# PART 4 -- SUBSET TO KO vs IC1 AND SET THE REFERENCE LEVEL
# =============================================================================

meta <- meta_all[meta_all$group %in% groups_keep &
                 meta_all$gsm %in% colnames(counts_mat), ]

# Note the trailing comma -- df[rows, ] (PREP 2). Without it you subset COLUMNS.

# --- The single most important line in this script (PREP 3) ------------------
meta$condition <- factor(meta$group, levels = groups_keep)
levels(meta$condition)                  # confirm: "IC1" "KO" -- IC1 is reference

# If this were left to R's alphabetical default it would be "IC1" "KO" anyway,
# by luck. Set it explicitly regardless -- the habit is what protects you.

meta$passage <- factor(meta$passage)    # for the ~ passage + condition extension

rownames(meta) <- meta$gsm              # rownames are load-bearing (PREP 2)

cat("Design after subsetting:\n")
print(table(meta$condition))
cat("Reference level:", levels(meta$condition)[1], "\n")
cat("Expect: IC1 = 2, KO = 4 (authors had 3 vs 4)\n\n")


# =============================================================================
# PART 5 -- ALIGN COUNTS TO METADATA
# =============================================================================

# DESeq2 matches counts columns to metadata rows BY POSITION, not by name.
# Wrong order = wrong analysis, no warning, no error. Reorder, then assert.

counts <- counts_mat[, meta$gsm, drop = FALSE]

# Indexing a matrix by a character vector subsets AND reorders in one step.

stopifnot(all(colnames(counts) == rownames(meta)))   # PREP 4, keep this forever
cat("Counts columns aligned to metadata rows. OK\n\n")


# =============================================================================
# PART 6 -- SANITY CHECKS
# =============================================================================

cat("=============================================================\n")
cat("SANITY CHECKS\n")
cat("=============================================================\n\n")

# --- 6.1 Counts are raw non-negative integers --------------------------------

stopifnot(all(counts == round(counts)))
stopifnot(!any(is.na(counts)))
stopifnot(min(counts) >= 0)
cat("[6.1] Non-negative integers, no NAs. OK\n\n")


# --- 6.2 Library sizes -------------------------------------------------------

lib_sizes <- colSums(counts)            # vectorized; never loop for this

cat("[6.2] Library size, millions of assigned reads:\n")
print(round(lib_sizes / 1e6, 1))
cat("  max/min ratio:", round(max(lib_sizes) / min(lib_sizes), 2),
    " -- above ~3x, look for a failed library before proceeding\n\n")


# --- 6.3 Detected genes ------------------------------------------------------

cat("[6.3] Genes with non-zero counts per sample:\n")
print(colSums(counts > 0))              # TRUE counts as 1 (PREP 1.5)
cat("  of", nrow(counts), "genes in the matrix\n\n")


# --- 6.4 Annotation coverage -------------------------------------------------

in_annot <- rownames(counts) %in% annot$GeneID
cat("[6.4] Count-matrix genes present in the GEO annotation:",
    sum(in_annot), sprintf("(%.1f%%)\n", 100 * mean(in_annot)))
cat("  -> write this into METHODS_LOG §9, Mapping\n\n")


# --- 6.5 CPM for eyeball checks ----------------------------------------------

# Raw counts are not comparable across samples with different library sizes.
# sweep() divides every column by its own library size -- pandas' df.div(s, axis=1).
# This is for LOOKING only. DESeq2 does its own normalization via size factors.

cpm <- sweep(counts, 2, lib_sizes, FUN = "/") * 1e6

is_ic1 <- meta$condition == "IC1"       # logical vector, reused below
is_ko  <- meta$condition == "KO"


# --- 6.6 PRIMARY CONTROL: HTT ------------------------------------------------
# The one check that stops the project if it fails.

cat("[6.6] PRIMARY CONTROL -- HTT (Entrez", HTT_ENTREZ, ")\n")
stopifnot(HTT_ENTREZ %in% rownames(counts))

htt_cpm <- cpm[HTT_ENTREZ, ]

print(data.frame(
  sample    = meta$label,
  condition = as.character(meta$condition),
  raw       = as.integer(counts[HTT_ENTREZ, ]),
  cpm       = round(htt_cpm, 1)
), row.names = FALSE)

mean_ic1 <- mean(htt_cpm[is_ic1])
mean_ko  <- mean(htt_cpm[is_ko])

cat("\n  mean CPM  IC1:", round(mean_ic1, 1), "  KO:", round(mean_ko, 1), "\n")
cat("  approx log2(KO/IC1):", round(log2((mean_ko + 1) / (mean_ic1 + 1)), 2), "\n")

if (mean_ko > mean_ic1) {
  cat("  -> HTT higher in KO, matching Table S3 (+", HTT_EXPECTED_L2FC, "). OK\n", sep = "")
  cat("     IC1's corrected allele is silenced (METHODS_LOG §5.9), so the\n")
  cat("     control is monoallelic. The KO frameshift still transcribes.\n")
  cat("     Protein is gone; message is not.\n\n")
} else {
  cat("  -> *** FAIL: HTT is NOT higher in KO ***\n")
  cat("     STOP. Check group assignments in PART 3 against the GEO sample\n")
  cat("     titles, then factor levels in PART 4.\n\n")
}


# --- 6.7 SECONDARY CONTROLS --------------------------------------------------

sec <- data.frame(
  symbol   = c(genes_up_expected, genes_down_expected),
  expected = c(rep("up",   length(genes_up_expected)),
               rep("down", length(genes_down_expected)))
)

# match() is the vectorized lookup -- pandas' Series.map(). Returns the POSITION
# of each symbol in annot$Symbol, or NA when absent.
sec$entrez <- annot$GeneID[match(sec$symbol, annot$Symbol)]

found <- !is.na(sec$entrez) & sec$entrez %in% rownames(counts)

sec$ic1 <- sec$ko <- sec$log2r <- NA
sec$ic1[found]   <- round(rowMeans(cpm[sec$entrez[found], is_ic1, drop = FALSE]), 1)
sec$ko[found]    <- round(rowMeans(cpm[sec$entrez[found], is_ko,  drop = FALSE]), 1)
sec$log2r[found] <- round(log2((sec$ko[found] + 1) / (sec$ic1[found] + 1)), 2)

sec$ok <- ifelse(is.na(sec$log2r), NA,
                 ifelse(sec$expected == "up", sec$log2r > 0, sec$log2r < 0))

cat("[6.7] SECONDARY CONTROLS (mean CPM -- eyeball, not a test)\n")
print(sec[, c("symbol", "entrez", "ic1", "ko", "log2r", "expected", "ok")],
      row.names = FALSE)

cat("\n  correct direction:", sum(sec$ok, na.rm = TRUE), "of", sum(!is.na(sec$ok)), "\n")
cat("  These are lowly expressed TFs, so one or two misses is noise. A\n")
cat("  MAJORITY going the wrong way means go back to PART 4 factor levels.\n\n")


# =============================================================================
# PART 7 -- SAVE FOR 02_deseq2.R
# =============================================================================

if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

saveRDS(
  list(
    counts     = counts,
    meta       = meta,
    annot      = annot,
    symbol_of  = symbol_of,
    lib_sizes  = lib_sizes,
    provenance = list(
      counts_file = path_counts,
      annot_file  = path_annot,
      contrast    = "KO vs IC1",
      loaded_at   = Sys.time(),
      session     = utils::sessionInfo()      # PREP 10, goes in the README
    )
  ),
  file = file.path(dir_out, "01_loaded.rds")
)

# saveRDS writes ONE object with full fidelity -- types, names, factor levels
# all survive. readRDS() in the next script gets it back exactly. Unlike
# save()/load(), it does not silently clobber variables in your workspace.

cat("Saved:", file.path(dir_out, "01_loaded.rds"), "\n")
cat("Next: scripts/02_deseq2.R\n")


# =============================================================================
# SELF-CHECK -- answer these before moving to 02
# =============================================================================
#
# 1. counts_raw[, -1] -- what does that do, and what would Python have done?
# 2. Why does this script convert the data.frame to a matrix at all?
# 3. PART 5 reorders columns and then asserts they match. Why is the assert
#    still worth having if you just did the reordering?
# 4. annot$EnsemblGeneID has 12,299 blanks. What breaks if you skip the
#    "" -> NA line and later run sum(is.na(annot$EnsemblGeneID))?
# 5. HTT comes back at 40 CPM in KO and 380 in IC1. Is that a pass?
#
# Answers: 1) drops the first column, not "everything but the last" (PREP 1.3).
#          2) DESeq2 needs a homogeneous integer matrix with gene rownames; a
#             data.frame with a GeneID column will not do.
#          3) It documents the invariant and catches the case where meta$gsm
#             later contains an ID absent from the counts -- the reorder itself
#             would silently produce an all-NA column.
#          4) You get 0. "" is a real value, so the blanks pass !is.na() and
#             your mapping-loss numbers in §9 come out wrong.
#          5) Yes -- log2 ratio approx -3.2. A KO transcript is degraded by NMD,
#             not eliminated, so reduced-but-present is the expected result.
# =============================================================================
