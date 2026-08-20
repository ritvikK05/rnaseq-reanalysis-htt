# =============================================================================
# PREP DAY: Python -> R, using the actual operations from Week 1
#
# Read top to bottom. Run each block. Total time: ~60-75 minutes.
#
# PARTS 1-5 need only base R. Run these first, they cost nothing.
# PARTS 6-10 need Bioconductor. Install runs in the background while you
#            work through 1-5.
#
# Every gotcha below is one I expect you to actually hit in Week 1.
# =============================================================================


# =============================================================================
# PART 0 -- INSTALL (start this FIRST, it takes 10-20 min, then keep reading)
# =============================================================================

# Bioconductor is a SEPARATE package ecosystem from CRAN. This is the single
# most common install confusion. install.packages() will NOT find DESeq2.

install.packages("BiocManager")          # CRAN: the Bioconductor installer
BiocManager::install(c(
  "DESeq2",            # differential expression
  "apeglm",            # the shrinkage estimator you'll use in lfcShrink()
  "EnhancedVolcano",   # one-call volcano plots -- do NOT hand-roll in ggplot2
  "clusterProfiler",   # GO / pathway enrichment
  "org.Hs.eg.db",      # human gene ID annotation database
  "tximport"           # imports Salmon output into R
))

# The `::` above means "call this function from that package without loading
# the whole package." Python's closest analogue is `from x import y`, used inline.

# Sanity check -- if this prints without error you are set up.
library(DESeq2)
packageVersion("DESeq2")


# =============================================================================
# PART 1 -- THE FIVE SYNTAX FACTS THAT WILL BITE YOU
# =============================================================================

# --- 1.1 Assignment uses <- not = --------------------------------------------
# `=` mostly works, but `<-` is the convention and you will read it everywhere.
x <- c(10, 20, 30, 40, 50)   # c() = "combine", the universal list constructor
x

# --- 1.2 R is 1-INDEXED and slicing is INCLUSIVE ------------------------------
x[1]        # 10   <- first element. Python would give you 20.
x[1:3]      # 10 20 30   <- THREE elements. Python's x[1:3] gives two.
# Rule of thumb: R's a:b is Python's range(a, b+1).

# --- 1.3 NEGATIVE INDEXING MEANS "DROP", NOT "FROM THE END" -------------------
# This one silently produces wrong answers instead of erroring. Watch:
x[-1]       # 20 30 40 50  <- dropped the first element!
            # In Python this would have returned 50.
# To get the last element:
x[length(x)]   # 50
tail(x, 1)     # 50, more readable

# --- 1.4 TRUE / FALSE / NULL / NA --------------------------------------------
TRUE; FALSE          # capitalized, no quotes. `True` is an undefined variable.
NA                   # missing value -- this is the one that matters (see 7.2)
NULL                 # absence of an object, closer to Python's None
is.na(NA)            # TRUE   <- you test for NA with is.na(), never with ==
NA == NA             # NA     <- NOT TRUE. Comparison with NA propagates NA.

# --- 1.5 Everything is vectorized by default ---------------------------------
# You almost never write loops in R. This is a feature, lean into it.
x * 2                # operates on all 5 elements, no map/comprehension needed
x > 25               # FALSE FALSE TRUE TRUE TRUE -- a logical vector
x[x > 25]            # 30 40 50 -- boolean masking, same idea as pandas
sum(x > 25)          # 3 -- TRUE counts as 1, handy for "how many genes pass?"
sum (x[x>25])

# =============================================================================
# PART 2 -- DATA FRAMES (your pandas intuition mostly transfers)
# =============================================================================

# Build a fake metadata table like the one you'll load in Week 1.
meta <- data.frame(
  sample    = c("S1", "S2", "S3", "S4"),
  condition = c("control", "control", "treated", "treated"),
  batch     = c("A", "B", "A", "B")
)

meta                 # print
str(meta)            # structure -- your df.info(). USE THIS CONSTANTLY.
head(meta, 2)        # df.head(2)
dim(meta)            # df.shape
nrow(meta); ncol(meta)
colnames(meta)       # df.columns

# --- Column access ---
meta$condition          # df['condition'] -- the $ is your workhorse
meta[["condition"]]     # identical, useful when the name is in a variable

# --- Row/column selection: df[rows, cols] -- note the COMMA ------------------
meta[1, ]               # first row, all columns   <- trailing comma required
meta[, "condition"]     # all rows, one column
meta[meta$condition == "treated", ]   # boolean row filter, like df[df.x == y]

# --- Row names are a real, load-bearing concept in R -------------------------
# DESeq2 uses them to match samples between your counts and metadata.
rownames(meta) <- meta$sample
meta


# =============================================================================
# PART 3 -- FACTORS: THE #1 SOURCE OF SILENTLY WRONG RESULTS
# =============================================================================
#
# A factor is a categorical variable with ordered "levels". DESeq2 uses the
# FIRST level as the reference (denominator) in every comparison.
#
# R sorts levels ALPHABETICALLY by default. So:
#
#   control / treated  -> "control" is first -> correct, by luck
#   KO / WT            -> "KO" is first      -> WT is now the numerator,
#                                               every log2FC sign is FLIPPED
#
# Nothing errors. Your volcano plot looks fine. Your biology is backwards.
# Papers have been corrected over this.

# Watch it happen:
bad <- factor(c("WT", "KO", "WT", "KO"))
levels(bad)      # "KO" "WT"  <- KO is the reference. Probably not what you want.

good <- c("WT", "KO", "WT", "KO")
good_f <- factor(good, levels = c("WT", "KO"))
levels(good_f)

# ALWAYS set levels explicitly. Reference goes FIRST.
meta$condition <- factor(meta$condition, levels = c("control", "treated"))
levels(meta$condition)     # "control" "treated" -- treated vs control. Correct.

meta$batch <- factor(meta$batch)

# Make this a reflex: every time you create a design variable, set its levels
# on the next line, then print levels() to confirm.


# =============================================================================
# PART 4 -- MAKE A TOY DATASET SO YOU CAN RUN THE REAL WORKFLOW TODAY
# =============================================================================
# No downloads, no SRA, no waiting. 500 genes, 4 samples, 50 genes truly DE.

set.seed(42)                       # np.random.seed()
n_genes <- 500

counts <- matrix(
  rnbinom(n_genes * 4, mu = 100, size = 1/0.2),   # negative binomial, like real
  nrow = n_genes,
  ncol = 4
)
# Spike a real signal into the first 50 genes for the treated samples:
counts[1:50, 3:4] <- counts[1:50, 3:4] * 4

rownames(counts) <- paste0("GENE", 1:n_genes)   # paste0 = f-string / string concat
colnames(counts) <- meta$sample

counts[1:5, ]      # peek at the top-left corner


# --- GOTCHA: DESeq2 requires INTEGER counts ----------------------------------
# Salmon gives you estimated (fractional) counts. tximport handles the
# conversion for you -- but if you ever build a matrix by hand:
counts <- round(counts)
storage.mode(counts) <- "integer"

# --- GOTCHA: sample order must match. DESeq2 does not reorder for you. -------
# Make this check a permanent part of your script. It costs nothing and it
# catches a class of bug that is otherwise nearly invisible.
stopifnot(all(colnames(counts) == rownames(meta)))


# =============================================================================
# PART 5 -- S4 OBJECTS: WHY dds$counts DOESN'T WORK
# =============================================================================
#
# Bioconductor objects are S4. You do NOT reach into them with a dot or $.
# You call ACCESSOR FUNCTIONS that take the object as an argument.
#
#   Python:  dds.counts        obj.attribute
#   R (S4):  counts(dds)       accessor(obj)
#
# Once you internalize "it's functions all the way down," S4 stops being weird.
# Common accessors you'll use: counts(), colData(), assay(), rowData(),
# sizeFactors(), design(), resultsNames()


# =============================================================================
# PART 6 -- THE ACTUAL DESeq2 WORKFLOW (this is the whole thing)
# =============================================================================

library(DESeq2)

# Build the object. The design formula is R's model syntax: ~ means "modeled by".
# ~ condition          -> test condition
# ~ batch + condition  -> control for batch, test condition (LAST term is tested)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ condition
)

# Pre-filter low-count genes. Not required, but speeds things up and reduces
# multiple-testing burden. keep is a logical vector -> boolean mask.
keep <- rowSums(counts(dds)) >= 10       # note: counts(dds), not dds$counts
dds  <- dds[keep, ]
nrow(dds)                                 # how many genes survived

# Run it. This one call does size-factor normalization, dispersion estimation,
# and the Wald test. It prints its progress -- that's normal, not an error.
dds <- DESeq(dds)

# What comparisons are available?
resultsNames(dds)     # "Intercept"  "condition_treated_vs_control"
# ^ Read that name out loud. "treated_vs_control" confirms your reference level
#   is correct. If it says "control_vs_treated", go back to PART 3.


# =============================================================================
# PART 7 -- RESULTS, AND THE padj NA TRAP
# =============================================================================

# lfcShrink gives you moderated fold changes -- more reliable ranking for
# low-count genes. Use this, not plain results(), for anything you plot or rank.
res <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")

head(res)
summary(res)          # quick tally of up/down at the default alpha

# --- 7.1 It's not a data.frame. Convert when you want to manipulate it. ------
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)         # promote rownames to a real column

# --- 7.2 THE TRAP: padj contains NA ------------------------------------------
# DESeq2 sets padj to NA for genes filtered out by independent filtering or
# flagged as outliers. NA is NOT FALSE. This bites you two ways:

sum(is.na(res_df$padj))                 # see how many

# WRONG -- returns NA rows full of NA, silently inflating your gene count:
# sig <- res_df[res_df$padj < 0.05, ]

# RIGHT -- which() drops NAs:
sig <- res_df[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1), ]
nrow(sig)

# Equally right, more explicit:
sig <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

# --- 7.3 Sorting: order() returns INDICES, not sorted values -----------------
sig <- sig[order(sig$padj), ]           # ascending; use decreasing = TRUE to flip
head(sig[, c("gene", "log2FoldChange", "padj")], 10)


# =============================================================================
# PART 8 -- CONCORDANCE ANALYSIS (your Day 3 deliverable, in 6 lines)
# =============================================================================

my_genes    <- sig$gene
paper_genes <- paste0("GENE", c(1:40, 501:520))   # stand-in for the published list

shared    <- intersect(my_genes, paper_genes)     # recovered
missed    <- setdiff(paper_genes, my_genes)       # theirs, not mine
extra     <- setdiff(my_genes, paper_genes)       # mine, not theirs

length(shared); length(missed); length(extra)
round(100 * length(shared) / length(paper_genes), 1)   # % recovery

# The %in% operator is your `in` -- useful for tagging:
res_df$in_paper <- res_df$gene %in% paper_genes

# GOTCHA you WILL hit with real data: Ensembl IDs carry version suffixes
# (ENSG00000141510.16). The paper's list may not have them. Strip before matching:
#   res_df$gene <- sub("\\..*$", "", res_df$gene)
# A concordance of 0% almost always means an ID-format mismatch, not biology.


# =============================================================================
# PART 9 -- PLOTS (use the one-call versions, skip ggplot2 for now)
# =============================================================================

library(EnhancedVolcano)

EnhancedVolcano(
  res_df,
  lab             = res_df$gene,
  x               = "log2FoldChange",
  y               = "padj",
  pCutoff         = 0.05,
  FCcutoff        = 1,
  title           = "Treated vs Control"
)

# MA plot comes free with DESeq2:
plotMA(res, ylim = c(-5, 5))

# Save a figure -- note you must explicitly close the device with dev.off():
png("volcano.png", width = 1800, height = 1500, res = 200)
plotMA(res, ylim = c(-5, 5))
dev.off()

# Learn ggplot2 later. It is genuinely worth knowing, but it is a whole mental
# model (grammar of graphics) and it is NOT what Week 1 is testing.


# =============================================================================
# PART 10 -- REPRODUCIBILITY (goes in your README)
# =============================================================================

sessionInfo()        # full package versions -- paste into your README

# renv is R's virtualenv/poetry. Run these in your project directory:
#   renv::init()        # start tracking
#   renv::snapshot()    # write renv.lock  <- commit this file
#   renv::restore()     # what someone else runs to reproduce you


# =============================================================================
# SELF-CHECK -- if you can answer these without scrolling up, you're ready
# =============================================================================
#
# 1. x <- c(5,6,7,8). What is x[-2]? What is x[2:3]?
# 2. Your conditions are "KO" and "WT". Write the line that makes WT the
#    reference, and the line that verifies it worked.
# 3. Why does res_df[res_df$padj < 0.05, ] give you rows of NA?
# 4. How do you get the count matrix out of a DESeqDataSet called dds?
# 5. Your concordance with the paper is 0%. What do you check first?
#
# Answers: 1) 5 7 8 (drops 2nd); 6 7 (inclusive, 1-indexed).
#          2) meta$cond <- factor(meta$cond, levels=c("WT","KO")); levels(meta$cond)
#          3) padj has NAs; NA < 0.05 is NA, which indexes a row of NAs. Use which().
#          4) counts(dds) -- S4 accessor function, not dds$counts.
#          5) Gene ID format: Ensembl version suffixes, or symbol vs Ensembl ID.
# =============================================================================
