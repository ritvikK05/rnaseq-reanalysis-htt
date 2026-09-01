# RNA-seq Reanalysis — Methods Log

Running record of decisions, parameters, and known discrepancy sources.
Append as you go; this becomes the README on Days 4–5.

---

## 1. Dataset

| Field | Value |
|---|---|
| GEO accession | GSE270472 (total RNA-seq) |
| Companion | GSE270473 (small RNA-seq — not used) |
| Paper | Kozłowska et al., *Cell & Bioscience* 15:100 (2025) |
| DOI | 10.1186/s13578-025-01443-5 |
| PMID | 40635054 |
| License | Open access (CC BY) |
| Organism | *Homo sapiens* |
| System | iPSC-derived neural stem cells (NSCs), isogenic lines |
| Platform | Illumina NovaSeq 6000, paired-end 2×100, ~50M reads/sample |
| Library prep | KAPA RNA HyperPrep + RiboErase (HMR) — total RNA, rRNA-depleted, **not** poly-A selected |

**Comparison chosen:** KO-NSC vs IC1-NSC (HTT knockout vs isogenic control)

Rationale: larger effect size than the HD comparison (1,401 vs 822 DEGs as
reported; 1,464 vs 835 as recomputed — see §4), and it is the paper's cleanest
loss-of-function contrast.

Note on library prep: RiboErase total-RNA libraries capture non-coding
transcripts that poly-A selection would miss. This explains why their results
table contains a large ncRNA fraction, and is relevant to gene-universe
comparisons.

---

## 2. Sample table

| GSM | Label | Group | Passage | In NCBI counts? |
|---|---|---|---|---|
| GSM8343537 | IC1_1 | IC1 | 4 | **NO — absent** |
| GSM8343538 | IC1_2 | IC1 | 5 | yes |
| GSM8343539 | IC1_3 | IC1 | 6 | yes |
| GSM8343540 | HD_1 | HD | 4 | yes |
| GSM8343541 | HD_2 | HD | 5 | yes |
| GSM8343542 | HD_3 | HD | 6 | yes |
| GSM8343543 | HD_4 | HD | 7 | yes |
| GSM8343544 | KO_1 | KO | 4 | yes |
| GSM8343545 | KO_2 | KO | 5 | yes |
| GSM8343546 | KO_3 | KO | 6 | yes |
| GSM8343547 | KO_4 | KO | 7 | yes |

**Effective design: 2 control vs 4 KO** (authors used 3 vs 4).

---

## 3. Input files

### Counts
NCBI-generated raw count matrix from GEO ("Download RNA-seq counts").

- Identifiers: **Entrez Gene IDs**
- Annotation basis: **NCBI RefSeq, GRCh38.p13**
- Verified integer: `all(counts == round(counts))` → TRUE
- 10 sample columns (GSM8343537 missing)

### Annotation
`Human.GRCh38.p13.annot.tsv.gz` (GEO, GRCh38.p13)

- 39,376 rows
- Columns: `GeneID` (Entrez), `Symbol`, `Description`, `Synonyms`, `GeneType`,
  `EnsemblGeneID`, `Status`, coordinates, GO term columns
- Version-matched to the counts file (both GRCh38.p13)
- `EnsemblGeneID` empty for 12,299 genes (31.2%) — composition in §9
- Note: blanks import as `""`, not `NA`. Normalize on load:
  `annot$EnsemblGeneID[annot$EnsemblGeneID == ""] <- NA`
- Note: NCBI names this file with periods, not underscores. Do not rename it —
  the original filename is provenance.

### Truth table
Supplementary Table 3 (`13578_2025_1443_MOESM5_ESM.xlsx`)

- Sheet `KO-NSC vs IC1-NSC` → **`skip = 1`** — 31,264 rows
- Sheet `HD-NSC vs IC1-NSC` → **`skip = 2`** — ~28,057 rows
  (extra title row — the offset differs between sheets!)
- Columns: `gene_ID`, `Gene_name`, `log2FoldChange`, `Fold_change`, `P-value`,
  `P_adjusted`, `Description`
- Identifiers: **Ensembl (v102), no version suffix**
- Full results, not just significant genes

### Not used
`GSE270472_HD_KO_NSC.xlsx` contains TPM values only. TPM is normalized and
fractional — DESeq2 requires raw integer counts. Retained only for eyeballing
individual genes.

### Provenance
- Counts file: `GSE270472_raw_counts_GRCh38.p13_NCBI.tsv.gz`
- Annotation file: `Human.GRCh38.p13.annot.tsv.gz`
- Both downloaded 2026-08-20 from GEO GSE270472 → "Download RNA-seq counts"
- Supplementary Table 3 (MOESM5) and Supplementary Methods (MOESM13)
  downloaded 2026-08-20 from
  https://link.springer.com/article/10.1186/s13578-025-01443-5

---

## 4. Original authors' parameters

### Pipeline (confirmed, Supplementary Methods / MOESM13)

| Step | Tool |
|---|---|
| QC | FastQC v0.11.9 |
| Trim / filter | BBDuk 2 v37.02 |
| rRNA removal | bowtie2 v2.3.5.1 (`-X 1000`, `--un-conc`); rRNA-aligned reads discarded |
| Quantification | **RSEM v1.3.1** (bowtie2 backend, default settings) |
| Reference | **Ensembl v102**, GRCh38 |
| Differential expression | **DESeq2 v1.30.0** |
| GO (main text) | PANTHER v17, Fisher's exact + FDR, p < 0.05 |
| GO (supplementary time-course) | clusterProfiler v4.10.0 `enrichGO`, BH correction |

Bioinformatics performed by IDEAS4BIOLOGY (external vendor).

### Thresholds and reported results
- **Thresholds:** |log2FC| > 1.5 AND padj < 0.05
  (note: 1.5, *not* 1 — stricter than the common default)
- **Reported DEGs:** 1,401 (KO vs IC1), 822 (HD vs IC1); 331 overlap
  - KO: 1,231 up / 170 down
  - HD: 387 up / 435 down

### Validation of published thresholds
Applying the stated cutoffs (|log2FC| > 1.5, padj < 0.05) to Supplementary
Table 3 yields **1,464** DEGs for KO vs IC1 (1,290 up / 174 down) and **835**
for HD vs IC1, versus 1,401 and 822 reported in the text.

Ruled out as causes: duplicate gene IDs (n=0) and unannotated genes (n=0).
A stricter cutoff of 1.55 gives 1,422, and only 42 genes fall between 1.5 and
1.55 — so no single unstated threshold explains the gap. Both comparisons
overshoot by a similar proportion (4.5% and 1.6%), suggesting a systematic
rather than incidental cause.

**Decision:** the 1,464-gene list derived directly from the published table is
used as the reference set for this reanalysis.

### DE framework — inference and correction
Supplementary Table 3 reports adjusted p-values for all 31,264 genes with
**zero NAs**, which is inconsistent with DESeq2 defaults. My initial inference
was that a different package (edgeR or limma-voom) had been used.

**The Supplementary Methods refute this: they used DESeq2 v1.30.0.** The
absence of NAs therefore indicates the analysis was run with
`independentFiltering = FALSE` and `cooksCutoff = FALSE`, not a different
package. Still a real methodological difference from the defaults used here,
but a parameter choice rather than a framework change.

---

## 5. Known discrepancy sources

The core of the writeup. Each is a candidate explanation for imperfect concordance.

1. **Missing control replicate.** NCBI's pipeline excluded GSM8343537, leaving
   n=2 controls vs the authors' n=3. Reduces power; expect to miss borderline
   genes.

2. **Different annotation source.** Theirs: Ensembl v102. Mine: NCBI RefSeq
   GRCh38.p13. These are different annotation *sources*, not merely different
   releases — different gene models, different gene boundaries, different gene
   universes. Likely a major contributor, and the root cause of the Entrez ↔
   Ensembl mapping problem in #4.

3. **Different quantifier.** RSEM v1.3.1 (bowtie2 backend, EM-based
   transcript-level estimation, probabilistic assignment of multi-mapping
   reads) vs NCBI's count-based pipeline, which typically discards
   multi-mappers. Expect systematic differences for paralogues and repeat-rich
   loci.

4. **ID mapping loss.** Counts use Entrez, truth table uses Ensembl.
   `AnnotationDbi::select()` returns 1:many mappings; record the final mapping
   rate and how duplicates were resolved. Loss is concentrated in non-coding
   and pseudogene classes — only 178 protein-coding genes lack an Ensembl ID,
   so impact on a protein-coding-dominated DEG list should be minimal.

5. **rRNA depletion step.** They removed rRNA-aligned reads with bowtie2 before
   quantification. NCBI's pipeline does not include this step, so library
   composition and therefore size factors differ.

6. **Independent filtering disabled.** Same package (DESeq2), but they ran with
   `independentFiltering = FALSE` and `cooksCutoff = FALSE`. Genes assigned NA
   in this reanalysis were still tested in theirs — **unrecoverable by
   construction**. Quantify in §9. Note this cuts both ways. Independent 
   filtering removes low-count genes from the multiple-testing burden before BH 
   correction, so this reanalysis tests fewer hypotheses (29,590 vs their 31,264) 
   and each surviving gene carries a lighter correction. That is a power *gain* 
   offsetting the loss in §5.1 — and is the most likely reason this reanalysis 
   calls 2,252 DEGs against their 1,464 despite having one fewer control replicate.

7. **Shrinkage.** Using `lfcShrink(type="apeglm")`; their table contains
   unshrunken estimates. Evidence: top KO hit FKBPL has log2FC ≈ 22
   (fold change ~4.1 million) — a near-zero-denominator artifact that apeglm
   will pull hard toward zero. **Expect substantially different top-ranked
   genes even where significance calls agree.**

8. **DESeq2 version drift.** Theirs v1.30.0 (2020); mine current (Bioconductor
   3.23). Defaults and shrinkage behavior have changed across releases.

9. **Control line caveat (authors' own).** IC1 was originally intended as the
   isogenic control, but the corrected allele turned out to be silenced,
   producing monoallelic *HTT* expression. The authors switched to IC2 for
   later experiments — but the RNA-seq comparison still uses IC1. So the
   "control" is not a clean wild-type.

   This is also what produces the *positive* HTT log2FC (+1.34) in Table S3 —
   see §7. A monoallelic control against a frameshift-edited KO that still
   transcribes gives more HTT message in the knockout, not less.

10. **Passage as unmodeled covariate.** Passages are matched across lines
    (IC1: 4–6; KO: 4–7), so passage functions like a batch variable.
    Baseline model is `~ condition`; test `~ passage + condition` as an
    extension.

    PCA (§9) does not support an additive passage term. KO_1/2/3 (p4–p6)
    cluster within ~1.5 units on PC2 while KO_4 (p7) sits ~28 units away,
    and the two IC1 samples trend in the opposite direction. An additive
    term assumes one shared passage effect across both lines, which the data
    contradict; the interaction is unfittable at n=6. Run as an exploratory
    sensitivity check using *continuous* passage (1 df, leaving 3 residual)
    rather than a factor (3 df, leaving 1), and report it as exploratory
    rather than as a corrected model.

11. **Enrichment tool.** clusterProfiler v4.20.0 vs their PANTHER v17 (main
    text). Different gene universes and algorithms — GO results will differ for
    tool reasons alone. Note they used clusterProfiler v4.10.0 for the
    supplementary time-course analysis, so a closer comparison is possible
    there.

    **Outcome (04).** The recovered DEG set reproduces the paper's biology
    convincingly — skeletal system development at p.adj 2.5e-27, plus ECM
    organization, embryonic morphogenesis, ossification and pattern
    specification. Agreement therefore holds at the level of biology and not
    only of gene identifiers, which is the check that matters. Term-by-term
    comparison against their PANTHER output is not attempted: different
    universe, different algorithm, and their published term list is
    summarised rather than exhaustive.

12. **Reference DEG list is dominated by near-zero-baseline genes.**
    Prediction registered before running 02, from two observations: seven of
    the eight RT-qPCR-validated controls sit at 0.0–0.4 CPM in IC1 (§7), and
    the top reported hit FKBPL has log2FC ≈ 22 (§5.7). Their 1,290-up /
    174-down asymmetry is therefore likely driven by genes switching on from
    zero rather than by symmetric regulation. Consequence: apeglm shrinks
    exactly this class hardest, so recovery should be *lowest* in the
    category that makes up most of their DEG list. Test by stratifying
    recovery rate by baseMean in `03_concordance.R`.

    **Outcome: not supported.** Recovery by baseMean quartile came back flat
    (81.3 / 80.3 / 80.3 / 85.2, Q1 to Q4), so low-expression genes are
    recovered about as reliably as high-expression ones.

    Two things this doesn't rule out. The stratification only includes genes
    tested in both analyses, so the 244 unmappable reference DEGs are left
    out by construction, and those are where a low-expression effect would
    most likely show up. So I didn't really test the prediction on the group
    it was about. Also, the related test in §9 (gene-type composition of the
    extra DEGs) came back null too, and the up/down asymmetry that made me
    write this prediction mostly went away after ID mapping. See §5.13.

13. **What explains the DEG excess (revised after 05).** I call 2,252 DEGs
    against their 1,464. I registered four explanations and all four
    failed:

    - baseline expression (§5.12 — recovery flat across baseMean quartiles)
    - gene-class composition (§9 — 34.2% vs 34.3% ncRNA+pseudogene share)
    - independent filtering (05A — disabling it moved the count from 2,252
      to 2,226, closing 3.3% of the 788-gene gap)
    - the extras being real signal I had the power to find (04 — see below)

    So the excess has no explanation downstream of the count matrix. The
    sensitivity runs in 05 show this directly: recovery is 68.2 / 68.0 /
    62.6 / 67.3 across the four variants, so nothing I can change at the
    `results()` stage improves agreement with Table S3, and the variant that
    moves the DEG count most actually makes recovery worse.

    That leaves the upstream differences: §5.2 (RefSeq vs Ensembl v102 gene
    models), §5.3 (NCBI's count-based pipeline vs RSEM's EM assignment of
    multi-mapping reads), and §5.5 (rRNA depletion changing library
    composition and so the size factors). These give different counts from
    the same reads, and no `results()` parameter can undo that. To confirm
    it I would have to re-quantify from FASTQ with RSEM against Ensembl
    v102, which was outside what I could do here. That's the next experiment
    if I come back to this.

    **A second possibility I can't test.** 05B shows the DEG count moves 25%
    when I drop one of six samples (§9, Extensions). A design that sensitive
    to which samples are included could produce a gap this size from sample
    membership alone, and my dataset differs from the authors' by exactly
    one sample, the missing IC1_1 (§5.1). That wouldn't require any pipeline
    difference at all. I can't separate it from the quantification
    explanation without the excluded replicate.

    **04 supports the upstream explanation.** GO enrichment on the 786 extra
    DEGs gives 3 marginal terms (best p.adj 0.0086) against 256 for the
    recovered set, and the two sets share zero GO terms (Jaccard 0.000). So
    the extras aren't coherent biology, and they aren't the same biology as
    the shared set either, which makes it unlikely they're real signal I was
    sensitive enough to catch. The one term that does come up,
    miRNA-mediated post-transcriptional gene silencing, is non-coding RNA
    machinery — which is where §5.2 / §5.3 / §5.5 would predict the
    differences to land.

    **One thing did get confirmed.** Disabling filtering and Cook's took my
    padj NAs from 2,326 to exactly 0, matching Table S3. I had inferred in
    §5.6 that the authors disabled these settings purely from the fact that
    their table has no NAs, and this confirms it.

---

## 6. Analysis decisions

- **Reference level:** `factor(condition, levels = c("IC1", "KO"))`
  → verify with `resultsNames(dds)` = `condition_KO_vs_IC1`
- **Pre-filter:** `rowSums(counts(dds)) >= 10`
- **Shrinkage:** `lfcShrink(coef = "condition_KO_vs_IC1", type = "apeglm")`
- **NA handling:** use `which()` — `padj` contains NAs from independent filtering
- **Thresholds:** match the paper (|log2FC| > 1.5, padj < 0.05)
- **Namespacing:** `AnnotationDbi::select()`, `dplyr::filter()` — several
  packages mask these
- **Filtering:** keeping DESeq2 defaults (independent filtering ON) rather than
  matching their disabled setting. Rationale: defaults are the better practice,
  and the difference becomes a measurable discrepancy source (§5.6) rather than
  a hidden one. Optionally rerun with filtering off as a sensitivity check.
- **Which log2FC is thresholded:** the |log2FC| > 1.5 cutoff is applied to
  UNSHRUNKEN estimates for the matched-threshold DEG set, because Table S3
  contains unshrunken values (§5.7). Thresholding apeglm-shrunken values
  against their unshrunken ones would measure shrinkage rather than
  concordance. Shrunken estimates are used for ranking, plots, and effect-size
  reporting, and the shrunken DEG count is reported as a sensitivity analysis.
  The gap between the two counts quantifies §5.7.
- **Independent filtering target:** `results(alpha = 0.05)` to match the padj
  cutoff. The 0.1 default optimises the filter for a threshold not being used.

---

## 7. Sanity checks

### Primary positive control
**HTT** (Entrez 3064, ENSG00000197386) should show strong **UPregulation**
in KO vs IC1. Table S3 reports **log2FC = +1.34**.

This is counterintuitive for a knockout and was initially logged backwards.
Two things explain it:

- The KO is a frameshift edit, not a deletion. The locus still transcribes;
  the message is nonfunctional. Protein is absent, transcript is not.
- IC1's corrected allele is silenced (§5.9), so the *control* expresses HTT
  monoallelically. The comparison is effectively two alleles vs one.

Confirmed at the raw-count stage (01_load_data.R §6.6): mean CPM 57.5 in IC1
vs 163.6 in KO, approximate log2 ratio **+1.49** against their +1.34. Close
agreement on unnormalised CPM with n=2 vs n=4.

If HTT is not clearly higher in KO, stop — sample labelling or the count
matrix is wrong.

### Secondary controls
Expected **upregulated** in KO (validated by the authors via RT-qPCR):
`TWIST1`, `SIX1`, `TBX1`, `TBX15`, `MSX2`, `MEOX2`, `FOXD1`

Expected **downregulated**: `PAX6` (paper reports log2FC = −2.38, KO vs IC1)

If these don't appear with the right direction, check factor levels first.

All 8 secondary controls confirmed at the CPM stage, correct direction:
TWIST1 +2.58, SIX1 +3.64, TBX1 +0.72, TBX15 +4.24, MSX2 +3.73, MEOX2 +0.45,
FOXD1 +3.66, PAX6 −1.69 (paper reports −2.38).

All seven upregulated genes rise from near-zero baseline (0.0–0.4 CPM in IC1),
so log2 ratios are inflated by the +1 pseudocount and should not be compared
to DESeq2 output directly. See §5.12 for the consequence.

---

## 8. Environment

- macOS (Apple Silicon), R 4.6.1, Bioconductor 3.23
- DESeq2, apeglm, EnhancedVolcano, clusterProfiler 4.20.0, org.Hs.eg.db, tximport
- Capture `sessionInfo()` output for the README
- Environment captured in `notes/sessionInfo.txt` (R 4.6.1, Bioconductor 3.23)
- No `renv` lockfile; package versions recorded via sessionInfo instead

**Fallback rule:** if still fighting R syntax at end of Day 2, switch to
pydeseq2. No exceptions.

---

## 9. Results

### Completed
- [x] Reference DEG count recomputed from Table S3 (KO): **1,464**
      (paper reports 1,401 — see §4)
- [x] Reference DEG count recomputed from Table S3 (HD): **835**
      (paper reports 822)
- [x] Annotation coverage: 39,376 genes, 12,299 (31.2%) lacking an Ensembl ID
      — 9,206 ncRNA, 1,610 blank type, 666 pseudogene, 437 tRNA, 130 snoRNA,
      28 snRNA, 20 other, 19 rRNA, 5 unknown, and only **178 protein-coding**
- [x] Original pipeline identified: RSEM v1.3.1 / Ensembl v102 / DESeq2 v1.30.0
      (see §4)
- [x] Library sizes: 45.1–58.2M assigned reads, max/min ratio 1.29 —
      no sample excluded
- [x] Detected genes (non-zero): 29,213–31,602 of 39,376 (~75–80%).
      Detection saturates; depth and detection track together across groups
      (GSM8343539 at 47.6M/30,874 vs GSM8343545 at 47.5M/30,911), so the
      higher KO depth is sequencing, not biology
- [x] Secondary controls confirmed at CPM stage: 8/8 correct direction (§7)
- [x] Dispersion fit inspected: standard shape — trend descends then flattens
      to the biological floor; per-gene estimates scatter around the fitted
      trend with final values shrunk hard onto it, as expected at n=2 controls
- [x] PCA (VST, blind): PC1 85%, clean KO/IC1 separation — the knockout
      dominates total variance. PC2 9%, and does NOT show a passage gradient:
      KO_1/2/3 (p4–p6) cluster within ~1.5 units while KO_4 (p7) sits ~28
      units away, so PC2 reflects one divergent sample rather than a trend.
      The two IC1 replicates span nearly the full PC2 range (+14.7 to −14.2),
      making them the most dissimilar pair in the dataset — the sole estimate
      of within-control variability comes from an unusually spread pair,
      a concrete mechanism for the conservatism predicted in §5.1.

### Mapping
- [x] Coverage of count-matrix genes via GEO annotation: **100%** (39,376/39,376)
- [x] Additional genes recovered via `org.Hs.eg.db`: 984. Full sequence:
      22,013 of 29,590 tested genes carried an Ensembl ID from the GEO
      annotation (74.4%); org.Hs.eg.db recovered 984 more, giving 22,997
      (77.7%); collapsing 118 duplicate Ensembl IDs leaves **22,879 unique
      genes** entering the comparison. The remaining 6,593 cannot
      participate regardless of their biology — §5.2 from this side, with
      the 244 unmappable reference DEGs being the same problem from theirs.
- [x] Duplicate resolution method: 118 Ensembl IDs received more than one
      Entrez mapping; kept the entry with the lowest padj, tie-broken on
      highest baseMean. 1:many mappings from `AnnotationDbi::select()` were
      collapsed by keeping the first Ensembl ID per Entrez ID — arbitrary,
      but these are genes the primary annotation could not map at all,
      overwhelmingly ncRNA and pseudogenes (§9), so the choice has little
      leverage on a protein-coding-dominated DEG list.
- [x] Reference DEGs unmappable by construction: 244 of 1,464 (16.7%) —
      never entered the tested gene universe here, so unrecoverable
      regardless of statistical choices. This is the gap between the 68.2%
      and 81.8% recovery rates.

### My analysis
- [x] Genes tested after pre-filtering: 29590 (vs 31,264 in their table)
- [x] My DEG count at matched thresholds: 2252 (1483 up / 769 down)
- [x] DEG count on shrunken estimates: 1901 (1331 up / 570 down)
      — difference from the unshrunken count is the magnitude of §5.7
- [x] Genes with padj = NA: 2326 (7.9% of tested)
      — the §5.6 ceiling on recovery; their run had zero
- [x] DEG count on shrunken (apeglm) estimates: **1,901** (1,331 up /
      570 down) vs 2,252 unshrunken. The 351-gene gap is the magnitude of
      §5.7 — shrinkage alone removes 15.6% of the calls at matched
      thresholds. Reported as a sensitivity figure; the unshrunken set is
      used for concordance because Table S3 is unshrunken (§6).
- [x] HTT log2FC and padj: +1.44 (expect ~+1.34, Table S3; +1.49 at CPM stage)
- [x] PC1 variance explained: 85%  (VST, blind)
- [x] Up/down asymmetry differs sharply: 1,483 up / 769 down (1.9:1) vs their
      1,290 / 174 (7.4:1). Up counts agree within 15%; down count is 4.4x
      theirs. Direction-specific, so not explained by a uniform power
      difference. Candidates: §5.5 (rRNA depletion → size factors),
      §5.2/§5.3 (annotation and quantifier). Test in 03 by stratifying the
      excess downregulated genes by gene_type.
      **Resolved in 03:** after ID mapping and duplicate collapse the extra
      DEGs split 424 down / 362 up with near-identical gene-type composition
      (34.2% vs 34.3% ncRNA+pseudogene). The asymmetry was largely an
      artifact of the unmapped gene set, not a direction-specific mechanism.
      See §5.13.

### Concordance
- [x] Recovered (shared): 998
- [x] Missed (theirs, not mine): 466
- [x] Extra (mine, not theirs): 786
- [x] Recovery rate: 68.2% of all 1,464 reference DEGs; 81.8% of the 1,220
      that were mappable and tested here. The 13.6-point gap is §5.4/§5.2 —
      genes never testable in this reanalysis, not statistical disagreement.
- [x] My DEG count after mapping and duplicate collapse: 1,784 (from 2,252
      unmapped) — 468 lost to unmappable IDs or duplicate resolution
- [x] Of shared DEGs, fraction with matching log2FC sign: 100.0% (998/998)
      — rules out reference-level or labelling error entirely
- [x] Spearman correlation of log2FC across all shared genes: 0.617
      (19,190 genes tested in both); 0.896 restricted to the 998 recovered
      DEGs. The gap is expected: non-DE genes have little signal to agree
      about, so two independent pipelines rank them differently. High
      agreement where signal exists, moderate where it does not.
- [x] Of missed genes, fraction assigned NA by DESeq2 filtering: 10 of 222
      tested (4.5%). Much smaller than §5.6 anticipated — 2,326 genes have
      padj = NA overall, but almost none are genes the authors called
      significant. The filtering ceiling is real but nearly irrelevant here.
- [x] Of missed genes, fraction with padj between 0.05 and 0.15: 72 of 222
      (32.4%) — threshold sensitivity from n=2 controls (§5.1)
- [x] Missed-gene decomposition: 466 total = 244 never tested here (52%,
      structural — §5.2/§5.4) + 222 tested but not called. Of the 222:
      10 NA, 72 near-miss, 32 significant but below the fold-change cutoff,
      leaving ~108 genuine disagreements (~7% of the 1,464 reference set).
- [x] Recovery rate stratified by baseMean quartile: 81.3 / 80.3 / 80.3 / 85.2
      (Q1 lowest to Q4 highest). §5.12 predicted recovery would rise with
      expression level; it does not. Recovery is essentially flat, with Q1
      exceeding Q2 and Q3 and under five points separating the extremes.
      The prediction is not supported.
- [x] Gene-type composition of "extra" DEGs, up vs down: essentially identical
      (ncRNA+pseudogene share 34.2% down vs 34.3% up). The 786 extras split
      424 down / 362 up — close to balanced, so the raw up/down asymmetry
      noted on Day 2 largely dissolves after mapping and duplicate collapse
      (2,252 → 1,784 DEGs). No compositional signature; §5.5 and §5.2 are
      not implicated by this test. Only snoRNA (14/0) and snRNA (9/0) are
      one-sided, and at those counts that is suggestive at best.
      
### Enrichment
- [x] Set A (1,483 up): skeletal system development (p.adj 2.4e-23),
      extracellular matrix organization, vasculature and heart development,
      embryonic morphogenesis — the neural crest / mesoderm program the paper
      reports, recovered independently.
- [x] Set B (769 down): only 2 terms, but both on target — central nervous
      system neuron differentiation and regulation of Notch signaling
      (both p.adj 0.018). PAX6's reduction predicts loss of NSC identity, and
      Notch is the canonical pathway maintaining progenitors in an
      undifferentiated state.
- [x] Set D (998 recovered DEGs) reproduces the paper's biology: skeletal
      system development (p.adj 2.5e-27), extracellular matrix organization,
      embryonic morphogenesis, ossification, pattern specification. That's
      the neural crest / mesoderm program the authors describe, so the
      concordance in 03 is agreement on the biology and not just on gene
      IDs.
- [x] Set C (786 extra DEGs) is essentially null: 3 terms, best p.adj 0.0086,
      against set D's 1e-27. The extras are neither more of the paper's
      biology nor a coherent alternative signal. Their top term is
      miRNA-mediated post-transcriptional gene silencing (26/487) — non-coding
      RNA machinery, consistent with the 34% ncRNA/pseudogene composition of
      the extras and with §5.2's account of where RefSeq and Ensembl v102
      diverge most.
- [x] **GO term overlap between sets C and D is exactly zero.** Recovered:
      256 terms after `simplify()`. Extra: 3. Shared: 0, Jaccard 0.000. So
      it isn't just that the extras enrich more weakly — they share no
      biological theme at all with the genes both analyses agree on. This is
      the strongest evidence I have that they aren't real signal.
      (Raw counts before `simplify()`: 1,150 recovered, 6 extra, 1,249 up,
      2 down.)
- [x] The gap in enrichment strength between up and down (2.4e-23 vs 0.018)
      looks biological, not just statistical. The knockout gains a
      coordinated mesoderm program but loses neural identity diffusely, and
      ORA picks up coordinated programs much better than scattered loss.
      That's probably part of why the authors report a 7.4:1 up/down
      ratio.
- [x] Enrichment parameters, for comparability against their PANTHER run
      (§5.11): GO Biological Process, hypergeometric test, BH correction,
      p.adjust < 0.05, `simplify(cutoff = 0.7)` to collapse redundant terms.
      Universe = genes tested here with a non-NA padj, NOT all annotated
      genes — using the full annotation would inflate every p-value, since
      genes never testable here cannot be enriched. Their universe differs
      again, so GO results diverge for this reason before any biology.

### Extensions
- [x] `independentFiltering = FALSE`, `cooksCutoff = FALSE` (05A): 2,226 DEGs
      vs 2,252 baseline — a fall of 26, closing 3.3% of the 788-gene gap to
      their 1,464. Recovery essentially unchanged (68.0% vs 68.2%).
      **§5.13's filtering account is falsified**: the direction is right but
      the magnitude is negligible. Registered as the primary test of §5.13
      and reported as a failure.

      Confirmed separately: padj NAs went from 2,326 to exactly 0, matching
      Table S3's zero NAs. This verifies §5.6's inference about the authors'
      undocumented parameter choice, which had been made from the absence of
      NAs in their table alone.
- [x] Leave-one-out dropping KO_4: 2,826 DEGs vs 2,252 baseline — removing a
      sample INCREASED the count by 574. Jaccard overlap with the baseline
      DEG list is 0.662, so roughly a third of the union differs. KO_4 is
      the divergent sample on PC2 (§9); including it inflates within-KO
      dispersion and makes every test conservative. Dropping it tightens the
      KO group and hundreds of genes cross threshold.

      **Implication:** the DEG count is not a stable quantity at this sample
      size — it moves 25% on the removal of one of six samples, and the
      68.2% recovery figure inherits that instability. This also bears on
      §5.13: a design this sensitive to composition can produce a 788-gene
      swing from sample membership alone, and this reanalysis differs from
      the authors' by exactly one sample (the missing IC1_1, §5.1). Not
      testable here, but a candidate explanation that does not require any
      pipeline difference at all.

- [x] **KO_4's effect is direction-specific.** Dropping it barely changes
      the up count (1,483 → 1,424) but nearly doubles the down count
      (769 → 1,402), so all 574 genes I gain are downregulated. It looks
      like KO_4 inflates within-KO variance specifically for genes that go
      down in the knockout. Worth noting this runs opposite to the Day 2
      puzzle: without KO_4 my up/down ratio is 1.0:1, which is further from
      the authors' 7.4:1 than my baseline 1.9:1. So keeping KO_4 in actually
      moves me toward their result.
- [x] Exploratory `~ passage_num + condition` (continuous passage, §5.10):
      2,315 DEGs vs 2,252 baseline, Jaccard 0.858. Passage absorbs almost
      nothing, confirming the PCA's finding that no passage gradient
      structures the data. §5.10 is closed: passage is not a meaningful
      unmodeled covariate in this design.
- [x] **Recovery barely moves across all three variants** (68.2 / 68.0 /
      62.6 / 67.3). Nothing I can change downstream of the count matrix
      improves agreement with Table S3, and the variant that changes the DEG
      count most (05B, +25%) makes recovery worse by 5.6 points. Together
      with §5.13 failing, this puts the disagreement upstream, in
      quantification and annotation rather than in my statistical
      settings.

- [x] **Sensitivity summary.** Across the four variants the DEG count runs
      2,226 (filtering off) / 2,252 (baseline) / 2,315 (passage) / 2,826
      (drop KO_4). The two analytical choices I made move it by under 3%,
      but removing one of six samples moves it by 25%. At this sample size,
      which samples are included matters more than how I analyze them. The
      68.2% recovery figure carries that same instability, so I report the
      range alongside it rather than treating it as a fixed number.

      Full table: `results/tables/05_sensitivity_summary.csv`

---

## 10. Repo layout

Matches `git ls-files` as of Day 5. Data directories are gitignored.

```
rnaseq-reanalysis-htt/
├── .gitignore
├── METHODS_LOG.md
├── README.md
├── rnaseq-reanalysis-htt.Rproj
├── scripts/
│   ├── 01_load_data.R           counts + metadata import, sanity checks
│   ├── 02_deseq2.R              DE analysis
│   ├── 03_concordance.R         comparison against Table S3
│   ├── 04_enrichment.R          clusterProfiler GO
│   ├── 05_sensitivity.R         filtering off, leave-one-out, passage
│   └── 06_figures.R             volcano + concordance panel
├── data/                        gitignored
│   ├── raw/                     see Provenance in §3
│   └── processed/               regenerated by the scripts
├── notes/
│   └── prep_python_to_r.R       coding conventions for all scripts
└── results/
    ├── figures/
    │   ├── 02_dispersion.png
    │   ├── 02_ma_comparison.png
    │   ├── 02_pca.png
    │   ├── 03_lfc_scatter.png
    │   ├── 03_recovery_by_basemean.png
    │   ├── 04_go_up.png
    │   ├── 04_go_down.png
    │   ├── 04_go_recovered.png
    │   ├── 04_go_extra.png
    │   ├── 05_deg_counts.png
    │   ├── 05_recovery.png
    │   ├── 06_volcano.png
    │   └── 06_concordance_panel.png
    └── tables/
        ├── 02_results_full.csv
        ├── 03_shared_genes.csv
        ├── 03_missed_genes.csv
        ├── 03_extra_genes.csv
        ├── 04_go_up.csv
        ├── 04_go_down.csv
        ├── 04_go_recovered.csv
        ├── 04_go_extra.csv
        └── 05_sensitivity_summary.csv
```

Note the script numbering runs 01-06 but 05 was written and run before 04:
the sensitivity analysis tested §5.13, which the enrichment work then built
on. Running them in numeric order reproduces everything regardless.

Repo: https://github.com/ritvikK05/rnaseq-reanalysis-htt
