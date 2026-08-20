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

**Comparison chosen:** KO-NSC vs IC1-NSC (HTT knockout vs isogenic control)

Rationale: larger effect size than the HD comparison (1,401 vs 822 DEGs as
reported; 1,464 vs 835 as recomputed — see §4), and it is the paper's cleanest
loss-of-function contrast.

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
- Verified integer: `all(counts == round(counts))` → TRUE
- 10 sample columns (GSM8343537 missing)

### Annotation
`Human_GRCh38_p13_annot.tsv.gz` (GEO, GRCh38.p13)

- 39,376 rows
- Columns: `GeneID` (Entrez), `Symbol`, `Description`, `Synonyms`, `GeneType`,
  `EnsemblGeneID`, `Status`, coordinates, GO term columns
- `EnsemblGeneID` empty for 12,299 genes (31.2%) — composition in §9
- Note: blanks import as `""`, not `NA`. Normalize on load:
  `annot$EnsemblGeneID[annot$EnsemblGeneID == ""] <- NA`

### Truth table
Supplementary Table 3 (`13578_2025_1443_MOESM5_ESM.xlsx`)

- Sheet `KO-NSC vs IC1-NSC` → **`skip = 1`** — 31,264 rows
- Sheet `HD-NSC vs IC1-NSC` → **`skip = 2`** — ~28,057 rows
  (extra title row — the offset differs between sheets!)
- Columns: `gene_ID`, `Gene_name`, `log2FoldChange`, `Fold_change`, `P-value`,
  `P_adjusted`, `Description`
- Identifiers: **Ensembl, no version suffix**
- Full results, not just significant genes

### Not used
`GSE270472_HD_KO_NSC.xlsx` contains TPM values only. TPM is normalized and
fractional — DESeq2 requires raw integer counts. Retained only for eyeballing
individual genes.

### Provenance
- Counts file: `GSE270472_raw_counts_GRCh38.p13_NCBI.tsv.gz`
  ⚠️ confirm exact patch version from the filename on disk
- Annotation file: `Human_GRCh38_p13_annot.tsv.gz`
- Both downloaded 2026-08-20 from GEO GSE270472 → "Download RNA-seq counts"
- Supplementary Table 3 downloaded 2026-08-20 from
  https://link.springer.com/article/10.1186/s13578-025-01443-5

---

## 4. Original authors' parameters

### As stated in the paper
- **Thresholds:** |log2FC| > 1.5 AND padj < 0.05
  (note: 1.5, *not* 1 — stricter than the common default)
- **Reported DEGs:** 1,401 (KO vs IC1), 822 (HD vs IC1); 331 overlap
  - KO: 1,231 up / 170 down
  - HD: 387 up / 435 down
- **GO enrichment:** PANTHER v17, Fisher's exact test, FDR correction p < 0.05
- **Bioinformatics vendor:** IDEAS4BIOLOGY (external company)
- **Pipeline details:** ⚠️ TO CHECK — Supplementary Methods, not main text.
  Need: aligner, quantifier, annotation release, DE package.

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

### Inferred DE framework
Supplementary Table 3 reports adjusted p-values for all 31,264 genes with
**zero NAs**. DESeq2's default independent filtering and Cook's distance
outlier detection both produce NAs, so the original analysis likely used
edgeR, limma-voom, or DESeq2 with both filters disabled. This is a structural
difference from the reanalysis and is expected to account for part of the
non-concordance. ⚠️ Confirm against Supplementary Methods.

---

## 5. Known discrepancy sources

The core of the writeup. Each is a candidate explanation for imperfect concordance.

1. **Missing control replicate.** NCBI's pipeline excluded GSM8343537, leaving
   n=2 controls vs the authors' n=3. Reduces power; expect to miss borderline
   genes.

2. **Different quantification pipeline.** NCBI's standardized counts vs the
   vendor's pipeline. Different aligner, annotation release, and counting rules.

3. **Shrinkage.** Using `lfcShrink(type="apeglm")`; their table contains
   unshrunken estimates. Evidence: top KO hit FKBPL has log2FC ≈ 22
   (fold change ~4.1 million) — a near-zero-denominator artifact that apeglm
   will pull hard toward zero. **Expect substantially different top-ranked
   genes even where significance calls agree.**

4. **ID mapping loss.** Counts use Entrez, truth table uses Ensembl.
   `AnnotationDbi::select()` returns 1:many mappings; record the final mapping
   rate and how duplicates were resolved. Loss is concentrated in non-coding
   and pseudogene classes — only 178 protein-coding genes lack an Ensembl ID,
   so impact on a protein-coding-dominated DEG list should be minimal.

5. **Control line caveat (authors' own).** IC1 was originally intended as the
   isogenic control, but the corrected allele turned out to be silenced,
   producing monoallelic *HTT* expression. The authors switched to IC2 for
   later experiments — but the RNA-seq comparison still uses IC1. So the
   "control" is not a clean wild-type.

6. **Passage as unmodeled covariate.** Passages are matched across lines
   (IC1: 4–6; KO: 4–7), so passage functions like a batch variable.
   Baseline model is `~ condition`; test `~ passage + condition` as an extension.

7. **Enrichment tool.** clusterProfiler vs their PANTHER v17. Different gene
   universes and algorithms — GO results will differ for tool reasons alone.

8. **Different DE framework.** Their table has zero NAs across all 31,264
   genes, inconsistent with DESeq2 defaults. See §4. Genes DESeq2 assigns NA
   are unrecoverable by construction — quantify this in §9.

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

---

## 7. Sanity checks

### Primary positive control
**HTT** (Entrez 3064, ENSG00000197386) should show strong downregulation in
KO vs IC1. If HTT is not among the most significant downregulated genes,
something is wrong with sample labeling or the count matrix — stop and
investigate before anything else.

Note: CRISPR knockouts often still produce transcript (frameshifted, subject
to nonsense-mediated decay), so expect reduced rather than zero counts.

### Secondary controls
Expected **upregulated** in KO (validated by the authors via RT-qPCR):
`TWIST1`, `SIX1`, `TBX1`, `TBX15`, `MSX2`, `MEOX2`, `FOXD1`

Expected **downregulated**: `PAX6` (paper reports log2FC = −2.38, KO vs IC1)

If these don't appear with the right direction, check factor levels first.

---

## 8. Environment

- macOS (Apple Silicon), R 4.6.1, Bioconductor 3.23
- DESeq2, apeglm, EnhancedVolcano, clusterProfiler 4.20.0, org.Hs.eg.db, tximport
- Capture `sessionInfo()` output for the README
- Lockfile: `renv::snapshot()` → commit `renv.lock`

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
      20 other, 19 rRNA, 28 snRNA, 5 unknown, and only **178 protein-coding**

### Mapping
- [ ] Coverage of count-matrix genes via GEO annotation: _____%
- [ ] Additional genes recovered via `org.Hs.eg.db`: _____
- [ ] Duplicate resolution method: _____
- [ ] Reference DEGs unmappable by construction: _____

### My analysis
- [ ] Genes tested after pre-filtering: _____ (vs 31,264 in their table)
- [ ] My DEG count at matched thresholds: _____ (____ up / ____ down)
- [ ] HTT log2FC and padj: _____ (primary positive control)

### Concordance
- [ ] Recovered (shared): _____
- [ ] Missed (theirs, not mine): _____
- [ ] Extra (mine, not theirs): _____
- [ ] Recovery rate: _____%
- [ ] Of shared DEGs, fraction with matching log2FC sign: _____%
      (should be ~100%; below ~95% suggests a labeling or reference-level problem)
- [ ] Of missed genes, fraction assigned NA by DESeq2 filtering: _____
- [ ] Of missed genes, fraction with padj between 0.05 and 0.15: _____
      (threshold sensitivity vs genuine disagreement)
- [ ] Spearman correlation of log2FC across all shared genes: _____

### Extensions
- [ ] Did `~ passage + condition` change the result? _____

---

## 10. Repo layout

```
rnaseq-reanalysis-htt/
├── METHODS_LOG.md
├── README.md                    (Days 4–5)
├── scripts/
│   ├── 01_load_data.R           counts + metadata import, sanity checks
│   ├── 02_deseq2.R              DE analysis
│   ├── 03_concordance.R         comparison against Table S3
│   └── 04_enrichment.R          clusterProfiler GO
├── data/raw/                    gitignored — see Provenance in §3
├── notes/
└── results/
    ├── figures/
    └── tables/
```

Repo: https://github.com/ritvikK05/rnaseq-reanalysis-htt
