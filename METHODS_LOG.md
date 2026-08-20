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

Rationale: larger effect size than the HD comparison (1,401 vs 822 DEGs), and
it is the paper's cleanest loss-of-function contrast.

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

**Counts** — NCBI-generated raw count matrix from GEO ("Download RNA-seq counts").
- Identifiers: **Entrez Gene IDs**
- Verified integer: `all(counts == round(counts))` → TRUE
- 10 sample columns (GSM8343537 missing)

**Truth table** — Supplementary Table 3 (`13578_2025_1443_MOESM5_ESM.xlsx`)
- Sheet `KO-NSC vs IC1-NSC` → **`skip = 1`**
- Sheet `HD-NSC vs IC1-NSC` → **`skip = 2`** (extra title row — different offset!)
- Columns: `gene_ID`, `Gene_name`, `log2FoldChange`, `Fold_change`, `P-value`,
  `P_adjusted`, `Description`
- Identifiers: **Ensembl, no version suffix**
- ~31,000 rows (KO) / ~28,000 rows (HD) — full results, not just significant genes

**Not used:** `GSE270472_HD_KO_NSC.xlsx` contains TPM values only. TPM is
normalized and fractional — DESeq2 requires raw integer counts. Retained only
for eyeballing individual genes.

---

## 4. Original authors' parameters (match these)

- **Thresholds:** |log2FC| > 1.5 AND padj < 0.05
  - Note: 1.5, *not* 1. Stricter than the common default.
- **Reported DEGs:** 1,401 (KO vs IC1), 822 (HD vs IC1); 331 overlap
  - KO breakdown: 1,231 up / 170 down
  - HD breakdown: 387 up / 435 down
- **GO enrichment:** PANTHER v17, Fisher's exact test, FDR correction p < 0.05
- **Bioinformatics vendor:** IDEAS4BIOLOGY (external company)
- **Pipeline details:** ⚠️ TO CHECK — in Supplementary Methods, not main text.
  Need: aligner, quantifier, annotation release, DE package.

**Validation check:** filtering their table to the above thresholds should
return ~1,401 rows. If it doesn't, stop and resolve before proceeding.

---

## 5. Known discrepancy sources

The core of the writeup. Each is a candidate explanation for imperfect concordance.

1. **Missing control replicate.** NCBI's pipeline excluded GSM8343537, leaving
   n=2 controls vs the authors' n=3. Reduces power; expect to miss borderline genes.

2. **Different quantification pipeline.** NCBI's standardized counts vs the
   vendor's pipeline. Different aligner, annotation release, and counting rules.

3. **Shrinkage.** Using `lfcShrink(type="apeglm")`; their table contains
   unshrunken estimates. Evidence: top KO hit FKBPL has log2FC ≈ 22
   (fold change ~4.1 million) — a near-zero-denominator artifact that apeglm
   will pull hard toward zero. **Expect substantially different top-ranked genes
   even where significance calls agree.**

4. **ID mapping loss.** Counts use Entrez, truth table uses Ensembl.
   `AnnotationDbi::select()` returns 1:many mappings. Record the final mapping
   rate and how duplicates were resolved.

5. **Control line caveat (authors' own).** IC1 was originally intended as the
   isogenic control, but the corrected allele turned out to be silenced,
   producing monoallelic *HTT* expression. The authors switched to IC2 for later
   experiments — but the RNA-seq comparison still uses IC1. So the "control" is
   not a clean wild-type.

6. **Passage as unmodeled covariate.** Passages are matched across lines
   (IC1: 4–6; KO: 4–7), so passage functions like a batch variable.
   Baseline model is `~ condition`; test `~ passage + condition` as an extension.

7. **Enrichment tool.** clusterProfiler vs their PANTHER v17. Different gene
   universes and algorithms — GO results will differ for tool reasons alone.

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

Expected upregulated in KO (validated by the authors via RT-qPCR):
`TWIST1`, `SIX1`, `TBX1`, `TBX15`, `MSX2`, `MEOX2`, `FOXD1`

Expected downregulated: `PAX6` (paper reports log2FC = −2.38, KO vs IC1)

If these don't appear with the right direction, check factor levels before
anything else.

---

## 8. Environment

- macOS (Apple Silicon), R 4.6.1, Bioconductor 3.23
- DESeq2, apeglm, EnhancedVolcano, clusterProfiler 4.20.0, org.Hs.eg.db, tximport
- Capture `sessionInfo()` output for the README
- Lockfile: `renv::snapshot()` → commit `renv.lock`

**Fallback rule:** if still fighting R syntax at end of Day 2, switch to
pydeseq2. No exceptions.

---

## 9. Results to fill in

- [ ] Their DEG count reproduced from Table S3: _____ (target ~1,401)
- [ ] My DEG count at matched thresholds: _____
- [ ] Entrez → Ensembl mapping rate: _____%
- [ ] Genes recovered (shared): _____
- [ ] Missed (theirs, not mine): _____
- [ ] Extra (mine, not theirs): _____
- [ ] Recovery rate: _____%
- [ ] Of missed genes, fraction with padj between 0.05 and 0.15: _____
      (threshold sensitivity vs genuine disagreement)
- [ ] Spearman correlation of log2FC across all shared genes: _____
- [ ] Did `~ passage + condition` change the result? _____
