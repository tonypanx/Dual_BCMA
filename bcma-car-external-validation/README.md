# External validation of the dual-BCMA CAR T case in the Rade et al. atlas

Analysis code for Figure 4d-h and Supplementary Figure S15.

The manuscript reports a single patient treated with dual BCMA-targeted CAR T
cells. To test whether the phenotypes seen in that case also separate cilta-cel
from ide-cel at the cohort level, the same signatures and gene panels were run
against the published single-cell atlas of Rade et al. (2026).

The atlas profiles all post-infusion T cells rather than sorted CAR T cells, so
every analysis here is restricted to CAR-transgene-positive cells first, using the
atlas authors' own `CAR_BY_EXPRS` call. Comparisons are therefore CAR+ versus 
CAR+ throughout.

Sample sizes after gating: 25 ide-cel and 18 cilta-cel patients with
post-infusion CAR+ cells; 17 ide and 7 cilta patients clear the per-patient
cell-count threshold for the CD8 tests, and 7 ide and 6 cilta for CD4.

## Where each panel comes from

| Panel | Script | Output file |
|---|---|---|
| Fig. 4e | `02_gsea_and_module_scores.R` | `results/gsea/gsea_CD8_cilta_vs_ide_enrichment.pdf`, `gsea_CD4_cilta_vs_ide_enrichment.pdf` |
| Fig. 4f | `02_gsea_and_module_scores.R` | `results/gsea/module_scores_selected.pdf` |
| Fig. 4g | `03_regulon_pseudobulk.R` | `results/regulon/regulon_pseudobulk_min10.pdf` |
| Fig. 4h | `05_key_gene_heatmap.R` | `results/heatmap/key_genes_avg_4group.pdf` |
| Fig. S15a | `02_gsea_and_module_scores.R` | `results/gsea/gsea_CD8_cilta_vs_ide_enrichment.pdf` |
| Fig. S15b | `04_key_gene_violins.R` | `results/key_genes/key_genes_pb_cilta_vs_ide_CD8.pdf` |
| Fig. S15c | `02_gsea_and_module_scores.R` | `results/gsea/gsea_CD4_cilta_vs_ide_enrichment.pdf` |
| Fig. S15d | `04_key_gene_violins.R` | `results/key_genes/key_genes_pb_cilta_vs_ide_CD4.pdf` |

The scripts write every panel that clears significance, or every gene in a
panel; the published figures are the relevant subsets, cropped and arranged in
Illustrator. Specifically:

- Fig. 4e is the Kaech naive-vs-day-8-effector curve from each of the two
  enrichment PDFs.
- Fig. S15a and S15c are the RUNX3 regulon, KLF2 regulon and Hallmark IFN-gamma
  curves from the same two files, plus Hallmark inflammatory response for CD4.
- Fig. 4f is the two CD8 panels of `module_scores_selected.pdf`.
- Fig. S15b is *GNLY, NKG7, BAX, CD27, TNFRSF9* and *FYN* from the CD8 violin
  sheet; Fig. S15d is *TNFRSF18, CX3CR1, BAX, FOXO1* and *ICOS* from the CD4 one.

## Method notes

**Pseudobulk everywhere.** The cohort is 43 patients, so all tests are on
per-patient means rather than per-cell values. Treating cells as independent
replicates would inflate significance by two orders of magnitude here.

**Per-patient cell minimum.** A patient must contribute at least ten CAR+ cells
to a comparison before their mean is used (`CAR_MIN` in `config.R`). CAR+ cells
are sparse at the later timepoints, and means from a handful of cells are noisy
enough to swing a Wilcoxon on this sample size. The threshold is applied to each
comparison separately, so a patient can qualify for the CD8 test and not the CD4
one. `03_regulon_pseudobulk.R` repeats its tests at a threshold of twenty as a
sensitivity check; the direction and significance of Fig. 4g are unchanged.

**GSEA ranking.** Genes are ranked by `sign(log2FC) * -log10(p)` computed on the
per-patient means, then passed to `fgseaMultilevel` with `eps = 0`. The earlier
`fgseaSimple` runs put a `1/(nperm+1)` floor under the p-values, which left
several sets stacked on the BH boundary.

**Module scores versus GSEA.** The curated gene sets are eight to twelve genes
each; too small for a stable rank-based test. Those are scored per cell with
`AddModuleScore`, averaged per patient, and compared with a two-sided Wilcoxon.
The larger signature and MSigDB sets go through GSEA. Sets smaller than ten
genes are dropped from GSEA entirely.

**Multiple testing.** GSEA p-values are BH-adjusted within each comparison.
Module scores are BH-adjusted within lineage. The per-gene violin figures
display raw p-values, as noted in their subtitles; adjusted values for the same
tests are in `results/key_genes/key_genes_pb_stats.csv`.

**Responder definition.** Best response of CR or VGPR counts as a responder.
Response-stratified comparisons are computed and written out but do not appear
in Figure 4 or S15.

## Requirements

R 4.3 or later.

CRAN: `Seurat`, `ggplot2`, `dplyr`, `patchwork`, `scales`, `msigdbr`,
`showtext` (optional, used only to embed Arial).

Bioconductor: `fgsea`.

```r
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork", "scales",
                   "msigdbr", "showtext"))
BiocManager::install("fgsea")
```

Random seeds are fixed at 42 for `AddModuleScore` and `fgseaMultilevel`, so
repeated runs on the same input give identical numbers.

## Layout

```
config.R                      paths and the per-patient cell minimum
R/common.R                    column names, cohort definitions, small helpers
R/car_gating.R                CAR+ gating and the per-patient cell filter
00_export_car_barcodes.R
01_signature_scoring.R
02_gsea_and_module_scores.R
03_regulon_pseudobulk.R
04_key_gene_violins.R
05_key_gene_heatmap.R
```
