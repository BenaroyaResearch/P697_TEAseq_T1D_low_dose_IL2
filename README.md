# P697 TEA-seq T1D Low-Dose IL-2

This repository contains the analysis workflow for a pilot TEA-seq study of NK and Treg cells from people with T1D sampled longitudinally during a clinical study. The main biological goal is to test whether low-dose IL-2 and subsequent rapamycin exposure produce detectable RNA, ATAC, and multimodal state changes, and whether those changes resolve back toward baseline or persist at followup.

Important context:

- This is TEA-seq with RNA plus ATAC plus Feature Barcoding hashtags for demultiplexing.
- There is no CITE-seq antibody panel in this project.
- The main longitudinal phases used downstream are Baseline, Post-IL2, Post-RAPA, and Followup.

## Repository Layout

- `code/`: analysis scripts, notebooks, helper functions, and environment notes.
- `data/inputData/`: raw or semi-processed inputs used by the pipeline.
- `data/outputData/`: bundle files, intermediate exports, motif or trajectory inputs, and derived analysis objects.
- `figures/`: saved figures from preprocessing, WNN, differential analysis, and PHLOWER.
- `tables/`: exported differential results, enrichment tables, and summary CSVs.
- `documents/`: notes, drafts, presentations, and supporting project documents.

## Analysis Flow

The core workflow is bundle-based. Each major script saves a dated bundle, and the downstream scripts load the newest matching bundle automatically.

1. `code/P697_preprocessing.rmd`
   - Starts from the merged TEA-seq data and sample annotations.
   - Performs RNA and ATAC QC, hashtag-based donor and visit parsing, metadata joins, and CellTypist label integration.
   - Produces donor-by-visit coverage summaries and the preprocessing bundle.
   - Output: `P697.<date>_preprocessingBundle.rds`

2. `code/P697_NKPeakCalling.rmd`
   - Loads the preprocessing bundle.
   - Focuses on the NK compartment.
   - Applies ATAC depth filtering, performs MACS2 peak calling, rebuilds the ATAC assay, computes LSI and Harmony integration, and saves an NK-specific peak-calling bundle.
   - Output: `P697.<date>_NK_peakCallingBundle.rds`

3. `code/P697_TregPeakCalling.rmd`
   - Same structure as the NK peak-calling script, but for the Treg compartment.
   - Output: `P697.<date>_Treg_peakCallingBundle.rds`

4. `code/P697_WNN.rmd`
   - Loads the newest NK and Treg peak-calling bundles.
   - Creates compartment-specific WNN objects using RNA Harmony plus ATAC Harmony-LSI embeddings.
   - Evaluates clustering resolution, defines RNA-only, ATAC-only, and WNN clusters, and generates multimodal cluster characterization plots.
   - Output: `P697.<date>_WNN_bundle.rds`

5. `code/P697_DA.rmd`
   - Loads the newest WNN bundle.
   - Runs donor-paired pseudobulk RNA differential expression and ATAC differential accessibility globally and per WNN cluster.
   - Runs Hallmark, GO:BP, C7, and IL2-STAT5 roast analyses.
   - Exports SCENIC+/PHLOWER inputs and saves a dated DA bundle.
   - Also writes PI-facing coverage summaries, annotated ATAC tables, and multimodal peak-gene summaries.
   - Output: `P697.<date>_DA_bundle.rds`

6. `code/P697_phlower_NK.ipynb` and `code/P697_phlower_Treg.ipynb`
   - Consume the `scenic_export_*` directories written by `code/P697_DA.rmd`.
   - Run PHLOWER trajectory analyses on the combined RNA plus ATAC embedding.
   - These notebooks are best treated as supporting multimodal trajectory analyses rather than the primary source of treatment conclusions.

## Script Reference

### Core R Markdown scripts

| Script | Role | Reads | Writes |
|---|---|---|---|
| `code/P697_preprocessing.rmd` | Initial QC, metadata harmonization, CellTypist integration, donor and visit summaries | raw inputs, CellTypist CSVs | preprocessing bundle, QC figures, cell-count tables |
| `code/P697_NKPeakCalling.rmd` | NK ATAC peak calling and ATAC assay reconstruction | preprocessing bundle | NK peak-calling bundle, NK ATAC QC figures |
| `code/P697_TregPeakCalling.rmd` | Treg ATAC peak calling and ATAC assay reconstruction | preprocessing bundle | Treg peak-calling bundle, Treg ATAC QC figures |
| `code/P697_WNN.rmd` | Multimodal integration and cluster characterization | NK and Treg peak-calling bundles | WNN bundle, WNN figures |
| `code/P697_DA.rmd` | RNA and ATAC differential analysis, enrichment, multimodal summaries, exports | WNN bundle | DA bundle, differential tables, enrichment tables, multimodal summaries, scenic exports |

### Shared code and supporting files

| File | Purpose |
|---|---|
| `code/P697_functions.R` | Shared helpers such as `savePlot()`, dotplot utilities, and clustering optimization helpers |
| `code/P697_preprocessing_initialQCArchived.rmd` | Archived earlier preprocessing and QC workflow; not the active pipeline |
| `code/PHLOWER_PATCHES.md` | Notes about PHLOWER compatibility patches used in this project |
| `code/phlower_environment.yml` | Conda environment definition for the PHLOWER notebooks |

### Notebooks used for annotations or secondary analyses

| Notebook | Purpose | Used by |
|---|---|---|
| `code/CellTypistL1OnP697.ipynb` | CellTypist L1 annotation generation | `code/P697_preprocessing.rmd` reads resulting CSVs |
| `code/CellTypistL2OnP697.ipynb` | CellTypist L2 annotation generation | `code/P697_preprocessing.rmd` reads resulting CSVs |
| `code/CellTypistL3OnP697.ipynb` | CellTypist L3 annotation generation | `code/P697_preprocessing.rmd` reads resulting CSVs |
| `code/P697_phlower_NK.ipynb` | NK trajectory analysis on scenic export | supporting analysis |
| `code/P697_phlower_Treg.ipynb` | Treg trajectory analysis on scenic export | supporting analysis |

## Fastest Rerun Strategy

If the goal is to refresh the downstream figures and tables for the PI without recomputing peak calling, the current bundle structure supports that.

### Fast PI refresh

- Re-run `code/P697_DA.rmd` only.
- This script loads the newest `P697.*_WNN_bundle.rds` automatically.
- As long as the existing WNN bundle is acceptable, you do not need to rerun peak calling.

### Re-run WNN but not peak calling

- Re-run `code/P697_WNN.rmd` if you change clustering, modality weighting, marker characterization, or WNN-based plots.
- This still reuses the newest NK and Treg peak-calling bundles.

### Re-run peak calling only if these change

- ATAC depth threshold.
- MACS2 peak set.
- Rebuilt ATAC assay.
- ATAC LSI or Harmony integration choices.
- Any change to the underlying ATAC feature space.

If none of those change, the existing peak-calling bundles remain valid for downstream DA, enrichment, annotation, and multimodal summary updates.

## Key Outputs

### Bundles

- `data/outputData/P697.<date>_preprocessingBundle.rds`
- `data/outputData/P697.<date>_NK_peakCallingBundle.rds`
- `data/outputData/P697.<date>_Treg_peakCallingBundle.rds`
- `data/outputData/P697.<date>_WNN_bundle.rds`
- `data/outputData/P697.<date>_DA_bundle.rds`

### Common figure groups

- Preprocessing QC: `figures/P697.<date>_QC*`, `figures/P697.<date>_ATAC_*`
- WNN multimodal summaries: `figures/P697.<date>_*_WNN_*`
- Differential analysis: `figures/P697.<date>_*_RNA_DGE_*`, `figures/P697.<date>_*_ATAC_DA_*`
- Enrichment: `figures/P697.<date>_*_GSEA_*`, `figures/roast_IL2STAT5/`
- PHLOWER: `figures/phlower/`

### Common table groups

- RNA differential expression: `tables/P697.<date>_*_RNA_DGE_*.csv`
- ATAC differential accessibility: `tables/P697.<date>_*_ATAC_DA_*.csv`
- Enrichment outputs: `tables/P697.<date>_*_GSEA_*.csv`
- Roast summaries: `tables/P697.<date>_roast_IL2STAT5_*.csv`
- Coverage and multimodal summaries: `tables/P697.<date>_*_coverage*.csv`, `tables/P697.<date>_*_multimodal_peak_gene_*.csv`

## Notes on Interpretation

- RNA, ATAC, and WNN conclusions should be interpreted together.
- WNN is useful for state structure and cluster composition.
- RNA is strongest for immediate biological interpretation and enrichment.
- ATAC is strongest for durable remodeling hypotheses, especially when linked to nearby genes and motifs.
- PHLOWER is best used as supporting evidence for longitudinal state movement, not as the only basis for a treatment conclusion.

## Practical Recommendation

For a time-constrained PI update, start from the newest WNN bundle and run only `code/P697_DA.rmd`. That path preserves the existing peak calls and WNN state definitions while refreshing the downstream RNA, ATAC, and multimodal summaries.