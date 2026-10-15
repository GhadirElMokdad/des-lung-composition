# Desmin (DES) in lung adenocarcinoma: tissue-composition analyses

Scripts and result tables behind the article *Desmin (DES) in Lung Adenocarcinoma: Multi-Database Characterization Points to a Tissue-Composition Origin Rather Than a Tumor-Intrinsic Mechanism*. Archive DOI: **10.5281/zenodo.23398863
**.

## Reproduce
```
pip install -r requirements.txt      # Python 3.11
./run_all.sh                          # downloads about 2 GB of public data, then runs 01 to 08 (about 15 min)
```
`python 00_download_data.py --skip-large` fetches only the small files. Data go to `./data` (change with `DES_DATA=/path`). All scripts write to `./results`. A full re-run on a clean checkout reproduced every result file byte for byte, except `literature_replication.json` (depends on the run date).

## Scripts, outputs and article items
| Script | What it does | Outputs | Article |
|---|---|---|---|
| `00_download_data.py` | downloads all inputs | `data/` | Methods |
| `01_tcga_luad_bulk.py` | extracts genes from the TCGA-LUAD TPM matrix | `data/derived/luad_genes.pkl` | – |
| `02_bulk_strengthen.py` | ROC with DeLong tests, Cox models (HR per SD), cut-offs, fold changes, correlations | `luad_roc_comparators.csv`, `luad_cox_models.csv`, `luad_fold_changes.csv`, `luad_des_correlations_tumour.csv`, `luad_summary.json` | Tables 4–5, Fig. 6B |
| `03a_extract_luca_census.py` | counts of 22 genes of the Lung Cancer Atlas through the Census SOMA API | `data/census/luca_*.pkl` | – |
| `03b_luca_des.py` | DES by cell type, donor-level tests, composition decomposition, origin of the eight retained genes | `luca_*.csv/json` | Fig. 6A,C |
| `03c_krasnow_normal_lung.py` | DES in a normal-lung atlas | `des_krasnow_*.csv` | Supplement |
| `04_figure6.py` | Figure 6 | `Figure6_DES_composition.png` | Fig. 6 |
| `05_cptac_ccle.py` | CPTAC proteomics vs purity; CCLE expression | `cptac_ccle_summary.json` | Results |
| `06_ccle_luad_prad_lusc.py` | CCLE adenocarcinoma lines; TCGA-PRAD and TCGA-LUSC | `ccle_luad_prad_lusc_summary.json` | Results |
| `07_literature_replication.py` | Europe PMC replication of the literature search | `literature_replication.json` | Results |
| `08_spatial_visium.py` | Visium NSCLC sections: DES by spot class, tumour vs adjacent | `spatial_visium_*.csv/json`, `Figure7_spatial_DES.png` | Fig. 7 |

## Data sources (open access, retrieved 2026-10-07)
- **TCGA-LUAD/LUSC/PRAD**, STAR TPM, survival and clinical tables: GDC hub on S3 (UCSC Xena, Goldman et al. 2020), gencode v36 probe map.
- **CPTAC LUAD (Gillette et al. 2020)** and **CCLE (Ghandi et al. 2019)**: cBioPortal DataHub, studies `luad_cptac_2020` and `ccle_broad_2019`.
- **Lung Cancer Atlas, extended (Salcher et al. 2022)**: CZ CELLxGENE Census release 2025-11-04, dataset `1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1`; normal-lung atlas (Krasnow lab) dataset `8c42cfd0-0b0a-46d5-910c-fc833d83c45e`.
- **Visium NSCLC (De Zuani et al. 2024)**: BioStudies E-MTAB-13530 (files are served from `ftp.ebi.ac.uk`).
- **Literature**: Europe PMC REST API.

## Method notes
- A cell is DES-positive if it has at least one raw count; compartments come from the atlas annotation. Donor, not cell, is the unit for tests (at least 30 cells per compartment, 500 for pseudobulk values).
- Survival: one tumour per patient (n = 503, 182 deaths); Cox models use DES as log2(TPM+1) standardised (HR per SD).
- Visium: smooth-muscle spots are the top 10% of the smooth-muscle score; epithelial spots have a score in the top 40% with smooth-muscle and vascular scores at or below the median.

## Limitations
CPTAC, HPA, DepMap, STRING and Enrichr values quoted in the article were not recomputed. The atlas labels only 30,906 tumour cells as malignant. Donor-level stromal proportions did not differ significantly between tumour and normal lung. Visium spots are multicellular (55 µm). TCGA tumour numbers differ slightly from UALCAN (528 vs 515).

## Licence
Code MIT, results CC BY 4.0 (see `LICENSE`). Fill in authors in `LICENSE`, `CITATION.cff` and `.zenodo.json` before release.
