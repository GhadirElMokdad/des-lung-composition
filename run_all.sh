#!/usr/bin/env bash
# Reproduce every number and figure used in the article. Requires Python 3.11, requirements.txt, internet access
# (hosts: gdc-hub.s3.us-east-1.amazonaws.com, media.githubusercontent.com, www.ebi.ac.uk, ftp.ebi.ac.uk, cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com).
set -euo pipefail
python 00_download_data.py            # add --skip-large to fetch only small files
python 01_tcga_luad_bulk.py           # TCGA-LUAD gene table
python 02_bulk_strengthen.py          # ROC, Cox, correlations      -> Tables 4-5, Figure 6B
python 03a_extract_luca_census.py     # LuCA counts via SOMA        (about 3 min)
python 03b_luca_des.py                # DES by cell type            -> Figure 6A,C
python 03c_krasnow_normal_lung.py     # normal-lung atlas
python 04_figure6.py                  # Figure 6
python 05_cptac_ccle.py               # CPTAC proteomics/purity, CCLE
python 06_ccle_luad_prad_lusc.py      # lines, PRAD, LUSC
python 07_literature_replication.py   # Europe PMC counts (results depend on the run date)
python 08_spatial_visium.py           # Visium                      -> Figure 7
