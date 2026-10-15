"""Download every public input (about 2 GB). Usage: python 00_download_data.py [--skip-large] [--dry]
Sources (all open access, retrieved 2026-10-07):
  GDC hub (UCSC Xena, S3):  TCGA-LUAD/LUSC/PRAD STAR TPM, TCGA-LUAD survival and clinical tables, gencode v36 probe map
  cBioPortal DataHub (GitHub, Git LFS): CPTAC LUAD (luad_cptac_2020) and CCLE (ccle_broad_2019)
  BioStudies E-MTAB-13530: 10x Visium NSCLC sections (files are served from ftp.ebi.ac.uk)
  CELLxGENE Census S3: Krasnow-lab normal-lung atlas (h5ad). The Lung Cancer Atlas counts are read directly by 03a_extract_luca_census.py."""
import sys, os, csv, json, subprocess, urllib.request
from common import *
DRY = '--dry' in sys.argv; SKIP = '--skip-large' in sys.argv
def get(url, dest, large=False):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000: print('have', dest); return
    if large and SKIP: print('skip (large)', dest); return
    print('get ', url, '->', dest)
    if DRY: return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    subprocess.run(['curl', '-sSL', '--fail', '--retry', '3', '-o', dest, url], check=True)
G = 'https://gdc-hub.s3.us-east-1.amazonaws.com/download'
for f, big in [('TCGA-LUAD.star_tpm.tsv.gz', True), ('TCGA-LUAD.survival.tsv.gz', False), ('TCGA-LUAD.clinical.tsv.gz', False), ('TCGA-LUSC.star_tpm.tsv.gz', True), ('TCGA-PRAD.star_tpm.tsv.gz', True), ('gencode.v36.annotation.gtf.gene.probemap', False)]:
    get(f'{G}/{f}', f'{GDC}/{f}', big)
M = 'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public'
for study, fs in [('luad_cptac_2020', ['data_clinical_sample.txt', 'data_protein_quantification.txt']), ('ccle_broad_2019', ['data_clinical_sample.txt', 'data_mrna_seq_rpkm.txt'])]:
    for f in fs: get(f'{M}/{study}/{f}', f'{CBIO}/{study}_{f}', f == 'data_mrna_seq_rpkm.txt')
get('https://cellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com/cell-census/2025-11-04/h5ads/8c42cfd0-0b0a-46d5-910c-fc833d83c45e.h5ad', f'{CENSUS}/krasnow10x.h5ad', True)
# Visium: file list from the BioStudies API, then every filtered matrix and spatial archive
B = 'https://www.ebi.ac.uk/biostudies'
if not DRY:
    st = json.load(urllib.request.urlopen(f'{B}/api/v1/studies/E-MTAB-13530', timeout=60)); files = []
    def rec(s):
        for f in s.get('files', []):
            for g in (f if isinstance(f, list) else [f]): files.append(g['path'])
        for ss in s.get('subsections', []):
            for x in (ss if isinstance(ss, list) else [ss]): rec(x)
    rec(st['section'])
    for f in files:
        if f.endswith('filtered_feature_bc_matrix.h5') or f.endswith('spatial.tar'): get(f'{B}/files/E-MTAB-13530/{f}', f'{VISIUM}/{f}', True)
    get(f'{B}/files/E-MTAB-13530/E-MTAB-13530.sdrf.txt', f'{VISIUM}/sdrf.txt')
    rows = list(csv.reader(open(f'{VISIUM}/sdrf.txt'), delimiter='\t')); h = rows[0]; ix = {k: i for i, k in enumerate(h)}
    meta = {r[0]: [r[ix['Characteristics[individual]']], r[ix['Characteristics[disease]']]] for r in rows[1:]}
    json.dump(meta, open(f'{VISIUM}/samples.json', 'w'), indent=1); print('samples.json:', len(meta), 'sections')
