"""Counts of 22 genes for the extended Lung Cancer Atlas (Salcher et al. 2022; CELLxGENE Census 2025-11-04, dataset 1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1)
read through the SOMA API (anonymous S3 access). Output: data/census/luca_{obs,var,counts}.pkl. Takes about 3 minutes."""
import os, tiledbsoma as t, pandas as pd, numpy as np, time
from common import *
os.makedirs(CENSUS, exist_ok=True)
ctx = t.SOMATileDBContext(tiledb_config={"vfs.s3.region": "us-west-2", "vfs.s3.no_sign_request": "true"})
exp = t.open("s3://cellxgene-census-public-us-west-2/cell-census/2025-11-04/soma/census_data/homo_sapiens/", context=ctx)
did = '1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1'
GENES = ['DES', 'ACTA2', 'MYH11', 'CNN1', 'TAGLN', 'PDGFRB', 'RGS5', 'CDH5', 'PECAM1', 'SFTPC', 'AGER', 'NAPSA', 'EPCAM', 'KRT7', 'CA4', 'ADH1B', 'BTNL9', 'SPP1', 'FAM83A', 'PTPRC', 'COL1A1', 'CLDN18']
t0 = time.time()
obs = exp.obs.read(value_filter=f"dataset_id == '{did}'", column_names=['soma_joinid', 'donor_id', 'cell_type', 'disease', 'tissue', 'tissue_general', 'raw_sum', 'assay', 'suspension_type']).concat().to_pandas()
obs.to_pickle(f'{CENSUS}/luca_obs.pkl'); print(len(obs), time.time() - t0)
var = exp.ms['RNA'].var.read(value_filter="feature_name in %s" % GENES, column_names=['soma_joinid', 'feature_name', 'feature_id']).concat().to_pandas(); var.to_pickle(f'{CENSUS}/luca_var.pkl')
X = exp.ms['RNA'].X['raw']; oid = obs.soma_joinid.values; vids = var.soma_joinid.values; rows, cols, vals = [], [], []
for i in range(0, len(oid), 60000):
    r = X.read((np.sort(oid[i:i + 60000]), vids)).tables().concat().to_pandas(); rows.append(r.soma_dim_0.values); cols.append(r.soma_dim_1.values); vals.append(r.soma_data.values); print(i, round(time.time() - t0), flush=True)
pd.DataFrame({'obs': np.concatenate(rows), 'var': np.concatenate(cols), 'v': np.concatenate(vals)}).to_pickle(f'{CENSUS}/luca_counts.pkl')
