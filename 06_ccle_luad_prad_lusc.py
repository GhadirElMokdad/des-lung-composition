"""(a) DES in CCLE lung adenocarcinoma lines; (b) TCGA-PRAD and TCGA-LUSC: DES fold change and correlation with smooth-muscle/epithelial genes."""
import pandas as pd, numpy as np, gzip, json, sys
from scipy import stats
from common import *
D = CBIO; res = {}
cl = pd.read_csv(f'{D}/ccle_broad_2019_data_clinical_sample.txt', sep='\t', comment='#').set_index('SAMPLE_ID')
want = ['DES', 'MYH11', 'EPCAM', 'SFTPC']; rows = {}
for ch in pd.read_csv(f'{D}/ccle_broad_2019_data_mrna_seq_rpkm.txt', sep='\t', chunksize=2000, index_col=0):
    for g, r in ch[ch.index.isin(want)].iterrows(): rows[g] = r
E_ = pd.DataFrame(rows)
lung = E_[E_.index.str.endswith('_LUNG')].join(cl[['HIST_SUBTYPE1']])
luad = lung[lung.HIST_SUBTYPE1.isin(['adenocarcinoma', 'bronchioloalveolar_adenocarcinoma'])]
res['ccle_luad'] = dict(n=len(luad), median=float(luad.DES.median()), n_gt1=int((luad.DES > 1).sum()), n_gt0p1=int((luad.DES > 0.1).sum()), max=float(luad.DES.max()), epcam_median=float(luad.EPCAM.median()))
rh = E_[E_.index.str.endswith('_SOFT_TISSUE')].DES; res['ccle_soft_tissue_top5_median'] = float(rh.sort_values(ascending=False).head(5).median())
res['ccle_luad_fold_below_top_sarcoma'] = float(rh.sort_values(ascending=False).head(5).median() / max(luad.DES.median(), 1e-3))
# TCGA PRAD / LUSC
pm = pd.read_csv(f'{GDC}/gencode.v36.annotation.gtf.gene.probemap', sep='\t'); mp = dict(zip(pm.id, pm.gene))
G = ['DES', 'MYH11', 'CNN1', 'ACTA2', 'TAGLN', 'EPCAM', 'KLK3', 'NKX3-1', 'KRT5', 'KRT7', 'AGER', 'SFTPC']
for c in ['PRAD', 'LUSC']:
    rows = {}
    with gzip.open(f'{GDC}/TCGA-{c}.star_tpm.tsv.gz', 'rt') as f:
        cols = f.readline().rstrip('\n').split('\t')[1:]
        for line in f:
            k, rest = line.split('\t', 1); g = mp.get(k)
            if g in G and g not in rows: rows[g] = np.array(rest.rstrip('\n').split('\t'), float)
    X = pd.DataFrame(rows, index=cols); code = X.index.str[13:15]; tum = X[code == '01']; nor = X[code == '11']; tpm = lambda x: 2 ** x - 1
    r = dict(n_tumour=len(tum), n_normal=len(nor), median_tpm_tumour=float(tpm(tum.DES).median()), median_tpm_normal=float(tpm(nor.DES).median()),
             fold_drop=float(tpm(nor.DES).median() / tpm(tum.DES).median()), log2fc=float(np.log2(tpm(tum.DES).median() + 1) - np.log2(tpm(nor.DES).median() + 1)))
    y = np.r_[np.ones(len(tum)), np.zeros(len(nor))]; s = -np.r_[tum.DES, nor.DES]; r['auc_DES'] = float(stats.mannwhitneyu(s[y == 1], s[y == 0]).statistic / (len(tum) * len(nor)))
    for g in ['MYH11', 'CNN1', 'ACTA2', 'TAGLN', 'EPCAM'] + (['KLK3', 'NKX3-1'] if c == 'PRAD' else ['KRT5']):
        if g in X: r[f'fold_drop_{g}'] = float(tpm(nor[g]).median() / max(tpm(tum[g]).median(), 1e-9)); r[f'rho_DES_{g}_tumour'] = float(stats.spearmanr(tum.DES, tum[g])[0]); r[f'rho_DES_{g}_normal'] = float(stats.spearmanr(nor.DES, nor[g])[0])
    res[c] = r
json.dump(res, open(f'{OUT}/ccle_luad_prad_lusc_summary.json', 'w'), indent=1); print(json.dumps(res, indent=1))
