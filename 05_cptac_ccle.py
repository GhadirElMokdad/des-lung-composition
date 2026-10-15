"""CPTAC-LUAD proteomics vs tumour purity, and CCLE expression of DES in lung adenocarcinoma lines (cBioPortal DataHub, LFS)."""
import pandas as pd, numpy as np, json, sys
from scipy import stats
from common import *
D = sys.argv[1] if len(sys.argv) > 1 else CBIO
# ---------- CPTAC LUAD: protein vs purity (110 tumours)
c = pd.read_csv(f'{D}/luad_cptac_2020_data_clinical_sample.txt', sep='\t', comment='#').set_index('PATIENT_ID')
p = pd.read_csv(f'{D}/luad_cptac_2020_data_protein_quantification.txt', sep='\t'); p.index = p.iloc[:, 0].str.split('|').str[0]; p = p.drop(columns=p.columns[0]).T
G = ['DES', 'MYH11', 'CNN1', 'ACTA2', 'TAGLN', 'AGER', 'SFTPC', 'NAPSA', 'PECAM1', 'CDH5', 'EPCAM', 'VIM']
X = p[G].join(c[['TUMOR_PURITY_BYESTIMATE_RNASEQ', 'TSNET PURITY', 'ESTIMATE STROMALSCORE', 'ESTIMATE IMMUNESCORE']], how='inner')
res = {'n': len(X)}
for pur in ['TUMOR_PURITY_BYESTIMATE_RNASEQ', 'TSNET PURITY', 'ESTIMATE STROMALSCORE']:
    d = X[['DES', pur]].dropna(); r, pv = stats.spearmanr(d.DES, d[pur]); res[f'DES_vs_{pur}'] = dict(rho=r, p=pv, n=len(d))
for g in G[1:]:
    d = X[['DES', g]].dropna(); r, pv = stats.spearmanr(d.DES, d[g]); res[f'DES_vs_{g}_protein'] = dict(rho=r, p=pv, n=len(d))
# purity vs markers (reference)
for g in ['SFTPC', 'AGER', 'MYH11', 'EPCAM']:
    d = X[[g, 'TUMOR_PURITY_BYESTIMATE_RNASEQ']].dropna(); r, pv = stats.spearmanr(d[g], d['TUMOR_PURITY_BYESTIMATE_RNASEQ']); res[f'{g}_vs_purity'] = dict(rho=r, p=pv)
# partial: DES vs MYH11 given purity
def rpart(a, b, cz):
    ra, rb, rc = [stats.rankdata(v) for v in (a, b, cz)]
    f = lambda u, w: u - np.polyval(np.polyfit(w, u, 1), w); return np.corrcoef(f(ra, rc), f(rb, rc))[0, 1]
d = X[['DES', 'MYH11', 'SFTPC', 'TUMOR_PURITY_BYESTIMATE_RNASEQ']].dropna()
res['partial_DES_MYH11_given_purity'] = rpart(d.DES.values, d.MYH11.values, d.TUMOR_PURITY_BYESTIMATE_RNASEQ.values)
res['partial_DES_SFTPC_given_purity'] = rpart(d.DES.values, d.SFTPC.values, d.TUMOR_PURITY_BYESTIMATE_RNASEQ.values)
res['partial_DES_purity_given_MYH11'] = rpart(d.DES.values, d.TUMOR_PURITY_BYESTIMATE_RNASEQ.values, d.MYH11.values)
# ---------- CCLE: DES mRNA in lung lines
cl = pd.read_csv(f'{D}/ccle_broad_2019_data_clinical_sample.txt', sep='\t', comment='#')
cl.columns = [x.strip() for x in cl.columns]
hdr = pd.read_csv(f'{D}/ccle_broad_2019_data_mrna_seq_rpkm.txt', sep='\t', nrows=0).columns.tolist()
rows = {}
want = ['DES', 'MYH11', 'ACTA2', 'SFTPC', 'AGER', 'EPCAM', 'VIM', 'GAPDH', 'KRT19']
for ch in pd.read_csv(f'{D}/ccle_broad_2019_data_mrna_seq_rpkm.txt', sep='\t', chunksize=2000, index_col=0):
    sel = ch[ch.index.isin(want)]
    for g, r in sel.iterrows(): rows[g] = r
E_ = pd.DataFrame(rows)
E_['line'] = E_.index; E_['tissue'] = [x.split('_', 1)[1] if '_' in x else '' for x in E_.index]
res['ccle_n_lines'] = len(E_)
lung = E_[E_.tissue == 'LUNG']; res['ccle_lung_n'] = len(lung)
allDES = E_.DES.astype(float)
res['ccle_DES_lung_rpkm'] = dict(median=float(lung.DES.median()), max=float(lung.DES.max()), n_gt1=int((lung.DES > 1).sum()), n_gt0p1=int((lung.DES > 0.1).sum()), q75=float(lung.DES.quantile(.75)))
res['ccle_DES_all_median'] = float(allDES.median())
# reference: smooth muscle lines/tissues with high DES
top = E_.sort_values('DES', ascending=False).head(8)[['DES']]
res['ccle_top_DES_lines'] = {k: float(v) for k, v in top.DES.items()}
res['ccle_lung_SFTPC_median'] = float(lung.SFTPC.median()); res['ccle_lung_EPCAM_median'] = float(lung.EPCAM.median()); res['ccle_lung_VIM_median'] = float(lung.VIM.median())
res['ccle_lung_top_DES'] = {k: float(v) for k, v in lung.sort_values('DES', ascending=False).head(5).DES.items()}
json.dump(res, open(f'{OUT}/cptac_ccle_summary.json', 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
