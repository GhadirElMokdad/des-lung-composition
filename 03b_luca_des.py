"""DES by cell type in the LuCA extended atlas (CELLxGENE Census 2025-11-04, dataset 1e6a6ef9): LUAD tumour vs normal lung, donor as unit."""
import pandas as pd, numpy as np, json, sys, warnings
from scipy import stats
warnings.filterwarnings('ignore')
from common import *
D = sys.argv[1] if len(sys.argv) > 1 else CENSUS
o = pd.read_pickle(f'{D}/luca_obs.pkl'); v = pd.read_pickle(f'{D}/luca_var.pkl'); c = pd.read_pickle(f'{D}/luca_counts.pkl')
for k in ['disease', 'cell_type', 'tissue_general', 'donor_id']: o[k] = o[k].astype(str)
o = o[o.disease.isin(['lung adenocarcinoma', 'normal']) & (o.tissue_general == 'lung')].copy()
o['group'] = o.disease.map({'lung adenocarcinoma': 'tumour', 'normal': 'normal'})
def comp(ct):
    if ct == 'malignant cell': return 'malignant'
    if ct in ('smooth muscle cell',): return 'smooth muscle'
    if ct == 'pericyte': return 'pericyte'
    if 'fibroblast' in ct or ct == 'stromal cell': return 'fibroblast'
    if 'endothelial' in ct: return 'endothelial'
    if ct in ('pulmonary alveolar type 2 cell', 'pulmonary alveolar type 1 cell', 'club cell', 'multiciliated epithelial cell', 'epithelial cell of lung'): return 'epithelium'
    if ct == 'mesothelial cell': return 'mesothelial'
    return 'immune'
o['comp'] = o.cell_type.map(comp)
fname = dict(zip(v.soma_joinid, v.feature_name)); c['gene'] = c['var'].map(fname)
o = o.set_index('soma_joinid'); c = c[c.obs.isin(o.index)]
GENES = ['DES', 'AGER', 'SFTPC', 'CA4', 'ADH1B', 'BTNL9', 'SPP1', 'FAM83A', 'MYH11', 'ACTA2', 'PTPRC']
cnt = {g: c[c.gene == g].set_index('obs').v.reindex(o.index).fillna(0).values for g in GENES}
for g in GENES: o['c_' + g] = cnt[g]
tot = o.raw_sum.values
# ---- 1. DES by cell type, tumour and normal
def tab(g, grp, by='cell_type'):
    d = o[o.group == grp]; r = d.groupby(by).agg(n=('c_' + g, 'size'), pct=('c_' + g, lambda x: 100 * (x > 0).mean()), counts=('c_' + g, 'sum'), lib=('raw_sum', 'sum'))
    r['cpm'] = r.counts / r.lib * 1e6; r['share_of_counts'] = 100 * r.counts / r.counts.sum(); r['donors'] = d.groupby(by).donor_id.nunique(); return r
rows = []
for grp in ['normal', 'tumour']:
    for by in ['cell_type', 'comp']:
        r = tab('DES', grp, by).reset_index().rename(columns={by: 'label'}); r['level'] = by; r['group'] = grp; rows.append(r)
pd.concat(rows).to_csv(f'{OUT}/luca_des_by_celltype.csv', index=False)
for grp in ['normal', 'tumour']:
    print('==', grp); print(tab('DES', grp, 'comp').round(2).sort_values('share_of_counts', ascending=False).to_string())
    print(tab('DES', grp).round(2).query('n>=300').sort_values('pct', ascending=False).head(10).to_string())
# ---- 2. donor-level: stromal fraction and pseudobulk DES
dn = o.groupby(['group', 'donor_id']).agg(n=('c_DES', 'size'), des=('c_DES', 'sum'), lib=('raw_sum', 'sum'))
dn['cpm'] = dn.des / dn.lib * 1e6
fr = o.groupby(['group', 'donor_id', 'comp']).size().unstack(fill_value=0)
dn['n_nonimmune'] = fr.drop(columns='immune').sum(1)
dn['stromal_frac_all'] = (fr['smooth muscle'] + fr['pericyte'] + fr['fibroblast']) / fr.sum(1) * 100
dn['stromal_frac_nonimmune'] = (fr['smooth muscle'] + fr['pericyte'] + fr['fibroblast']) / dn.n_nonimmune * 100
dn['malignant_frac'] = fr['malignant'] / fr.sum(1) * 100
dn = dn[dn.n >= 500]; dn.to_csv(f'{OUT}/luca_donor_level.csv')
res = {}
for col in ['cpm', 'stromal_frac_all', 'stromal_frac_nonimmune']:
    a = dn.loc['tumour', col].dropna(); b = dn.loc['normal', col].dropna()
    res[col] = dict(n_tumour=len(a), n_normal=len(b), med_tumour=float(a.median()), med_normal=float(b.median()), p=float(stats.mannwhitneyu(a, b).pvalue))
# DES % positive within stromal cell types per donor (>=30 cells)
for comp_ in ['smooth muscle', 'fibroblast', 'pericyte', 'endothelial', 'epithelium', 'malignant']:
    d = o[o.comp == comp_].groupby(['group', 'donor_id']).agg(n=('c_DES', 'size'), pct=('c_DES', lambda x: 100 * (x > 0).mean()))
    d = d[d.n >= 30]
    try:
        a = d.loc['tumour', 'pct']; b = d.loc['normal', 'pct']; res['pct_' + comp_] = dict(n_tumour=len(a), n_normal=len(b), med_tumour=float(a.median()), med_normal=float(b.median()), p=float(stats.mannwhitneyu(a, b).pvalue))
    except KeyError:
        try: a = d.loc['tumour', 'pct']; res['pct_' + comp_] = dict(n_tumour=len(a), med_tumour=float(a.median()), n_normal=0)
        except KeyError: pass
# composition decomposition: pseudobulk DES expected from normal per-cell-type rates and tumour composition
nm = o[o.group == 'normal']; tm = o[o.group == 'tumour']
rate_n = nm.groupby('comp').apply(lambda d: d.c_DES.sum() / d.raw_sum.sum())          # DES counts per count in each compartment, normal lung
w_t = tm.groupby('comp').raw_sum.sum() / tm.raw_sum.sum(); w_n = nm.groupby('comp').raw_sum.sum() / nm.raw_sum.sum()
exp_t = (rate_n * w_t).sum() * 1e6; obs_n = (rate_n * w_n).sum() * 1e6; obs_t = tm.c_DES.sum() / tm.raw_sum.sum() * 1e6
res['decomposition_cpm'] = dict(normal_observed=float(obs_n), tumour_observed=float(obs_t), tumour_expected_from_normal_rates_and_tumour_composition=float(exp_t))
print(json.dumps(res, indent=1))
# ---- 3. cell-type origin of the eight retained genes (pooled over cells, by compartment)
rows = []
for grp in ['normal', 'tumour']:
    for g in ['DES', 'AGER', 'SFTPC', 'CA4', 'ADH1B', 'BTNL9', 'SPP1', 'FAM83A']:
        r = tab(g, grp, 'cell_type'); r = r[r.n >= 300]
        top = r.sort_values('share_of_counts', ascending=False).head(2)
        rows.append(dict(group=grp, gene=g, top1=top.index[0], top1_share=top.share_of_counts.iloc[0], top1_pct_pos=top.pct.iloc[0], top2=top.index[1], top2_share=top.share_of_counts.iloc[1],
                         max_pct_celltype=r.pct.idxmax(), max_pct=r.pct.max()))
g8 = pd.DataFrame(rows); g8.to_csv(f'{OUT}/luca_eight_genes_origin.csv', index=False); print(g8.round(1).to_string())
json.dump(res, open(f'{OUT}/luca_summary.json', 'w'), indent=1)
