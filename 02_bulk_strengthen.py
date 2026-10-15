"""TCGA-LUAD (GDC hub) checks for the DES article: ROC comparators, adjusted Cox models, composition correlations.
Input: luad_genes.pkl produced by 01_tcga_luad_bulk.py (log2(TPM+1) of selected genes), TCGA-LUAD survival/clinical."""
import numpy as np, pandas as pd, json, sys
from scipy import stats
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import survdiff
from common import *
D = sys.argv[1] if len(sys.argv) > 1 else DERIVED
df = pd.read_pickle(f'{DERIVED}/luad_genes.pkl'); code = df.index.str[13:15]
X = df[(code == '01') | (code == '11')]; y = (X.index.str[13:15] == '01').astype(int)
def delong(y, preds):
    pos = y == 1; m, n = pos.sum(), (~pos).sum(); k = len(preds); V10 = np.zeros((k, m)); V01 = np.zeros((k, n)); aucs = np.zeros(k)
    for i, s in enumerate(preds):
        s = np.asarray(s); sp_, sn = s[pos], s[~pos]; r = stats.rankdata(np.concatenate([sp_, sn])); rp = stats.rankdata(sp_); rn = stats.rankdata(sn)
        aucs[i] = (r[:m].sum() - m * (m + 1) / 2) / (m * n); V10[i] = (r[:m] - rp) / n; V01[i] = 1 - (r[m:] - rn) / m
    return aucs, np.cov(V10) / m + np.cov(V01) / n
rows = []
genes = ['DES', 'MYH11', 'CNN1', 'TAGLN', 'ACTA2', 'CDH5', 'PECAM1', 'AGER', 'SFTPC', 'NAPSA', 'CLDN18']
for g in genes:
    sc = -X[g].values                       # score oriented so that higher = tumour (all are lower in tumour)
    au, S = delong(y, [sc, -X['DES'].values]); d = au[0] - au[1]
    p = 1.0 if g == 'DES' else 2 * stats.norm.sf(abs(d) / np.sqrt(S[0, 0] + S[1, 1] - 2 * S[0, 1]))
    rng = np.random.default_rng(1); bs = []
    pos = np.where(y == 1)[0]; neg = np.where(y == 0)[0]
    for _ in range(1000):
        i = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        bs.append(delong(y[i], [sc[i]])[0][0])
    rows.append(dict(gene=g, auc=au[0], ci_lo=np.percentile(bs, 2.5), ci_hi=np.percentile(bs, 97.5), delong_p_vs_DES=p))
pd.DataFrame(rows).to_csv(f'{OUT}/luad_roc_comparators.csv', index=False)
print(pd.DataFrame(rows).round(4).to_string())
# --- survival
sv = pd.read_csv(f'{GDC}/TCGA-LUAD.survival.tsv.gz', sep='\t', index_col=0)
cl = pd.read_csv(f'{GDC}/TCGA-LUAD.clinical.tsv.gz', sep='\t', index_col=0)
T = df[code == '01'].join(sv, how='inner'); T = T[~T.index.str[:12].duplicated()]
T['pid'] = T.index.str[:12]
cl2 = cl.copy(); cl2['pid'] = cl2.index.str[:12]; cl2 = cl2[~cl2.pid.duplicated()].set_index('pid')
T['age'] = pd.to_numeric(cl2.loc[T.pid, 'age_at_index.demographic']).values
st = cl2.loc[T.pid, 'ajcc_pathologic_stage.diagnoses'].astype(str).str.extract(r'Stage (I{1,3}|IV)')[0].map({'I': 1, 'II': 2, 'III': 3, 'IV': 4}).values
T['stage'] = st
T['male'] = (cl2.loc[T.pid, 'gender.demographic'] == 'male').astype(int).values
print('survival cohort', len(T), 'events', int(T.OS.sum()), 'missing age', T.age.isna().sum(), 'missing stage', int(np.isnan(T.stage).sum()))
res = []
def cox(cols, label, data):
    d = data.dropna(subset=cols + ['OS', 'OS.time']); d = d[d['OS.time'] > 0]
    Z = d[cols].copy()
    for c in cols:
        if c not in ('male',): Z[c] = (Z[c] - Z[c].mean()) / Z[c].std() if c in genes else Z[c]
    m = PHReg(d['OS.time'].values, Z.values, status=d.OS.values).fit()
    ci = np.exp(m.conf_int())
    for i, c in enumerate(cols):
        res.append(dict(model=label, n=len(d), events=int(d.OS.sum()), term=c, HR=np.exp(m.params[i]), lo=ci[i, 0], hi=ci[i, 1], p=m.pvalues[i]))
cox(['DES'], 'DES alone (per SD)', T)
cox(['DES', 'age', 'stage', 'male'], 'DES + age + stage + sex', T)
cox(['DES', 'age', 'stage', 'male', 'SFTPC'], '+ SFTPC (normal-lung content)', T)
cox(['DES', 'age', 'stage', 'male', 'AGER'], '+ AGER', T)
cox(['DES', 'age', 'stage', 'male', 'MYH11'], '+ MYH11', T)
R = pd.DataFrame(res); R.to_csv(f'{OUT}/luad_cox_models.csv', index=False); print(R.round(4).to_string())
# cut-offs
cut = {}
for nm, q in [('top quartile', .75), ('median', .5), ('top tertile', 2 / 3)]:
    hi = (T.DES >= T.DES.quantile(q)).astype(int); s, p = survdiff(T['OS.time'], T.OS, hi); cut[nm] = dict(n_high=int(hi.sum()), p=float(p))
print(cut)
# --- fold changes and composition correlations
tum = df[code == '01']; nor = df[code == '11']; tpm = lambda x: 2 ** x - 1
fc = pd.DataFrame({'median_TPM_normal': tpm(nor).median(), 'median_TPM_tumour': tpm(tum).median()}); fc['fold_drop'] = fc.median_TPM_normal / fc.median_TPM_tumour
fc.to_csv(f'{OUT}/luad_fold_changes.csv')
cor = {g: stats.spearmanr(tum.DES, tum[g])[0] for g in ['MYH11', 'CNN1', 'TAGLN', 'ACTA2', 'PDGFRB', 'RGS5', 'CDH5', 'PECAM1', 'COL1A1', 'AGER', 'SFTPC', 'NAPSA', 'NKX2-1', 'EPCAM', 'KRT7', 'PTPRC']}
pd.Series(cor).to_csv(f'{OUT}/luad_des_correlations_tumour.csv')
# partial correlation of DES with smooth muscle score controlling for normal-lung epithelial score
z = lambda s: (s - s.mean()) / s.std()
sm = z(tum.MYH11) + z(tum.CNN1); lung = z(tum.SFTPC) + z(tum.AGER)
def rpartial(a, b, c):
    ra, rb, rc = [stats.rankdata(v) for v in (a, b, c)]
    def resid(u, w): bb = np.polyfit(w, u, 1); return u - np.polyval(bb, w)
    return np.corrcoef(resid(ra, rc), resid(rb, rc))[0, 1]
pc = dict(DES_vs_smoothmuscle=stats.spearmanr(tum.DES, sm)[0], DES_vs_lungscore=stats.spearmanr(tum.DES, lung)[0], sm_vs_lung=stats.spearmanr(sm, lung)[0],
          DES_sm_given_lung=rpartial(tum.DES.values, sm.values, lung.values), DES_lung_given_sm=rpartial(tum.DES.values, lung.values, sm.values))
print(pc); json.dump(dict(cutoffs=cut, partial=pc, cor=cor), open(f'{OUT}/luad_summary.json', 'w'), indent=1)
