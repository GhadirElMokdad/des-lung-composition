"""Spatial test of the DES tissue-composition hypothesis on 10x Visium NSCLC sections (BioStudies E-MTAB-13530, Hu et al. Nat Commun 2024).
Input directory <D> holds, per section S (e.g. P10_T1, P10_B1, D1_1): S-filtered_feature_bc_matrix.h5 and S-spatial.tar.
Usage: python 08_spatial_visium.py <D> [sdrf]   (sdrf optional: histology per patient; otherwise inferred from NAPSA/NKX2-1 vs KRT5/TP63)"""
import sys, os, glob, tarfile, io, json, re
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy import stats
from common import *
D = sys.argv[1] if len(sys.argv) > 1 else VISIUM
META = json.load(open(sys.argv[2] if len(sys.argv) > 2 else f'{VISIUM}/samples.json'))
SM = ['MYH11', 'CNN1', 'TAGLN', 'ACTA2']; EN = ['CDH5', 'PECAM1', 'VWF', 'CLDN5']; EP = ['EPCAM', 'KRT7', 'KRT19', 'NAPSA', 'NKX2-1', 'SFTPC']; AD = ['NAPSA', 'NKX2-1']; SQ = ['KRT5', 'TP63', 'SOX2']
GEN = ['DES'] + SM + EN + EP + SQ + ['PTPRC']
def load(sec):
    with h5py.File(f'{D}/{sec}-filtered_feature_bc_matrix.h5', 'r') as f:
        m = f['matrix']; shape = tuple(m['shape'][:]); X = sp.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=shape)
        names = np.array([x.decode() for x in m['features']['name'][:]]); bc = np.array([x.decode() for x in m['barcodes'][:]])
    idx = {g: np.where(names == g)[0] for g in GEN}
    cnt = pd.DataFrame({g: np.asarray(X[idx[g][0], :].todense()).ravel() if len(idx[g]) else np.zeros(len(bc)) for g in GEN}, index=bc)
    tot = np.asarray(X.sum(0)).ravel(); cnt['total'] = tot
    with tarfile.open(f'{D}/{sec}-spatial.tar') as t:
        mem = [x for x in t.getmembers() if re.search(r'tissue_positions.*\.csv$', x.name)][0]
        pos = pd.read_csv(io.BytesIO(t.extractfile(mem).read()), header=None if 'list' in mem.name else 0)
    pos = pos.iloc[:, :6]; pos.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'px_row', 'px_col']; pos = pos.set_index('barcode')
    return cnt.join(pos[['in_tissue', 'array_row', 'array_col', 'px_row', 'px_col']], how='inner').query('in_tissue == 1')
def score(c, genes):
    lg = np.log1p(c[genes].div(c.total, axis=0) * 1e4); z = (lg - lg.mean()) / lg.std().replace(0, 1); return z.mean(1)
secs = sorted({os.path.basename(p).split('-filtered')[0] for p in glob.glob(f'{D}/*-filtered_feature_bc_matrix.h5')})
rows = []; keep = {}
for s in secs:
    try: c = load(s)
    except Exception as e: print('skip', s, e); continue
    if len(c) < 200: continue
    grp = 'tumour' if '_T' in s else 'adjacent' if '_B' in s else 'healthy'
    c['sm'] = score(c, SM); c['en'] = score(c, EN); c['ep'] = score(c, EP); c['ad'] = score(c, AD); c['sq'] = score(c, SQ)
    hist = {'lung adenocarcinoma': 'LUAD', 'lung squamous cell carcinoma': 'LUSC', 'non-small cell carcinoma': 'NSCLC-NOS', 'normal': 'healthy'}.get(META.get(s, ['', ''])[1], 'LUAD-like' if c.ad.mean() > c.sq.mean() else 'LUSC-like')
    des = c.DES.values; pos = des > 0
    smhi = c.sm >= c.sm.quantile(.90); enhi = (c.en >= c.en.quantile(.90)) & ~smhi
    epi = (c.ep >= c.ep.quantile(.60)) & (c.sm <= c.sm.median()) & (c.en <= c.en.median())     # epithelial/tumour-nest spots without smooth-muscle or vascular signal
    r = dict(section=s, patient=s.split('_')[0], group=grp, histology=hist, n_spots=len(c), des_counts=int(des.sum()), des_cpm=des.sum() / c.total.sum() * 1e6, pct_spots_des=100 * pos.mean(),
             pct_des_in_smhi=100 * (smhi & pos).sum() / max(pos.sum(), 1), pct_spots_des_smhi=100 * pos[smhi.values].mean() if smhi.sum() else np.nan, pct_spots_des_epi=100 * pos[epi.values].mean() if epi.sum() else np.nan,
             pct_spots_des_enhi=100 * pos[enhi.values].mean() if enhi.sum() else np.nan, share_counts_epi=100 * des[epi.values].sum() / max(des.sum(), 1), share_counts_smhi=100 * des[smhi.values].sum() / max(des.sum(), 1), share_counts_enhi=100 * des[enhi.values].sum() / max(des.sum(), 1),
             n_epi=int(epi.sum()), n_smhi=int(smhi.sum()))
    for nm in ['sm', 'en', 'ep']: r['rho_des_' + nm] = stats.spearmanr(c.DES, c[nm])[0] if pos.sum() > 5 else np.nan
    rows.append(r); keep[s] = c
R = pd.DataFrame(rows); R.to_csv(f'{OUT}/spatial_visium_sections.csv', index=False); print(R.round(2).to_string())
summ = {}
R_f = R[(R.histology == 'LUAD') | (R.group == 'healthy')] if (R.histology == 'LUAD').any() else R
for g, d in R_f.groupby('group'):
    summ[g] = dict(n_sections=len(d), median_des_cpm=float(d.des_cpm.median()), median_rho_sm=float(d.rho_des_sm.median()), median_rho_ep=float(d.rho_des_ep.median()), median_pct_des_in_smhi=float(d.pct_des_in_smhi.median()),
                   median_share_counts_epi=float(d.share_counts_epi.median()), median_pct_spots_des_smhi=float(d.pct_spots_des_smhi.median()), median_pct_spots_des_epi=float(d.pct_spots_des_epi.median()))
R_all = R; R = R[(R.histology == 'LUAD') | (R.group == 'healthy')] if (R.histology == 'LUAD').any() else R
t = R[R.group == 'tumour']; a = R[R.group == 'adjacent']
summ['LUSC_tumour_median_des_cpm'] = float(R_all[(R_all.histology == 'LUSC') & (R_all.group == 'tumour')].des_cpm.median()) if (R_all.histology == 'LUSC').any() else None
summ['LUSC_adjacent_median_des_cpm'] = float(R_all[(R_all.histology == 'LUSC') & (R_all.group == 'adjacent')].des_cpm.median()) if (R_all.histology == 'LUSC').any() else None
if len(t) > 2 and len(a) > 2: summ['tumour_vs_adjacent_des_cpm_MWU_p'] = float(stats.mannwhitneyu(t.des_cpm, a.des_cpm).pvalue)
pt = t.groupby('patient').des_cpm.median(); pa = a.groupby('patient').des_cpm.median(); both = pt.index.intersection(pa.index)
if len(both) > 3: summ['paired_patients'] = int(len(both)); summ['paired_wilcoxon_p'] = float(stats.wilcoxon(pt[both], pa[both]).pvalue); summ['paired_n_tumour_lower'] = int((pt[both] < pa[both]).sum())
tt = t.dropna(subset=['pct_spots_des_epi', 'pct_spots_des_smhi'])
if len(tt) > 2: summ['tumour_pct_spots_des_epi_vs_smhi_wilcoxon_p'] = float(stats.wilcoxon(tt.pct_spots_des_epi, tt.pct_spots_des_smhi).pvalue); summ['tumour_n_sections_epi_vs_smhi'] = int(len(tt)); summ['tumour_n_epi_lower'] = int((tt.pct_spots_des_epi < tt.pct_spots_des_smhi).sum())
summ['median_enrichment_des_pos_spots_in_top10pct_sm_tumour'] = float((t.pct_des_in_smhi / 10).median()); summ['median_enrichment_adjacent'] = float((a.pct_des_in_smhi / 10).median())
json.dump(summ, open(f'{OUT}/spatial_visium_summary.json', 'w'), indent=1); print(json.dumps(summ, indent=1))
# figure: one tumour and one adjacent section of the same patient
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
pick = [s for s in R[(R.group == 'tumour')].sort_values('des_counts', ascending=False).section][:1]
if pick:
    pat = pick[0].split('_')[0]; adj = [s for s in R[(R.group == 'adjacent') & (R.patient == pat)].section][:1]; show = pick + adj
    fig, ax = plt.subplots(len(show), 3, figsize=(10, 3.4 * len(show)), squeeze=False)
    for i, s in enumerate(show):
        c = keep[s]
        for j, (nm, lab, cm) in enumerate([('DES', 'DES counts', 'Reds'), ('sm', 'Smooth-muscle score', 'Blues'), ('ep', 'Epithelial/tumour score', 'Greens')]):
            v = c[nm].values.astype(float); a_ = ax[i, j]; sc = a_.scatter(c.px_col, -c.px_row, c=v, s=6, cmap=cm, vmax=np.percentile(v, 99) if nm != 'DES' else max(2, np.percentile(v, 99.5)), vmin=0 if nm == 'DES' else None); a_.set_aspect('equal'); a_.axis('off'); a_.set_title(f'{s}: {lab}', fontsize=8); fig.colorbar(sc, ax=a_, shrink=.6)
    fig.tight_layout(); fig.savefig(f'{OUT}/Figure7_spatial_DES.png', dpi=200)
