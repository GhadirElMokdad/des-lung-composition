"""DES by cell type in the normal-lung atlas of the Krasnow lab (10x, CELLxGENE Census dataset 8c42cfd0-0b0a-46d5-910c-fc833d83c45e, 65,662 cells).
Output: results/des_krasnow_celltype.csv, results/des_krasnow_share.csv."""
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
from common import *
a = ad.read_h5ad(f'{CENSUS}/krasnow10x.h5ad'); names = a.raw.var['feature_name'].values
X = a.raw.X.tocsc() if sp.issparse(a.raw.X) else a.raw.X; tot = np.asarray(X.sum(1)).ravel()
c = np.asarray(X[:, np.where(names == 'DES')[0][0]].todense()).ravel()
d = pd.DataFrame({'ct': a.obs.cell_type.values, 'pos': c > 0, 'cpm': c / tot * 1e6})
r = d.groupby('ct', observed=True).agg(n=('pos', 'size'), pct_DES=('pos', 'mean'), cpm_DES=('cpm', 'mean')); r['pct_DES'] *= 100; r = r[r.n >= 100].sort_values('pct_DES', ascending=False)
share = (d.groupby('ct', observed=True).cpm.sum() / d.cpm.sum() * 100).sort_values(ascending=False).round(1).head(8)
r.to_csv(f'{OUT}/des_krasnow_celltype.csv'); share.to_csv(f'{OUT}/des_krasnow_share.csv'); print(r.head(10)); print(share)
