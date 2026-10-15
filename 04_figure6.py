import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
BL, OR, GR, AQ = '#2a78d6', '#eb6834', '#8a8f98', '#1baf7a'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .8})
ct = pd.read_csv('results/luca_des_by_celltype.csv'); ct = ct[ct.level == 'comp']
fc = pd.read_csv('results/luad_fold_changes.csv', index_col=0).fold_drop
g8 = pd.read_csv('results/luca_eight_genes_origin.csv')
fig = plt.figure(figsize=(7.4, 9.8)); gs = fig.add_gridspec(3, 1, hspace=.55, height_ratios=[1.1, 1, 1.25])
# A: % DES-positive cells by compartment, normal vs tumour
ax = fig.add_subplot(gs[0]); order = ['smooth muscle', 'fibroblast', 'pericyte', 'endothelial', 'malignant', 'epithelium', 'immune']
lab = {'smooth muscle': 'Smooth muscle', 'fibroblast': 'Fibroblast', 'pericyte': 'Pericyte', 'endothelial': 'Endothelial', 'malignant': 'Malignant', 'epithelium': 'Non-malignant epithelium', 'immune': 'Immune'}
y = np.arange(len(order)); h = .36
for k, (grp, col, off) in enumerate([('normal', BL, -h / 2), ('tumour', OR, h / 2)]):
    d = ct[ct.group == grp].set_index('label').reindex(order)
    ax.barh(y + off, d.pct, height=h - .04, color=col, label='Normal lung' if grp == 'normal' else 'LUAD tumour')
    for yy, v, n in zip(y + off, d.pct, d.n): ax.text(v + 1, yy, f'{v:.1f}%  (n = {int(n):,})' if v >= .1 else f'<0.1%  (n = {int(n):,})', va='center', fontsize=7.5, color='#333')
ax.set_yticks(y); ax.set_yticklabels([lab[o] for o in order]); ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_xlabel('Cells with at least one DES count (%)'); ax.legend(frameon=False, loc='lower right')
fig.text(.01, {'A':.985,'B':.655,'C':.335}['A'], 'A', fontsize=13, fontweight='bold')
# B: bulk fold drop
ax = fig.add_subplot(gs[1]); b = fc[['DES', 'MYH11', 'CNN1', 'TAGLN', 'ACTA2', 'CDH5', 'PECAM1', 'AGER', 'SFTPC', 'CLDN18']].sort_values(ascending=True)
names = {'DES': 'DES', 'MYH11': 'MYH11 (smooth muscle)', 'CNN1': 'CNN1 (smooth muscle)', 'TAGLN': 'TAGLN (smooth muscle)', 'ACTA2': 'ACTA2 (smooth muscle)', 'CDH5': 'CDH5 (endothelium)', 'PECAM1': 'PECAM1 (endothelium)', 'AGER': 'AGER (alveolar type 1)', 'SFTPC': 'SFTPC (alveolar type 2)', 'CLDN18': 'CLDN18 (alveolar epithelium)'}
cols = [OR if g == 'DES' else GR for g in b.index]
ax.barh(range(len(b)), b.values, color=cols, height=.62); ax.set_xscale('log'); ax.set_yticks(range(len(b))); ax.set_yticklabels([names[g] for g in b.index])
for i, v in enumerate(b.values): ax.text(v * 1.08, i, f'{v:.1f}×', va='center', fontsize=8)
ax.set_xlim(1, 600); ax.set_xlabel('Fall in median TPM, normal to tumour (TCGA-LUAD, 59 normal, 528 tumours)')
fig.text(.01, {'A':.985,'B':.655,'C':.335}['B'], 'B', fontsize=13, fontweight='bold')
# C: eight genes, share of counts from the top cell type in normal lung
ax = fig.add_subplot(gs[2]); d = g8[g8.group == 'normal'].set_index('gene').loc[['DES', 'AGER', 'SFTPC', 'CA4', 'ADH1B', 'BTNL9', 'SPP1', 'FAM83A']]
t = g8[g8.group == 'tumour'].set_index('gene').loc[d.index]
ypos = np.arange(len(d))
ax.barh(ypos - h / 2, d.top1_share, height=h - .04, color=BL, label='Normal lung'); ax.barh(ypos + h / 2, t.top1_share, height=h - .04, color=OR, label='LUAD tumour')
short = lambda s: s.replace('pulmonary alveolar type 1 cell', 'alveolar type 1').replace('pulmonary alveolar type 2 cell', 'alveolar type 2').replace('capillary endothelial cell', 'capillary endothelium').replace('fibroblast of lung', 'fibroblast').replace('smooth muscle cell', 'smooth muscle').replace('multiciliated epithelial cell', 'ciliated cell').replace('epithelial cell of lung', 'epithelial cell').replace('malignant cell', 'malignant')
for i, g in enumerate(d.index):
    ax.text(d.top1_share.iloc[i] + 1, i - h / 2, short(d.top1.iloc[i]), va='center', fontsize=7.5); ax.text(t.top1_share.iloc[i] + 1, i + h / 2, short(t.top1.iloc[i]), va='center', fontsize=7.5)
ax.set_yticks(ypos); ax.set_yticklabels(d.index, fontstyle='italic'); ax.invert_yaxis(); ax.set_xlim(0, 135); ax.set_xticks([0, 25, 50, 75, 100]); ax.set_xlabel('Share of the gene’s counts from its main cell type (%)'); ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(.5, 1.0), ncol=2)
fig.text(.01, {'A':.985,'B':.655,'C':.335}['C'], 'C', fontsize=13, fontweight='bold')
fig.subplots_adjust(left=.30, right=.97, top=.98, bottom=.06); fig.savefig('results/Figure6_DES_composition.png', dpi=220); print('ok')
