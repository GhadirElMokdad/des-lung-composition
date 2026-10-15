"""Extract the genes of interest from the TCGA-LUAD STAR TPM matrix (GDC hub) -> data/derived/luad_genes.pkl"""
import pandas as pd,numpy as np,gzip,json
from common import *
from scipy import stats
pm=pd.read_csv(f'{GDC}/gencode.v36.annotation.gtf.gene.probemap',sep='\t')
ids=pm.set_index('id')['gene']
want=['DES','ACTA2','MYH11','CNN1','TAGLN','CDH5','PECAM1','AGER','SFTPC','NAPSA','NKX2-1','EPCAM','PTPRC','COL1A1','PDGFRB','RGS5','KRT7','CLDN18']
sel=[i for i,g in ids.items() if g in want]
rows={}
with gzip.open(f'{GDC}/TCGA-LUAD.star_tpm.tsv.gz','rt') as f:
    hdr=f.readline().rstrip('\n').split('\t')
    for line in f:
        k=line.split('\t',1)[0]
        if k in sel:
            p=line.rstrip('\n').split('\t'); rows[k]=np.array(p[1:],float)
df=pd.DataFrame({ids[k]:v for k,v in rows.items()},index=hdr[1:])  # log2(TPM+1)
df.to_pickle(f'{DERIVED}/luad_genes.pkl'); print(df.shape, df.columns.tolist())
