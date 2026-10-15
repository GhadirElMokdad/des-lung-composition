"""Shared paths. Set DES_DATA to move the data directory (default: ./data)."""
import os
DATA = os.environ.get('DES_DATA', 'data')
GDC = f'{DATA}/gdc'; CBIO = f'{DATA}/cbio'; CENSUS = f'{DATA}/census'; VISIUM = f'{DATA}/visium'; DERIVED = f'{DATA}/derived'
OUT = 'results'
for d in (OUT, DERIVED): os.makedirs(d, exist_ok=True)
