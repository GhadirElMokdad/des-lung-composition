"""Europe PMC replication of the targeted literature search (gene symbol + lung adenocarcinoma / LUAD in title or abstract)."""
import json, urllib.parse, urllib.request, datetime
from common import *
def hits(q):
    u = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode({'query': q, 'format': 'json', 'pageSize': 1, 'resultType': 'lite'})
    return json.load(urllib.request.urlopen(u, timeout=60))['hitCount']
out = {'date': datetime.date.today().isoformat(), 'queries': {}}
lu = '("lung adenocarcinoma" OR LUAD)'
for g in ['FAM83A', 'SPP1', 'BTNL9', 'CA4', 'AGER', 'ADH1B', 'SFTPC']:
    out['queries'][g] = dict(q=f'TITLE_ABS:"{g}" AND TITLE_ABS:{lu}', n=hits(f'TITLE_ABS:"{g}" AND TITLE_ABS:{lu}'))
for name, q in [('DES (symbol, title/abstract)', f'TITLE_ABS:"DES" AND TITLE_ABS:{lu}'), ('desmin', f'TITLE_ABS:"desmin" AND TITLE_ABS:{lu}'), ('desmin, title only + LUAD/NSCLC', 'TITLE:"desmin" AND (TITLE_ABS:"lung adenocarcinoma" OR TITLE_ABS:"non-small cell lung")'), ('desmin gene expression + lung cancer', 'TITLE_ABS:"desmin" AND TITLE_ABS:"lung cancer" AND (TITLE_ABS:expression OR TITLE_ABS:prognos* OR TITLE_ABS:biomarker)')]:
    out['queries'][name] = dict(q=q, n=hits(q))
json.dump(out, open('results/literature_replication.json', 'w'), indent=1); print(json.dumps(out, indent=1))
