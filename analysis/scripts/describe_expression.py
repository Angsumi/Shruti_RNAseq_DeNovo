from pathlib import Path
import csv,gzip,statistics,json
from itertools import combinations
r=Path(__file__).resolve().parents[2];o=r/'analysis/results';samples=['AZ1','AZ2','AZ3'];data={s:{} for s in samples}
with gzip.open(o/'delivered_counts_union.tsv.gz','rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  for s in samples:
   if x[s]!='NA':data[s][x['transcript_ID']]=int(x[s])
with (o/'expression_descriptive_summary.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['sample','delivered_rows','observed_count_sum','median_delivered_count','minimum','maximum','zero_count_rows','interpretation'])
 for s in samples:
  vals=list(data[s].values());w.writerow([s,len(vals),sum(vals),statistics.median(vals),min(vals),max(vals),vals.count(0),'Delivered selected rows only; not complete library totals or TPM'])
with (o/'expression_membership_overlap.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['sample_1','sample_2','shared_delivered_IDs','union_delivered_IDs','Jaccard_membership','interpretation'])
 for a,b in combinations(samples,2):
  x,y=set(data[a]),set(data[b]);w.writerow([a,b,len(x&y),len(x|y),round(len(x&y)/len(x|y),6),'Workbook membership overlap; absence is not zero expression'])
