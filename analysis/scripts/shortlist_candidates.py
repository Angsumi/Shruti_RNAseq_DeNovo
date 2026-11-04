"""Evidence-based review tiers, not an AMP activity classifier."""
from pathlib import Path
import csv,json
r=Path(__file__).resolve().parents[2];o=r/'analysis/results'
def read(name):return list(csv.DictReader((o/name).open(),delimiter='\t'))
e=read('amp_candidate_evidence.tsv');p={x['peptide_ID']:x for x in read('phobius_candidate_summary.tsv')};hits={}
for x in read('amp_full_sequence_realignments.tsv'):
 if x['status']=='aligned':hits.setdefault(x['query_ID'],[]).append(x)
rows=[]
for x in e:
 q=x['peptide_ID'];sp=p[q]['signal_peptide_predicted']=='True';tm=int(p[q]['transmembrane_helices']);complete=x['ORF_type']=='complete';observed=sum(x[s+'_delivered_count']!='NA' for s in ['AZ1','AZ2','AZ3'])
 # Select one actual alignment by reference coverage then score; do not combine maxima from different hits.
 h=max(hits.get(q,[]),key=lambda t:(float(t['target_coverage_percent']),float(t['local_score'])),default={})
 tier='A: complete ORF, signal peptide, no TM' if complete and sp and tm==0 else 'B: signal peptide, no TM; partial ORF' if sp and tm==0 else 'C: review other evidence'
 rows.append({**x,'signal_peptide_predicted':sp,'TM_helices':tm,'samples_with_delivered_counts':observed,'review_tier':tier,'selected_DRAMP_ID':h.get('DRAMP_ID','NA'),'selected_target_coverage_percent':h.get('target_coverage_percent','NA'),'selected_identity_percent':h.get('identity_percent','NA'),'selected_query_coverage_percent':h.get('query_coverage_percent','NA'),'claim_limit':'Review priority only; no functional validation or novelty claim'})
rows.sort(key=lambda x:(x['review_tier'],-x['samples_with_delivered_counts'],x['peptide_ID']))
for name,data in [('amp_integrated_review.tsv',rows),('amp_priority_review.tsv',[x for x in rows if x['review_tier'].startswith(('A:','B:'))])]:
 with (o/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
summary={t:sum(x['review_tier']==t for x in rows) for t in sorted(set(x['review_tier'] for x in rows))};(o/'amp_review_tiers.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
