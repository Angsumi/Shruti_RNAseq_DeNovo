"""Local UniVec BLAST screen using NCBI VecScreen scoring rules; no automatic trimming."""
from pathlib import Path
import csv,gzip,json,hashlib
from collections import defaultdict,Counter
r=Path(__file__).resolve().parents[2];out=r/'analysis/results';hits=defaultdict(list)
for line in (r/'.runtime/univec_hits.tsv').open():
 f=line.rstrip('\n').split('\t');assert len(f)==13;hits[f[0]].append(f)
counts=Counter();byquery=defaultdict(set);strength={'below_weak':0,'weak':1,'moderate':2,'strong':3}
with gzip.open(out/'univec_alignment_review.tsv.gz','wt',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['transcript_ID','UniVec_ID','identity_percent','alignment_length','transcript_length','start','end','raw_score','position_category','match_strength','UniVec_title'])
 for q,hs in hits.items():
  spans=[(min(int(h[5]),int(h[6])),max(int(h[5]),int(h[6])),int(h[9])) for h in hs];terminal=[a<=25 or int(hs[i][4])-b<=25 for i,(a,b,s) in enumerate(spans)]
  def category(score,t):
   return 'strong' if score>=(24 if t else 30) else 'moderate' if score>=(19 if t else 25) else 'weak' if score>=(16 if t else 23) else 'below_weak'
  # Treat proximity to a currently significant match as terminal; repeat until stable.
  changed=True
  while changed:
   changed=False;significant=[i for i,(_,_,s) in enumerate(spans) if category(s,terminal[i])!='below_weak']
   for i,(a,b,s) in enumerate(spans):
    if terminal[i]:continue
    if any(i!=j and max(0,a-spans[j][1],spans[j][0]-b)<=25 for j in significant):terminal[i]=True;changed=True
  for i,h in enumerate(hs):
   a,b,score=spans[i];cat=category(score,terminal[i]);counts[cat]+=1
   if cat!='below_weak':byquery[q].add(cat)
   w.writerow([q,h[1],h[2],h[3],h[4],a,b,score,'terminal_or_near_match' if terminal[i] else 'internal',cat,h[12]])
with (out/'univec_flagged_transcripts.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['transcript_ID','strongest_match','action'])
 for q in sorted(byquery):w.writerow([q,max(byquery[q],key=lambda x:strength[x]),'Review only; no sequence removed or trimmed'])
summary={'method':'BLASTn 2.17.0+ vs UniVec; VecScreen-style score categories. Local reimplementation, not official FCS-adaptor output.','parameters':'-task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -max_target_seqs 100 -num_threads 6','category_source':'https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/','reference_url':'https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec','reference_sha256':hashlib.file_digest((r/'.runtime/references/UniVec').open('rb'),'sha256').hexdigest(),'alignment_categories':dict(counts),'transcripts_with_weak_or_stronger_match':len(byquery),'transcripts_with_moderate_or_stronger_match':sum(bool(x&{'moderate','strong'}) for x in byquery.values()),'transcripts_with_strong_match':sum('strong' in x for x in byquery.values()),'limits':'UniVec includes biological sequences. Match strength alone is not proof of laboratory contamination. Whole-body microbes can have legitimate homologues. Review source entries, boundaries and library adapters before any removal; output capped at 100 reference targets per query.'}
(out/'univec_methods.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
