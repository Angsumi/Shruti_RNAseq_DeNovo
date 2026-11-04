"""Join retained AMP matches, peptide headers and delivered counts without activity prediction."""
import runpy,re,csv
from pathlib import Path
from collections import defaultdict
r=Path(__file__).resolve().parents[2]
a=runpy.run_path(str(r/'analysis/scripts/audit_deliverables.py'))
matches=defaultdict(list)
for row in a['out']:matches[row[0]].append(row)
headers={n:h for n,h,s in a['fasta'](a['SRC']/'Transcriptome DeNovo Assembly/Dinoponera gigantea Assembly.pep')}
rows=[]
for q in sorted(a['queries']):
 seq=a['pep'].get(q,'');tid=re.sub(r'\.p\d+$','',q);m=re.search(r'ORF type:(\S+)',headers.get(q,''));orf=m.group(1) if m else 'unknown'
 hits=matches[q];longest=max(x[3] for x in hits);observed=[a['counts'][s].get(tid,'NA') for s in ['AZ1','AZ2','AZ3']]
 flags=[]
 if orf!='complete':flags.append('partial_ORF')
 if longest<20:flags.append('all_match_fragments_under_20aa')
 if any(x=='NA' for x in observed):flags.append('missing_delivered_count')
 rows.append([q,tid,len(seq),orf,len(hits),longest,round(longest/len(seq),4) if seq else 'NA',*observed,';'.join(flags) or 'review_full_alignment','Not validated; longest fragment/full precursor is not mature-peptide coverage'])
a['table']('amp_candidate_evidence.tsv',['peptide_ID','transcript_ID','full_length_aa','ORF_type','match_rows','longest_query_fragment_aa','longest_fragment_to_full_peptide_ratio','AZ1_delivered_count','AZ2_delivered_count','AZ3_delivered_count','review_flags','interpretation'],rows)
print('CANDIDATES',len(rows),'COMPLETE_ORF',sum(x[3]=='complete' for x in rows),'ALL_FRAGMENTS_UNDER20',sum(x[5]<20 for x in rows),'COUNTS_ALL3',sum(all(y!='NA' for y in x[7:10]) for x in rows))
