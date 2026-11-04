"""Local BLOSUM62 re-alignment of full candidates against their original DRAMP targets."""
from pathlib import Path
import sys,csv,json,hashlib
from Bio import SeqIO,Align
from Bio.Align import substitution_matrices
root=Path(__file__).resolve().parents[2];out=root/'analysis/results';ref=Path(sys.argv[1])
queries={x.id:str(x.seq) for x in SeqIO.parse(out/'amp_candidate_full_sequences.faa','fasta')};targets={entry.splitlines()[0].split()[0]:''.join(''.join(entry.splitlines()[1:]).split()).upper() for entry in ref.read_text().split('>')[1:]}
a=Align.PairwiseAligner();a.mode='local';a.substitution_matrix=substitution_matrices.load('BLOSUM62');a.open_gap_score=-11;a.extend_gap_score=-0.5
rows=[]
with (out/'amp_sequence_audit.tsv').open() as f:
 for hit in csv.DictReader(f,delimiter='\t'):
  q,t=hit['query_ID'],hit['DRAMP_ID'];qs=queries[q].rstrip('*');ts=targets.get(t,'').rstrip('*')
  if not ts:rows.append([q,t,'target_absent_from_download']+['NA']*9);continue
  if set(qs+ts)-set(a.substitution_matrix.alphabet):rows.append([q,t,'unsupported_residue']+['NA']*9);continue
  alignments=a.align(qs,ts)
  if not alignments:rows.append([q,t,'no_positive_alignment']+['NA']*9);continue
  aln=alignments[0];qc,tc=aln.coordinates;identical=0;columns=0
  for i in range(len(qc)-1):
   n,m=qc[i+1]-qc[i],tc[i+1]-tc[i];columns+=max(n,m)
   if n and m:identical+=sum(x==y for x,y in zip(qs[qc[i]:qc[i+1]],ts[tc[i]:tc[i+1]]))
  rows.append([q,t,'aligned',len(qs),len(ts),round(aln.score,2),columns,round(100*identical/columns,2),round(100*(qc[-1]-qc[0])/len(qs),2),round(100*(tc[-1]-tc[0])/len(ts),2),int(qc[0])+1,int(qc[-1])])
with (out/'amp_full_sequence_realignments.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['query_ID','DRAMP_ID','status','query_length','target_length','local_score','alignment_columns','identity_percent','query_coverage_percent','target_coverage_percent','query_start_1based','query_end_1based']);w.writerows(rows)
(out/'amp_realignment_methods.json').write_text(json.dumps({'method':'Biopython local pairwise alignment, BLOSUM62, gap open -11, gap extend -0.5; first optimal alignment retained','scope':'Original reported query-target pairs only; not all-vs-database search','reference_url':'https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.fasta','reference_sha256':hashlib.file_digest(ref.open('rb'),'sha256').hexdigest(),'reference_records':len(targets),'pairs':len(rows),'status_counts':{s:sum(r[2]==s for r in rows) for s in set(r[2] for r in rows)},'limits':'No E-values or statistical significance calculated; similarity does not validate antimicrobial activity; full precursor coverage differs from mature peptide coverage'},indent=2)+'\n')
print('Alignment statuses:',{s:sum(r[2]==s for r in rows) for s in set(r[2] for r in rows)})
