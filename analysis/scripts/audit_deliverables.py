"""Audit retained vendor files; no read processing or taxonomic classification."""
from pathlib import Path
import csv, gzip, hashlib, json, re, statistics, zipfile
from collections import Counter
from datetime import datetime, timezone
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'vendor_analysis/extracted_results/OPL_mRNA_3_24201_Analysis_Results'
OUT=ROOT/'analysis/results';OUT.mkdir(parents=True,exist_ok=True)
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def workbook(path):
 with zipfile.ZipFile(path) as z:
  strings=[]
  if 'xl/sharedStrings.xml' in z.namelist():
   strings=[''.join(t.text or '' for t in x.findall('.//s:t',NS)) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS)]
  rows=[]
  for row in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('s:sheetData/s:row',NS):
   vals={}
   for c in row.findall('s:c',NS):
    col=0
    for ch in re.match('[A-Z]+',c.get('r')).group():col=col*26+ord(ch)-64
    value=c.find('s:v',NS);text=value.text if value is not None else ''.join(t.text or '' for t in c.findall('.//s:t',NS))
    vals[col-1]=strings[int(text)] if c.get('t')=='s' and text else text
   rows.append([vals.get(i,'') for i in range(max(vals,default=-1)+1)])
  return rows

def fasta(path):
 name=None;parts=[];header=''
 with path.open() as f:
  for line in f:
   if line.startswith('>'):
    if name is not None:yield name,header,''.join(parts).upper()
    header=line[1:].strip();name=header.split()[0];parts=[]
   else:parts.append(line.strip())
 if name is not None:yield name,header,''.join(parts).upper()

def table(name,header,rows):
 with (OUT/name).open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)

nt={};lengths=[];bases=Counter();digests=Counter();duplicate_ids=0
fa=SRC/'Transcriptome DeNovo Assembly/Dinoponera gigantea Assembly.fa'
for name,header,seq in fasta(fa):
 duplicate_ids+=name in nt;nt[name]=len(seq);lengths.append(len(seq));bases.update(seq);digests[hashlib.sha256(seq.encode()).digest()]+=1
ordered=sorted(lengths,reverse=True);total=sum(lengths);cumulative=0;n50=0
for length in ordered:
 cumulative+=length
 if cumulative>=total/2:n50=length;break
pep={};orf=Counter()
for name,header,seq in fasta(SRC/'Transcriptome DeNovo Assembly/Dinoponera gigantea Assembly.pep'):
 pep[name]=seq
 m=re.search('ORF type:(\\S+)',header)
 orf[m.group(1) if m else 'unspecified']+=1
assembly={'records':len(lengths),'unique_ids':len(nt),'duplicate_id_records':duplicate_ids,'total_bases':total,'minimum_length':min(lengths),'median_length':statistics.median(lengths),'maximum_length':max(lengths),'N50':n50,'GC_percent':100*(bases['G']+bases['C'])/total,'N_bases':bases['N'],'exact_duplicate_sequence_records_beyond_first':sum(n-1 for n in digests.values()),'peptide_records':len(pep),'peptide_ids_without_transcript':sum(re.sub(r'\.p\d+$','',x) not in nt for x in pep),'orf_types':dict(orf)}
counts={};species=[];samples=[];invalid=[]
for sample in ['AZ1','AZ2','AZ3']:
 rows=workbook(SRC/f'Expressed Genes/{sample}_Expressed_Genes.xlsx');data=rows[1:];seen=Counter();valid={};taxa=Counter();bad=0
 for row in data:
  row+=['']*max(0,4-len(row));tid=row[0].strip()
  if not tid:continue
  seen[tid]+=1
  try:
   val=float(row[3]);assert val>=0 and val.is_integer()
   valid[tid]=int(val)
  except (ValueError,AssertionError):bad+=1;invalid.append([sample,tid,row[3]])
  m=re.search(r'\[([^\[\]]+)\]\s*$',row[2]);taxa[m.group(1) if m else 'unassigned description']+=1
 counts[sample]=valid
 samples.append([sample,sum(seen.values()),len(seen),sum(n-1 for n in seen.values()),sum(x in nt for x in seen),sum(x not in nt for x in seen),bad,sum(valid.values())])
 species.extend([sample,t,n,'Description label only; not organism assignment or contamination evidence'] for t,n in taxa.most_common())
table('expression_audit.tsv',['sample','populated_ID_rows','unique_IDs','duplicate_ID_rows','IDs_in_transcript_FASTA','IDs_not_in_FASTA','invalid_count_rows','sum_delivered_counts'],samples)
table('annotation_species_labels.tsv',['sample','protein_description_species_label','rows','interpretation'],species)
table('invalid_expression_counts.tsv',['sample','ID','value'],invalid)
# Missing values remain NA; no zero imputation and no conversion to TPM.
ids=sorted(set().union(*(set(x) for x in counts.values())))
with (OUT/'delivered_counts_union.tsv.gz').open('wb') as raw:
 with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as compressed:
  compressed.write(('transcript_ID\tAZ1\tAZ2\tAZ3\n'+''.join(x+'\t'+'\t'.join(str(counts[s].get(x,'NA')) for s in counts)+'\n' for x in ids)).encode())
shared=set.intersection(*(set(x) for x in counts.values()))
amp=workbook(SRC/'AMP Analysis/AMP_Results.xlsx');headers=amp[0];out=[];queries=set()
for row in amp[1:]:
 row+=['']*max(0,len(headers)-len(row));d=dict(zip(headers,row));qid=d.get('Query ID','').strip()
 if not qid:continue
 queries.add(qid);fragment=d.get('Query seq','').replace(' ','').upper();full=pep.get(qid,'')
 out.append([qid,d.get('DRAMP AMP ID',''),d.get('Percentage Identity',''),len(fragment),len(full),qid in pep, bool(fragment) and fragment in full,'Short similarity fragments do not establish AMP activity'])
table('amp_sequence_audit.tsv',['query_ID','DRAMP_ID','reported_percent_identity','query_fragment_length_aa','full_predicted_peptide_length_aa','query_ID_in_peptide_FASTA','fragment_is_exact_substring','interpretation'],out)
with (OUT/'amp_candidate_full_sequences.faa').open('w') as f:
 for q in sorted(queries):
  if q in pep:f.write('>'+q+'\n'+pep[q]+'\n')
summary={'audit_date_utc':datetime.now(timezone.utc).isoformat(),'assembly':assembly,'expression':{'samples':samples,'union_IDs':len(ids),'shared_IDs_all_three':len(shared),'missing_cell_policy':'NA means absent from delivered workbook; not assumed zero','normalization':'None; original delivered counts preserved','design':'Same group; no treatment comparison; biological independence unconfirmed'},'AMP':{'match_rows':len(out),'distinct_queries':len(queries),'queries_found_in_peptide_FASTA':len(queries&set(pep)),'match_rows_query_fragment_under_20_aa':sum(x[3]<20 for x in out),'minimum_query_fragment_length_aa':min(x[3] for x in out)},'contamination_status':'Not screened. Species names in homology descriptions are not taxonomic classifications. No search database or classification tool available in the current environment.'}
(OUT/'deliverable_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
