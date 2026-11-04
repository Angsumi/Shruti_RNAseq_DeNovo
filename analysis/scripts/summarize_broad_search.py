"""Summarize local Swiss-Prot searches; homology flags are not definitive contamination calls."""
from pathlib import Path
import csv,gzip,re,json,hashlib
from collections import Counter,defaultdict
ROOT=Path(__file__).resolve().parents[2];RUN=ROOT/'.runtime';OUT=ROOT/'analysis/results';REF=RUN/'references'
parents={}
with (REF/'nodes.dmp').open() as f:
 for line in f:
  a=line.split('|');parents[int(a[0])]=int(a[1])
groups=[(2,'Bacteria'),(2157,'Archaea'),(10239,'Viruses'),(4751,'Fungi'),(33090,'Viridiplantae'),(6656,'Arthropoda'),(33208,'Other_Metazoa')]
cache={}
def lineage(t):
 if t in cache:return cache[t]
 path=set();x=t
 while x not in path and x in parents:
  path.add(x)
  if x==1:break
  x=parents[x]
 group=next((name for tax,name in groups if tax in path),'Other_or_unresolved')
 value=(group,36668 in path,7399 in path);cache[t]=value;return value
reference={};reference_count=0
with gzip.open(REF/'uniprot_sprot.fasta.gz','rt') as f:
 for line in f:
  if line.startswith('>'):
   reference_count+=1;name=line[1:].split()[0];m=re.search(r' OX=(\d+)',line);tax=int(m.group(1)) if m else 0;reference[name]=(tax,*lineage(tax))
fields=['query_ID','review_status','best_reference_ID','taxon_ID','reference_group','best_hit_identity_percent','query_coverage_percent','reference_coverage_percent','E_value','bitscore','near_top_reference_groups','near_top_hit_count','best_reference_title']
summaries={}
def records(path):
 name=None
 with path.open() as f:
  for line in f:
   if line.startswith('>'):yield line[1:].split()[0]
def summarize(kind,raw,input_path):
 matches=defaultdict(list)
 for line in raw.open():
  f=line.rstrip('\n').split('\t');assert len(f)==13
  matches[f[0]].append(f)
 statuses=Counter();group_counts=Counter();best_rows={}
 with gzip.open(OUT/f'swissprot_{kind}_review.tsv.gz','wt',newline='') as handle:
  writer=csv.writer(handle,delimiter='\t',lineterminator='\n');writer.writerow(fields)
  for q in records(input_path):
   hs=matches.get(q,[])
   if not hs:row=[q,'no_reported_match']+['NA']*11
   else:
    best=max(hs,key=lambda f:float(f[11]));tax,group,ant,hym=reference.get(best[1],(0,'Other_or_unresolved',False,False));score=float(best[11]);near=[h for h in hs if float(h[11])>=score*.9];near_groups={reference.get(h[1],(0,'Other_or_unresolved',False,False))[1] for h in near}
    qcov=100*(abs(int(best[7])-int(best[6]))+1)/int(best[4]);scov=100*(abs(int(best[9])-int(best[8]))+1)/int(best[5]);strong=qcov>=50 and float(best[2])>=40 and float(best[10])<=1e-10
    status='bacterial_homology_review' if strong and near_groups=={'Bacteria'} else 'non_animal_homology_review' if strong and len(near_groups)==1 and group in ('Archaea','Viruses','Fungi','Viridiplantae') else 'arthropod_homology_support' if strong and near_groups=={'Arthropoda'} else 'mixed_or_limited_homology'
    row=[q,status,best[1],tax,group,best[2],round(qcov,2),round(scov,2),best[10],best[11],';'.join(sorted(near_groups)),len(near),best[12]];group_counts[group]+=1
   statuses[row[1]]+=1;writer.writerow(row);best_rows[q]=dict(zip(fields,row))
 summaries[kind]={'queries':sum(statuses.values()),'matched_queries':len(matches),'status_counts':dict(statuses),'best_reference_group_counts':dict(group_counts),'raw_alignment_rows':sum(len(x) for x in matches.values())}
 return best_rows
src=ROOT/'vendor_analysis/extracted_results/OPL_mRNA_3_24201_Analysis_Results/Transcriptome DeNovo Assembly'
protein=summarize('protein',RUN/'swissprot_protein_hits.tsv',src/'Dinoponera gigantea Assembly.pep')
transcript=summarize('transcript',RUN/'swissprot_transcript_hits.tsv',src/'Dinoponera gigantea Assembly.fa')
with (OUT/'amp_broad_annotation_review.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(fields)
 for x in csv.DictReader((OUT/'amp_integrated_review.tsv').open(),delimiter='\t'):w.writerow([protein[x['peptide_ID']][k] for k in fields])
summary={'software':'DIAMOND 2.2.8','reference':'Reviewed UniProtKB/Swiss-Prot','reference_release':(REF/'uniprot_reldate.txt').read_text(),'reference_records':reference_count,'reference_sha256':hashlib.file_digest((REF/'uniprot_sprot.fasta.gz').open('rb'),'sha256').hexdigest(),'taxonomy_sha256':hashlib.file_digest((REF/'taxdump.tar.gz').open('rb'),'sha256').hexdigest(),'parameters':'blastp and blastx; --sensitive --evalue 1e-5 --max-target-seqs 5 --threads 8','review_rule':'>=50% query coverage, >=40% identity, E<=1e-10, and same reference group among retained hits within 90% of best bitscore. Exploratory review flag, not a validated taxonomic classifier.','limits':'Biased reviewed database, at most five retained matches, homology cannot distinguish host contamination from symbionts, diet or other whole-body-associated organisms. No transcript removed. No-match does not mean contaminant-free. Not FCS-GX or comprehensive nucleotide taxonomy screening.','results':summaries}
(OUT/'broad_search_methods.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summaries,indent=2))
