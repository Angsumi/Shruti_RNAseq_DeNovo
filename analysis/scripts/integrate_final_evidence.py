"""Integrate independent predictions, broad homology and observed counts; retain conflicting evidence."""
from pathlib import Path
import csv,gzip,re,json,hashlib
from collections import Counter
r=Path(__file__).resolve().parents[2];o=r/'analysis/results'
def tsv(path):
 opener=gzip.open if str(path).endswith('.gz') else open
 with opener(path,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))
def save(name,rows):
 with (o/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
ampir={x['peptide_ID']:x for x in tsv(o/'ampir_all_proteins.tsv.gz')};vendor=tsv(o/'amp_integrated_review.tsv');annot={x['query_ID']:x for x in tsv(o/'swissprot_protein_review.tsv.gz')};newp={}
for line in (o/'ampir_high_score_phobius.txt').read_text().splitlines()[1:]:
 f=line.split();newp[f[0]]=(int(f[1]),f[2]=='Y')
seq={};headers={}
for e in (r/'vendor_analysis/extracted_results/OPL_mRNA_3_24201_Analysis_Results/Transcriptome DeNovo Assembly/Dinoponera gigantea Assembly.pep').read_text().split('>')[1:]:
 lines=e.splitlines();q=lines[0].split()[0];seq[q]=''.join(lines[1:]).rstrip('*');headers[q]=lines[0]
assert len(set(x['peptide_ID'] for x in vendor)&set(ampir))==370
final=[]
for x in vendor:
 q=x['peptide_ID'];a=annot[q];prob=ampir[q]['prob_AMP'];flags=[]
 if a['review_status'] in ('bacterial_homology_review','non_animal_homology_review'):flags.append('non_host_homology_requires_review')
 if float(prob)<.8:flags.append('not_supported_at_ampir_0.8_review_cutoff')
 if x['ORF_type']!='complete':flags.append('partial_ORF')
 final.append({**x,'ampir_precursor_probability':prob,'SwissProt_review_status':a['review_status'],'SwissProt_reference_group':a['reference_group'],'SwissProt_reference_title':a['best_reference_title'],'updated_review_flags':';'.join(flags) or 'manual_review_required','updated_interpretation':'Earlier secretion-based tier is provisional; activity and ant origin not established'})
save('amp_vendor_candidates_final_review.tsv',final)
new=[]
for q,x in ampir.items():
 if x['prob_AMP']=='NA' or float(x['prob_AMP'])<.8:continue
 tm,sp=newp[q];a=annot[q];orf=re.search(r'ORF type:(\S+)',headers[q]).group(1);tid=re.sub(r'\.p\d+$','',q)
 new.append({'peptide_ID':q,'transcript_ID':tid,'length_aa':len(seq[q]),'ampir_precursor_probability':x['prob_AMP'],'ORF_type':orf,'signal_peptide_predicted':sp,'TM_helices':tm,'sequence_cluster_sha256':hashlib.sha256(seq[q].encode()).hexdigest(),'SwissProt_review_status':a['review_status'],'SwissProt_reference_group':a['reference_group'],'SwissProt_reference_title':a['best_reference_title'],'interpretation':'Computational high score; mature boundaries, host origin and activity unvalidated'})
save('ampir_high_score_integrated_review.tsv',new)
summary={'predictor':'ampir 1.1.0 precursor SVM model','scope':'All 41,281 predicted proteins; terminal stop symbols removed; no manual mature-peptide cleavage inferred','review_cutoff':0.8,'cutoff_interpretation':'Suggested exploratory reporting threshold; not calibrated specificity or experimental validation','predictions':len(ampir),'NA_predictions':sum(x['prob_AMP']=='NA' for x in ampir.values()),'all_protein_high_scores':len(new),'unique_high_score_sequences':len({x['sequence_cluster_sha256'] for x in new}),'vendor_370_high_scores':sum(float(x['ampir_precursor_probability'])>=.8 for x in final),'vendor_old_tier_A_non_host_homology_flags':sum(x['review_tier'].startswith('A:') and x['SwissProt_review_status'] in ('bacterial_homology_review','non_animal_homology_review') for x in final),'high_score_with_SP':sum(x['signal_peptide_predicted'] for x in new),'high_score_complete_SP_no_TM':sum(x['ORF_type']=='complete' and x['signal_peptide_predicted'] and x['TM_helices']==0 for x in new),'high_score_homology_statuses':dict(Counter(x['SwissProt_review_status'] for x in new)),'training_reference':'https://ampir.marine-omics.net/','limitations':'Machine-learning model score does not validate activity, mature processing, host origin or novelty. Model training overlap with AMP reference databases has not been audited. No experimental results are available.'}
(o/'independent_amp_methods.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
