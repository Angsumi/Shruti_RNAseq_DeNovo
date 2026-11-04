"""Derive predictor inputs and score summaries without inventing mature processing boundaries."""
from pathlib import Path
import csv,gzip
r=Path(__file__).resolve().parents[2];o=r/'analysis/results'
with gzip.open(o/'ampir_all_proteins.tsv.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
vendor={x['peptide_ID'] for x in csv.DictReader((o/'amp_integrated_review.tsv').open(),delimiter='\t')}
sequences={}
for e in (r/'vendor_analysis/extracted_results/OPL_mRNA_3_24201_Analysis_Results/Transcriptome DeNovo Assembly/Dinoponera gigantea Assembly.pep').read_text().split('>')[1:]:
 lines=e.splitlines();sequences[lines[0].split()[0]]=''.join(lines[1:]).rstrip('*')
with (o/'ampir_vendor_candidates.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(x for x in rows if x['peptide_ID'] in vendor)
positive=[x for x in rows if x['prob_AMP']!='NA' and float(x['prob_AMP'])>=.8]
with (o/'ampir_high_score_sequences.faa').open('w') as f:
 for x in positive:f.write('>'+x['peptide_ID']+'\n'+sequences[x['peptide_ID']]+'\n')
with (o/'ampir_threshold_sensitivity.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['reporting_cutoff','all_protein_records_at_or_above','vendor_candidate_records_at_or_above','interpretation'])
 for c in [.5,.7,.8,.9,.95]:w.writerow([c,sum(float(x['prob_AMP'])>=c for x in rows if x['prob_AMP']!='NA'),sum(float(x['prob_AMP'])>=c for x in rows if x['peptide_ID'] in vendor and x['prob_AMP']!='NA'),'Exploratory score sensitivity; no cutoff validates biological activity'])
