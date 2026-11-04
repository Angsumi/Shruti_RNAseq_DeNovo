from pathlib import Path
import csv,json,hashlib
r=Path(__file__).resolve().parents[2];out=r/'analysis/results'
evidence={x['peptide_ID']:x for x in csv.DictReader((out/'amp_candidate_evidence.tsv').open(),delimiter='\t')};rows=[]
for line in (out/'phobius_predictions.txt').read_text().splitlines()[1:]:
 parts=line.split();q,tm,sp=parts[:3];assert q in evidence;assert tm.isdigit() and sp in ('0','Y')
 e=evidence[q];rows.append([q,e['ORF_type'],sp=='Y',int(tm),' '.join(parts[3:]),'N-terminal sequence incomplete; negative SP result inconclusive' if e['ORF_type'] in ('internal','5prime_partial') else 'Sequence-based prediction only; secretion/activity not experimentally established'])
assert len(rows)==len(evidence)==len(set(x[0] for x in rows))==370
assert (out/'phobius_raw_output.txt').read_text().count('\n>')+1==370
with (out/'phobius_candidate_summary.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['peptide_ID','ORF_type','signal_peptide_predicted','transmembrane_helices','topology','interpretation']);w.writerows(rows)
d={'software':'Phobius 1.01; decodeanhmm 1.2211d; 64-bit standalone','service':'https://phobius.sbc.su.se/','citation':'https://doi.org/10.1016/j.jmb.2004.03.016','command':'perl /path/to/phobius/phobius.pl -short analysis/results/amp_candidates_prediction_input.faa','input_preparation':'Only terminal translation-stop symbols removed; original candidate FASTA retained unchanged','input_sha256':hashlib.file_digest((out/'amp_candidates_prediction_input.faa').open('rb'),'sha256').hexdigest(),'candidates':len(rows),'SP_predicted':sum(x[2] for x in rows),'SP_predicted_no_TM':sum(x[2] and x[3]==0 for x in rows),'SP_predicted_complete_ORF':sum(x[2] and x[1]=='complete' for x in rows),'one_or_more_TM':sum(x[3]>0 for x in rows),'limits':'Predictions do not prove secretion or AMP activity. Negative SP calls on missing N termini are inconclusive. Wrapper emitted uninitialized string warnings; 370 unique prediction rows verified against raw decoder record count.'}
(out/'phobius_methods.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d,indent=2))
