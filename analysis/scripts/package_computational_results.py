from pathlib import Path
import csv,json,hashlib,gzip,shutil
r=Path(__file__).resolve().parents[2];o=r/'analysis/results';runtime=r/'.runtime'
# Preserve lineage-specific BUSCO summaries separately from the historical transcript-mode result.
bus=runtime/'busco/hymenoptera_proteins'
for ext in ['txt','json']:shutil.copyfile(next(bus.glob('short_summary.specific.*.'+ext)),o/('busco_hymenoptera_proteins_summary.'+ext))
shutil.copyfile(bus/'run_hymenoptera_odb12.2/full_table.tsv',o/'busco_hymenoptera_proteins_markers.tsv')
d=json.loads((o/'busco_hymenoptera_proteins_summary.json').read_text());(o/'busco_hymenoptera_methods.json').write_text(json.dumps({'software':'BUSCO 6.1.0; HMMER 3.4','mode':'proteins','lineage':d['lineage_dataset'],'input':'vendor_analysis/extracted_results/OPL_mRNA_3_24201_Analysis_Results/Transcriptome DeNovo Assembly/Dinoponera gigantea Assembly.pep','input_sha256':hashlib.file_digest((r/'vendor_analysis/extracted_results/OPL_mRNA_3_24201_Analysis_Results/Transcriptome DeNovo Assembly/Dinoponera gigantea Assembly.pep').open('rb'),'sha256').hexdigest(),'command':'busco --in <vendor peptides> --mode proteins --lineage_dataset hymenoptera_odb12.2 --cpu 8 --opt-out-run-stats','interpretation':'Completeness of delivered predicted proteins; different input and lineage from historical broad-eukaryote transcript-mode result. Do not treat as a controlled before/after comparison.'},indent=2)+'\n')
# Fresh database-wide DRAMP search; raw alignments are exported with a defined header.
header=['query_ID','DRAMP_ID','identity_percent','alignment_length','query_length','reference_length','query_start','query_end','reference_start','reference_end','E_value','bitscore','reference_title']
rows=[]
for line in (runtime/'dramp_all_proteins.tsv').open():
 f=line.rstrip('\n').split('\t');assert len(f)==13;rows.append(f)
with (o/'dramp_database_search.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(rows)
vendor={x['peptide_ID'] for x in csv.DictReader((o/'amp_integrated_review.tsv').open(),delimiter='\t')};independent={x['peptide_ID'] for x in csv.DictReader((o/'ampir_high_score_integrated_review.tsv').open(),delimiter='\t')};hit_ids={x[0] for x in rows}
clean=json.loads((runtime/'dramp_standardization.json').read_text());methods={'aligner':'DIAMOND 2.2.8','parameters':'blastp --ultra-sensitive --evalue 0.001 --max-target-seqs 10 --threads 6','reference':'DRAMP general AMP FASTA downloaded 4 November 2026; upstream path names DRAMP3.0_new; no unsupported release-number inference','reference_records_retained':clean['records_retained'],'reference_records_excluded_nonstandard_or_too_short':len(clean['excluded_IDs']),'reference_source_sha256':clean['source_sha256'],'alignment_rows':len(rows),'matched_protein_records':len(hit_ids),'vendor_candidate_records_with_match':len(hit_ids&vendor),'ampir_high_score_records_with_match':len(hit_ids&independent),'limits':'Reported alignments are similarity evidence only. E-values depend on this database and aligner. Short, modified and excluded reference peptides and masking/seed sensitivity limit detection. No-hit is not proof of no AMP activity or novelty.'}
(o/'dramp_database_search_methods.json').write_text(json.dumps(methods,indent=2)+'\n');(o/'dramp_reference_exclusions.json').write_text(json.dumps(clean,indent=2)+'\n');print(json.dumps(methods,indent=2))
# Final evidence tables join statistical search support without declaring biological activity.
for file in ['amp_vendor_candidates_final_review.tsv','ampir_high_score_integrated_review.tsv']:
 data=list(csv.DictReader((o/file).open(),delimiter='\t'))
 for x in data:
  hs=[h for h in rows if h[0]==x['peptide_ID']];x['fresh_DRAMP_search_match']=bool(hs);x['fresh_DRAMP_best_Evalue']=min((float(h[10]) for h in hs),default='NA')
 with (o/file).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
