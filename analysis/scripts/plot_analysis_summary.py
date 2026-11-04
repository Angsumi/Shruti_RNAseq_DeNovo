from pathlib import Path
import json,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path(__file__).resolve().parents[2];o=r/'analysis/results';p=r/'analysis/publication'
d=json.loads((o/'broad_search_methods.json').read_text())['results'];cats=['arthropod_homology_support','bacterial_homology_review','non_animal_homology_review','mixed_or_limited_homology','no_reported_match'];labels=['Arthropod support','Bacterial review flag','Other non-animal review flag','Mixed / limited homology','No reported match'];colors=['#377d62','#c78052','#b49b52','#889caf','#d9dedb']
fig,axs=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
for ax,kind in zip(axs,['protein','transcript']):
 counts=[d[kind]['status_counts'].get(c,0) for c in cats];ax.barh(labels[::-1],counts[::-1],color=colors[::-1]);ax.set_title(f"{kind.capitalize()} records (n={d[kind]['queries']:,})");ax.set_xlabel('Sequence records, not read abundance');ax.spines[['top','right']].set_visible(False)
fig.suptitle('Swiss-Prot homology review — exploratory flags, not contamination fractions',fontsize=12);fig.savefig(p/'homology_review.png',dpi=200);plt.close(fig)
rows=list(csv.DictReader((o/'ampir_threshold_sensitivity.tsv').open(),delimiter='\t'));fig,axs=plt.subplots(1,2,figsize=(9,4),layout='constrained')
for ax,col,title in [(axs[0],'all_protein_records_at_or_above','All predicted proteins'),(axs[1],'vendor_candidate_records_at_or_above','370 vendor candidates')]:
 ax.plot([float(x['reporting_cutoff']) for x in rows],[int(x[col]) for x in rows],marker='o',color='#25634e');ax.set_title(title);ax.set_xlabel('ampir precursor score cutoff');ax.set_ylabel('Records at or above cutoff');ax.spines[['top','right']].set_visible(False)
fig.suptitle('Threshold sensitivity — model scores do not validate activity',fontsize=12);fig.savefig(p/'ampir_score_sensitivity.png',dpi=200);plt.close(fig)
