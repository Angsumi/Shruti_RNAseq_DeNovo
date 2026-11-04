"""Export saved vendor BUSCO categories without rerunning BUSCO."""
import csv, json
from pathlib import Path
root = Path(__file__).resolve().parents[1]
summary = next((root / 'data/qc/busco_results').glob('short_summary.specific.*.json'))
data = json.loads(summary.read_text())['results']
out = root / 'results/tables/busco_summary.csv'
out.parent.mkdir(parents=True, exist_ok=True)
with out.open('w', newline='') as handle:
    writer = csv.writer(handle)
    writer.writerow(['Category', 'BUSCO_groups', 'Percentage'])
    for label, count, percent in [('Complete single-copy', 'Single copy BUSCOs', 'Single copy percentage'), ('Complete duplicated', 'Multi copy BUSCOs', 'Multi copy percentage'), ('Fragmented', 'Fragmented BUSCOs', 'Fragmented percentage'), ('Missing', 'Missing BUSCOs', 'Missing percentage')]:
        writer.writerow([label, data[count], data[percent]])
