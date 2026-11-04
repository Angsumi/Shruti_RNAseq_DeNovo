# Dinoponera gigantea transcriptome analysis

Open [index.html](index.html) for the analysis report and downloads for samples AZ1, AZ2 and AZ3.

1. Read quality control (vendor report).
2. Trinity de novo assembly and coding-region prediction.
3. Structural and functional annotation.
4. Read mapping and per-sample expression counts.
5. Antimicrobial peptide screening.
6. Post-delivery BUSCO assessment of the vendor assembly.

## Files

- `vendor_analysis/extracted_report/`: original vendor report and assets.
- `vendor_analysis/extracted_results/`: original assembly, peptides, annotations, expression workbooks and AMP results.
- `vendor_analysis/OPL_mRNA_3_24201_Analysis_Results.zip`: original delivery archive.
- `vendor_analysis/data/qc/`: saved BUSCO outputs and run provenance.
- `vendor_analysis/results/tables/`: BUSCO marker table and category summary.
- `vendor_analysis/scripts/`: BUSCO command, environment and summary exporter.
- `vendor_analysis/report_assets/`: images used by the report.
- `download_manifest.json`: downloadable file sizes and SHA-256 checksums.

BUSCO reports 35.2% complete, 30.4% fragmented and 34.4% missing across 125 groups. No control/treatment design is recorded. The page documents missing inputs and inconsistencies in the original report. Historical commands and saved results do not imply a new analysis was run during this rebuild.

Serve this directory using `python3 -m http.server 8000` and open `http://localhost:8000`.

## Continuing analysis

`analysis/results/` contains new retained-file audits and extracted candidate sequences. Run `python3 analysis/scripts/audit_deliverables.py` to regenerate them. See [research plan](analysis/RESEARCH_PLAN.md) for ordered next steps and publication requirements. The report labels completed, reported and pending work separately.
