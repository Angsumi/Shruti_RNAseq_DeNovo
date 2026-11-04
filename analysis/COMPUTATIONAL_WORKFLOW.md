# Computational workflow and reproducibility

All searches were run locally. Large tools, models, downloaded references and raw search intermediates are in ignored `.runtime/`; the report provides compressed per-sequence summaries and methods/checksums. Prediction input FASTAs retain original IDs, with only terminal translation-stop symbols removed where required.

1. `python3 analysis/scripts/audit_deliverables.py`
2. `python3 analysis/scripts/build_candidate_evidence.py`
3. Run Phobius 1.01 standalone as recorded in `phobius_methods.json`; use the 64-bit decoder. Its short wrapper emits string-initialization warnings, so retain the raw decoder output and check record counts.
4. Run `recheck_amp_alignments.py <downloaded DRAMP general FASTA>` with Biopython 1.88. Pairwise rechecking is separate from statistical database-wide search.
5. `python3 analysis/scripts/shortlist_candidates.py` creates the historical secretion-based review tier. This is not the final biological candidate ranking.
6. Download reviewed Swiss-Prot FASTA/release metadata and NCBI taxonomy from their official sources. Build a DIAMOND database and run blastp on vendor peptides and blastx on vendor transcripts. Use `--sensitive --evalue 1e-5 --max-target-seqs 5 --threads 8`; output fields must be `qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle`. Save raw tables as `.runtime/swissprot_protein_hits.tsv` and `.runtime/swissprot_transcript_hits.tsv`. Run `summarize_broad_search.py`.
7. Run BLASTn 2.17.0+ against official UniVec with parameters in `univec_methods.json`, saving `.runtime/univec_hits.tsv` with the recorded fields. Run `summarize_univec.py`. The output is a local VecScreen-style review and does not constitute official FCS certification. No sequences are removed.
8. Install ampir 1.1.0 in `.runtime/R-library` and run `Rscript analysis/scripts/run_ampir.R`. The precursor model is appropriate for full protein input. Use `macrel_eligibility.tsv` for the separate Macrel 1.6.1 mature-peptide model scope; never classify long precursors as if they were validated mature peptides.
9. Run `derive_ampir_outputs.py` to extract all ampir records at score >=0.8 with terminal stops removed, then run Phobius on them. Run `integrate_final_evidence.py` to combine models and broad homology. Thresholds are exploratory; retain all numeric scores and sensitivity analysis.
10. Run protein-mode BUSCO 6.1.0/HMMER 3.4 using hymenoptera_odb12.2, separate from the historical transcript-mode run.
11. Build a separate DIAMOND database from standard-residue DRAMP general sequences and run blastp `--ultra-sensitive --evalue 0.001 --max-target-seqs 10 --threads 6`. Save `.runtime/dramp_all_proteins.tsv`; excluded DRAMP records and checksums are preserved. Run `package_computational_results.py` after the lineage-specific BUSCO run completes.
12. Run `describe_expression.py`, `prepare_publication_metadata.py` and `plot_analysis_summary.py`. Refresh `download_manifest.json` after all outputs have been finalized.

Scientific limits and non-computational dependencies are recorded in `publication/READINESS_AND_CLAIMS.md`. No raw-read quantification, experimental assay, host-only assembly reconstruction or accession deposition is implied by these commands.

## Official reference sources

- Swiss-Prot: `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`
- Release: same directory, `reldate.txt`
- Taxonomy: `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`
- UniVec: `https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec`
- DRAMP: `https://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.fasta`
- BUSCO: official data server and lineage version captured in the exported summary.
