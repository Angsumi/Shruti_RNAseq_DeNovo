# Analysis plan toward a research manuscript

Updated 4 November 2026. AZ1, AZ2 and AZ3 belong to the same group. There is no control/treatment contrast. Tissue, caste/sex, pooling, colony, independence, collection details and library preparation remain to be documented. No differential-expression hypothesis testing is planned.

## Completed work and available evidence

- Vendor: read QC and trimming, Trinity assembly, protein/functional annotation, reported mapping, per-sample counts and DRAMP matches. Original files remain downloadable from the report.
- Post-delivery: saved BUSCO assessment of the vendor FASTA (35.2% complete across 125 broad eukaryotic markers).
- New retained-file audit: all expression IDs match transcript records, counts are nonnegative integers without duplicate IDs, and all AMP queries match peptide records. The gzip count table is a union of delivered rows, not a complete or normalized expression matrix. Missing entries remain NA.
- New sequence audit: 318,934 distinct transcript IDs, no exact duplicate nucleotide sequences, N50 419 bp, median 375 bp. N50 alone does not establish assembly quality. All peptide IDs resolve to transcripts; 11,850 ORFs are marked complete in original FASTA headers.
- New AMP audit: 576/688 match rows have query fragments shorter than 20 aa (minimum four). This descriptive cutoff is a review flag, not an AMP classifier or acceptance threshold. Full sequences for all 370 query IDs were extracted unchanged.

## Ordered next work

| Stage | Work | Required evidence / completion criterion |
|---|---|---|
| 1. Resolve provenance | Record sample metadata; request exact mapping reference, aligner, parameters, complete count outputs and original report corrections from vendor. | Document which counts are transcript-level (delivered IDs include isoform suffixes), selection rules, and why summed delivered counts differ from reported mapped reads. Missing rows must not be assumed zero. |
| 2. Contamination assessment | Screen nucleotide assembly for adapters/vectors and use a documented broad taxonomic reference search; inspect ambiguous hits, coverage and competing target/non-target matches. | Versioned database, commands, complete hit tables and a retained/flagged decision table. Homology-description species labels are not contamination assignments. Review symbionts and conserved genes before removing sequences. Currently pending: tools and databases are unavailable here. |
| 3. Assembly quality review | Run BUSCO with the most appropriate available insect/ant lineage, documenting its release; investigate the longest contigs and fragmentation. | Preserve broad-lineage result; report each lineage separately. Decide on any filtering only after contamination, coding potential and sequence evidence are reviewed. No automatic reassembly. |
| 4. Annotation consolidation | Reconcile transcript→gene→ORF IDs; rerun or extend functional searches with pinned protein/domain databases, retaining scores, coverage, thresholds and unannotated sequences. | A unified annotation table, coverage denominators, exact software/database versions and source-linked methods. |
| 5. AMP candidate validation | Re-evaluate full precursor sequences and match coverage, use documented similarity thresholds and independent AMP predictors, predict signal peptides/topology and review mature-peptide boundaries. Deduplicate candidates with a reversible mapping. | An evidence-ranked candidate table; predicted secretion/activity remain predictions. Partial ORFs require special caution. Short sequence matches alone cannot support novelty or function claims. |
| 6. Descriptive expression | Use the delivered counts to describe observed values and overlap. Obtain complete abundance outputs before TPM, full-library normalization or absence claims. Map final candidates to observed counts with NA preserved. | Reproducible descriptive tables; distinguish transcript and gene levels; no group contrast, no DE p-values, no unjustified cross-sample abundance ranking. |
| 7. Biological interpretation | Select a focused transcriptome-resource or computational AMP-candidate question based on tissue and final evidence. Compare candidates with relevant ant homologues, align curated sequences and assess phylogenetic support only where informative. | Literature-backed context, documented homolog selection and evidence supporting each biological claim. No assumption that samples are venom glands. |
| 8. Validation and publication package | Plan independent biological validation for selected claims; AMP activity claims require appropriate experiments. Freeze code, inputs, checksums, methods, figures, limitations and deposition records. | Assembly/annotation QC, final candidate table, sample metadata, data availability, accession numbers and a claim-to-evidence checklist reviewed before manuscript drafting. |

## Raw reads: conditional dependency

Assembly, annotation, mapping and counting are already reported as done. Raw reads are needed for independent remapping/requantification or reassembly only if evidence gaps justify those tasks. The supplied SFTP account authenticated, but the specified Rawdata folder was absent when checked on 4 November 2026; no reads were downloaded. Recovering raw data also matters for archival availability. No credentials are recorded here.

## Planned manuscript outputs

1. Sample and library metadata table.
2. Assembly/annotation QC table and BUSCO figure with clear lineage labels.
3. Taxonomic screening and filtering decision table.
4. Functional annotation summary with defined denominators.
5. Descriptive expression/overlap table with missing-value policy.
6. AMP evidence table and selected precursor diagrams/alignments.
7. Reproducible supplementary files, methods and data accessions.

## Method references

- [Trinity assembly quality assessment](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment)
- [Trinity N50 and expression-aware statistics](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats)
- [BUSCO user guide](https://busco.ezlab.org/busco_userguide)
- [SignalP 6.0 documentation](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
- [NCBI TSA submissions](https://submit.ncbi.nlm.nih.gov/about/tsa/)
- [DESeq2 design and replication guidance](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

## Further completed checks (4 November 2026)

Full candidates were locally re-aligned to the 688 original DRAMP target pairs using Biopython 1.88/BLOSUM62, gap open -11, extend -0.5. This is a pairwise recheck, not a database-wide search or significance test. The current public general reference FASTA has 12,784 records; URL and checksum are retained. Unicode whitespace was removed from reference sequences during parsing.

Phobius 1.01 standalone produced predictions for all 370 candidates: 21 signal-positive, 18 signal-positive without TM helices, 12 signal-positive with complete ORFs, 64 with one or more TM helices. Only terminal stop symbols were removed from predictor input. Original files remain unchanged. The full web submission was rejected for size, so predictions were run locally. Next: cross-reference these predictions with complete ORFs, alignment coverage and expression evidence for manual candidate review.
