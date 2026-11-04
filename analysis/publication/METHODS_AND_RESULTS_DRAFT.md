# Computational methods and results draft

This is a source-linked draft for author review. It contains no experimental validation claims and no invented library methods.

## Study design and source material

AZ1, AZ2 and AZ3 represent three individual worker ants collected whole-body from the same colony at Mizoram University campus on 1 October 2025 (author-provided metadata). The reported identification is *Diacamma assamense*. There was no control/treatment comparison. These data provide within-colony observations and do not support generalization across colonies or treatment-effect testing. Library preparation, strandedness and RNA extraction details were not available in the retained materials. The report describes paired-end sequencing with 150-bp reads.

Original vendor deliverables comprised a transcript FASTA, predicted peptides, structural and functional annotation, per-sample read-count spreadsheets, and DRAMP match tables. Historical raw reads could not be retrieved from the supplied server because the specified data directory was absent at the time of checking. The assembly and count outputs were analyzed directly; assembly, mapping and counting were not repeated.

## Retained-file integrity and expression

Transcript and peptide FASTA identifiers were cross-referenced with the vendor workbooks. All 318,934 transcript IDs were distinct; no exact duplicate nucleotide sequences were found. The assembly comprised 149,464,614 bases, with GC content 44.24%, median sequence length 375 bp and N50 419 bp. These metrics alone do not establish assembly quality. All 41,281 predicted peptide identifiers resolved to transcript identifiers. ORF header annotations included 11,850 complete, 15,821 internal, 8,933 5-prime partial and 4,677 3-prime partial predictions.

AZ1, AZ2 and AZ3 workbooks contained 27,922, 18,957 and 15,377 distinct populated identifiers, respectively. All identifiers resolved to the assembly and retained counts were nonnegative integers, without duplicate IDs. The union comprised 28,141 identifiers; 12,833 appeared in all three tables. Missing workbook rows were retained as NA rather than imputed as zero. Count sums refer only to retained rows and are not sequencing-library totals. No TPM, library-wide normalization, differential-expression p-values or treatment contrasts were inferred.

## Expanded homology annotation and taxonomic review

DIAMOND 2.2.8 searched the predicted proteins (blastp) and transcript assembly (blastx) against UniProtKB/Swiss-Prot release 2026_03 using sensitive mode, E-value <=1e-5 and at most five target sequences per query. Reference and NCBI taxonomy checksums are provided in `broad_search_methods.json`. Best-hit protein descriptions are inferred annotations, not experimentally established functions.

Matches were reported for 25,999/41,281 proteins and 51,062/318,934 transcripts. An exploratory review rule used at least 50% query coverage, 40% identity and E-value <=1e-10, together with concordant reference groups among retained hits within 90% of the best bit score. This flagged 12,398 proteins and 6,909 transcripts for bacterial-homology review. Counts refer to sequence records, not unique genes or read abundance. Reviewed-database bias, truncated hit lists, whole-body-associated microbes, diet and conserved homology prevent interpreting these flags as definitive contamination assignments. No sequences were removed and no host-only cleaned assembly is claimed.

## Vector/adaptor review

NCBI BLASTn 2.17.0+ searched the nucleotide assembly against UniVec using the published VecScreen search parameters. Local score-based categories and proximity handling were implemented in the supplied script, with at most 100 target sequences per query. This is a local VecScreen-style screen, not an official FCS-adaptor/FCS-GX report. Moderate-or-strong matches occurred in 329 transcripts, including strong matches in 99. Biological sequence content in UniVec and possible legitimate homologues require manual boundary and library-context review before trimming. No sequence was automatically edited.

## AMP evidence review

The vendor table contained 688 match rows across 370 query IDs, all present in the predicted peptide FASTA. In 576 match rows, the reported query fragment was shorter than 20 amino acids; for 347 candidates all reported fragments were shorter than 20 amino acids. These are descriptive review flags, not validated AMP thresholds.

Full candidate sequences were re-aligned locally to their original DRAMP targets using Biopython 1.88, BLOSUM62, gap-open penalty 11 and extension penalty 0.5. The first optimal local alignment was retained and identity and query/target coverage were reported. These pairwise alignments do not provide E-values or constitute a database-wide novelty search.

Phobius 1.01 standalone predicted signal peptides for 21 of the 370 records; 18 had no predicted transmembrane helix and 12 signal-positive records had complete predicted ORFs. Eleven met the combined complete-ORF, signal-positive, no-transmembrane-helix rule. This earlier secretion-based tier was subsequently reviewed against broader annotation and independent model scores; it must not be interpreted as a validated AMP shortlist.

ampir 1.1.0's precursor SVM model evaluated all 41,281 predicted proteins after removing terminal translation-stop symbols. A reporting cutoff of 0.8 was used for exploratory review, with a sensitivity table across multiple cutoffs. Fifty records (46 unique sequences) met this cutoff across the protein set; none of the original 370 vendor candidates met it. This does not demonstrate absence of biological AMP activity. Macrel 1.6.1 was also run on the nine original records within its 10–100-aa length/residue domain; all nine returned negative predictions. Most full precursors fall outside Macrel's intended domain.

Phobius analysis of the 50 high ampir-score records identified six signal-positive records. Three were complete ORFs with a signal peptide and no predicted transmembrane helix, but they encoded one exact peptide sequence whose best Swiss-Prot match was phage lambda spanin. Host origin and antimicrobial activity are not established. No experimentally validated, novel ant-origin AMP is claimed.

## Conclusions supported at present

The retained assembly permits reproducible annotation and candidate auditing, but short AMP matches, mixed taxonomic homology and incomplete provenance limit biological interpretation. The appropriate current conclusion is that candidate evidence requires more careful validation; neither antimicrobial activity nor absence of AMP activity has been demonstrated.

## Primary sources and reproducibility

- [DIAMOND](https://github.com/bbuchfink/diamond)
- [UniProt data downloads](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/getting-data-from-uniprot/)
- [NCBI VecScreen parameters and interpretation](https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/)
- [Phobius publication](https://doi.org/10.1016/j.jmb.2004.03.016)
- [ampir documentation](https://ampir.marine-omics.net/)
- [Macrel model scope](https://macrel.readthedocs.io/en/latest/faq/)
- [BUSCO user guide](https://busco.ezlab.org/busco_userguide)

Executable analysis scripts, reference release information, per-sequence results and downloadable-file checksums accompany the report. Database and homology screens are limited to the references and parameters explicitly recorded; no comprehensive nucleotide taxonomy certificate is implied.

## Lineage-specific completeness and fresh DRAMP search

BUSCO 6.1.0 protein mode with HMMER 3.4 and hymenoptera_odb12.2 (release 13 May 2026) assessed 5,651 ortholog groups in the delivered predicted protein set: 629 complete (11.1%; 532 single-copy and 97 duplicated), 404 fragmented (7.1%), and 4,618 missing (81.7%). The input and lineage differ from the saved eukaryote transcript-mode run; these percentages should not be interpreted as a controlled before/after assembly comparison.

A separate DIAMOND 2.2.8 ultra-sensitive search compared all predicted proteins against 10,951 standard-residue records retained from the downloaded DRAMP general set, with E-value <=0.001 and at most ten targets. Unicode whitespace was removed; 1,833 nonstandard or too-short records were excluded and listed. The search reported 254 alignments across 115 protein records, including 26 of the original 370 vendor candidates. One of the 50 high ampir-score records also had a DRAMP match; that record had a partial ORF and a best Swiss-Prot reference to human S100-A9. This does not establish ant origin or experimental activity. Short-peptide sensitivity and excluded modified sequences limit interpretation of negative matches.
