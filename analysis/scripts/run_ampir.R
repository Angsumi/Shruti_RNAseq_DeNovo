#!/usr/bin/env Rscript
.libPaths(c('.runtime/R-library', .libPaths()))
library(ampir)
input <- 'vendor_analysis/extracted_results/OPL_mRNA_3_24201_Analysis_Results/Transcriptome DeNovo Assembly/Dinoponera gigantea Assembly.pep'
d <- read_faa(input)
d$seq_aa <- sub('\\*+$', '', d$seq_aa)
p <- predict_amps(d, model='precursor', n_cores=4)
p$peptide_ID <- sub(' .*$', '', p$seq_name)
p$length_aa <- nchar(p$seq_aa)
p <- p[, c('peptide_ID','length_aa','prob_AMP')]
con <- gzfile('analysis/results/ampir_all_proteins.tsv.gz', 'wt')
write.table(p, con, sep='\t', quote=FALSE, row.names=FALSE, na='NA'); close(con)
writeLines(capture.output(sessionInfo()), 'analysis/results/ampir_session_info.txt')
cat('Proteins:', nrow(p), 'Predictions >=0.8:', sum(p$prob_AMP >= .8, na.rm=TRUE), '\n')
