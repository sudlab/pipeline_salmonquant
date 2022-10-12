library(tximport)

samples <- list.dirs("quantification.dir", full.names=FALSE, recursive = FALSE)
file_list <- file.path("quantification.dir", samples, "quant.sf")
names(file_list) <- samples

tx2gene <- read.delim("transcript2geneMap.tsv")
summarised_to_transcript <- tximport(files = file_list,
                                     type="salmon",
                                     txOut=TRUE)
save(summarised_to_transcript, file="transcripts_tximport.RData")
summarised_to_gene <- summarizeToGene(summarised_to_transcript, tx2gene = tx2gene)
save(summarised_to_gene, file="gene_tximport.RData")
