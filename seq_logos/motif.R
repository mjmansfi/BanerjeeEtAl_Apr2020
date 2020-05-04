#!/bin/Rscript

library(ggplot2)
library(ggseqlogo)

ace2_mammals <- read.delim('./ACE2.ordered.shared.motif_mammals.tsv', header=F)
ace2_nonmammals <- read.delim('./ACE2.ordered.shared.motif_others.tsv', header=F)

dpp4_mammals <- read.delim('./DPP4.ordered.shared.motif_mammals.tsv', header=F)
dpp4_nonmammals <- read.delim('./DPP4.ordered.shared.motif_others.tsv', header=F)

colnames(ace2_mammals) <- c('Species', 'Seq', 'Protein', 'Animal')
colnames(ace2_nonmammals) <- c('Species', 'Seq', 'Protein', 'Animal')
colnames(dpp4_mammals) <- c('Species', 'Seq', 'Protein', 'Animal')
colnames(dpp4_nonmammals) <- c('Species', 'Seq', 'Protein', 'Animal')

ace2_joined <- rbind(ace2_mammals, ace2_nonmammals)
dpp4_joined <- rbind(dpp4_mammals, dpp4_nonmammals)
ace2_joined$Seq <- as.character(ace2_joined$Seq)
dpp4_joined$Seq <- as.character(dpp4_joined$Seq)


ace2_split <- split(ace2_joined, f=ace2_joined$Animal)
dpp4_split <- split(dpp4_joined, f=dpp4_joined$Animal)

p1 <- ggseqlogo(lapply(ace2_split, function(x) x[,2]), seq_type='aa', nrow=2, method='bits')
p1 <- p1 + ggplot2::scale_x_continuous(breaks=1:15, labels=c("Q24", "T27", "F28", "D30", "K31", "H34", "E35", "E37", "D38", "Y41", "Q42", "Y83", "K353", "G354", "D355"))
ggsave(filename='motifs_ace2.pdf', p1)

p2 <- ggseqlogo(lapply(dpp4_split, function(x) x[,2]), seq_type='aa', nrow=2, method='bits')
p2 <- p2 + ggplot2::scale_x_continuous(breaks=1:11, labels=c("K267", "F269", "T288", "A289", "A291", "L294", "I295", "R317", "Y322", "R336", "Q344", "K392"))
ggsave(filename='motifs_dpp4.pdf', p2)

