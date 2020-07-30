#!/bin/Rscript

# Michael Mansfield, 2020
# michaeljamesmansfield@gmail.com
# Scripts for analyzing ACE2 and DPP4 phylogenies of NCBI orthologs.
# You probably shouldn't use this for anything since it's not very extensible.
# Send me an email if you want an explanation or a similar analysis.

library(ape)
library(phangorn)
library(viridis)

# Read in all the data. X.t is the tree, x.info is the NCBI .csv metadata,
# X.needle is the output of Needleman-Wunsch global alignment to the human
# gene. These are read in and processed before dealing with plotting the trees.
dpp4.t <- read.tree('./DPP4/RAxML/RAxML_bipartitions.DPP4.1')
dpp4.info <- read.csv('./DPP4/DPP4_orthologs.csv')
dpp4.needle <- read.delim('./DPP4/BLAST/needle.tb', header=F)
colnames(dpp4.needle) <- c("human", "seqB", "length", "pct_id", "pct_sim", "gaps", "score")
dpp4.info$human_pctid <- dpp4.needle[match(dpp4.info$RefSeq.Protein.accessions, dpp4.needle$seqB), 4]

ace2.t <- read.tree('./ACE2/RAxML/RAxML_bipartitions.ACE2.1')
ace2.info <- read.csv('./ACE2/ACE2_orthologs.csv')
ace2.needle <- read.delim('./ACE2/BLAST/needle.tb', header=F)
colnames(ace2.needle) <- c("human", "seqB", "length", "pct_id", "pct_sim", "gaps", "score")
ace2.info$human_pctid <- ace2.needle[match(ace2.info$RefSeq.Protein.accessions, ace2.needle$seqB), 4]

# make data frame for each gene with only overlapping species
ace2.info.shared <- ace2.info[which(ace2.info$Scientific.name %in% dpp4.info$Scientific.name),]
dpp4.info.shared <- dpp4.info[which(dpp4.info$Scientific.name %in% ace2.info$Scientific.name),]

# drop tree tips that are not shared
ace2.t.shared <- keep.tip(phy=ace2.t, tip=as.character(ace2.info.shared$RefSeq.Protein.accessions))
dpp4.t.shared <- keep.tip(phy=dpp4.t, tip=as.character(dpp4.info.shared$RefSeq.Protein.accessions))

# midpoint root and reorder tips by branch order
root_and_reorder <- function(tree){
  reformatted_tree <- ladderize(midpoint(tree))
  return(reformatted_tree)
}
ace2.t.2 <- root_and_reorder(ace2.t.shared)
dpp4.t.2 <- root_and_reorder(dpp4.t.shared)

# reorder the info data frames by tip order
reorder_by_phylogeny <- function(tree, tree_info){
  # filter out internal nodes from the second column of the edge matrix
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip,2]
  phylogenetic_order <- tree$tip.label[ordered_tips]
  reordered_info <- tree_info[match(phylogenetic_order, tree_info$RefSeq.Protein.accessions),]
  return(reordered_info)
}

dpp4.info.shared <- reorder_by_phylogeny(dpp4.t.2, dpp4.info.shared)
ace2.info.shared <- reorder_by_phylogeny(ace2.t.2, ace2.info.shared)

# add colour info by percent ID to human
brewer_ryb_10 <- c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")
# colour ramp the lowest %id is ~30, so 70 colours is enough
ryb <- inferno(70)

pct_id_to_colour <- function(percent_identities, colour_palette){
  percent_identities <- round(as.numeric(percent_identities))
  difference_from_min <- percent_identities - 30
  colour_vector <- ryb[difference_from_min]
  return(colour_vector)
}
ace2.info.shared$colour <- pct_id_to_colour(ace2.info.shared$human_pctid, ryb)
dpp4.info.shared$colour <- pct_id_to_colour(dpp4.info.shared$human_pctid, ryb)

# rename tip labels into more human-readable names
rename_tips <- function(tree, tree_info) {
  # by default, turns the NCBI accession into 
  # accesion + species
  # you can include more or less columns here if you want
  renamed_tree <- tree
  #renamed_tree$tip.label <- paste( tree$tip.label, tree_info[match(tree$tip.label, tree_info$RefSeq.Protein.accessions),4])
  renamed_tree$tip.label <- as.character(tree_info[match(tree$tip.label, tree_info$RefSeq.Protein.accessions),4])
  return(renamed_tree)
}
ace2.t.3 <- rename_tips(ace2.t.2, ace2.info.shared)
dpp4.t.3 <- rename_tips(dpp4.t.2, dpp4.info.shared)

# so you can save your work
#write.tree(ace2.t.3, file='ACE2.ordered.shared.newick')
#write.tree(dpp4.t.3, file='DPP4.ordered.shared.newick')
#write.table(ace2.info.shared, file='ACE2.ordered.shared.tsv', sep='\t', row.names=F, quote=F)
#write.table(dpp4.info.shared, file='DPP4.ordered.shared.tsv', sep='\t', row.names=F, quote=F)

# making plots. First draw the tree, then annotate it.
pdf('tree_plots.pdf')
plot(ace2.t.3, main='ACE2', cex=0.2, align.tip.label = T, edge.color = '#595959')
nodelabels( node=1:ace2.t.3$Nnode+Ntip(ace2.t.3), pie=cbind(as.numeric(ace2.t.3$node.label),100-as.numeric(ace2.t.3$node.label)), piecol=c('cyan','magenta'),cex=0.25)
tiplabels(pch=16, col=ace2.info.shared[match(ace2.t.3$tip.label, ace2.info.shared$Scientific.name),]$colour, cex=.5)

add.scale.bar()
#plot(rep(1, times=nrow(ace2.info.shared)), 1:nrow(ace2.info.shared), pch=16, col=ace2.info.shared$colour, cex=0.3)

plot(dpp4.t.3, main='DPP4', cex=0.2, align.tip.label = T, edge.color='#595959')
nodelabels( node=1:dpp4.t.3$Nnode+Ntip(dpp4.t.3), pie=cbind(as.numeric(dpp4.t.3$node.label),100-as.numeric(dpp4.t.3$node.label)), piecol=c('cyan','magenta'),cex=0.25)
tiplabels(pch=16, col=dpp4.info.shared[match(dpp4.t.3$tip.label, dpp4.info.shared$Scientific.name),]$colour, cex=.5)
add.scale.bar()

#plot(rep(1, times=nrow(dpp4.info.shared)), 1:nrow(dpp4.info.shared), pch=16, col=dpp4.info.shared$colour, cex=0.3)

plot(rep(1,70), 1:70, pch=15, col=inferno(70))

dev.off()
