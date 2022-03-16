library(data.table) 

setwd("~/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/03_MitochondrialProperties/StructuralInformation")

GeneOrder_aligned_with_tRNA <- read.delim("~/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/03_MitochondrialProperties/StructuralInformation/GeneOrder_aligned_with_tRNA.txt", row.names=1)
View(GeneOrder_aligned_with_tRNA)
GeneOrder_aligned_with_tRNA_NA <- na.omit (GeneOrder_aligned_with_tRNA, keep.rownames=TRUE)
GeneOrder_aligned_with_tRNA_NA_ID <- cbind(GeneOrder_aligned_with_tRNA_NA, "ID"=row.names(GeneOrder_aligned_with_tRNA_NA))

Counted_wtRNA <- setDT(GeneOrder_aligned_with_tRNA_NA)[, .N, by = c(names(GeneOrder_aligned_with_tRNA_NA))]
Counted_wtRNA <- Counted_wtRNA[order(N),]
View(Counted_wtRNA)
write.table(Counted_wtRNA, "Counts_wtRNAs.txt", sep="\t")

GeneOrder_aligned_without_tRNA <- read.delim("~/000_changed_documents/Analyses/AnnelidaMitochondrialGenomes/03_MitochondrialProperties/StructuralInformation/GeneOrder_aligned_without_tRNA.txt", row.names=1)
View(GeneOrder_aligned_without_tRNA)
GeneOrder_aligned_without_tRNA_NA <- na.omit (GeneOrder_aligned_without_tRNA, keep.rownames=TRUE)
GeneOrder_aligned_without_tRNA_NA_ID <- cbind(GeneOrder_aligned_without_tRNA_NA, "ID"=row.names(GeneOrder_aligned_without_tRNA_NA))

Counted_wotRNA <- setDT(GeneOrder_aligned_without_tRNA_NA)[, .N, by = c(names(GeneOrder_aligned_without_tRNA_NA))]
Counted_wotRNA <- Counted_wotRNA[order(N),]
View(Counted_wotRNA)
write.table(Counted_wotRNA, "Counts_wotRNAs.txt", sep="\t")

library(plyr)

TableMatch <- data.frame(Species=character(0), stringsAsFactors=FALSE)
for(i in 1:nrow(Counted_wtRNA)) {
  Found <- match_df(GeneOrder_aligned_with_tRNA_NA_ID, Counted_wtRNA[i,])
  Match <-apply( t(Found$ID) , 1, paste , collapse = " " )
  TableMatch <- rbind(TableMatch, Match)
}
write.table(TableMatch, "Counts_wtRNAs_MatchedSpecies.txt", sep="\t")

TableMatch <- data.frame(Species=character(0), stringsAsFactors=FALSE)
for(i in 1:nrow(Counted_wotRNA)) {
Found <- match_df(GeneOrder_aligned_without_tRNA_NA_ID, Counted_wotRNA[i,])
Match <-apply( t(Found$ID) , 1, paste , collapse = " " )
TableMatch <- rbind(TableMatch, Match)
}
write.table(TableMatch, "Counts_wotRNAs_MatchedSpecies.txt", sep="\t")
