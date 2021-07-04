if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ape')) install.packages('ape'); library('ape')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')

args <- commandArgs(trailingOnly=TRUE)
fasta_file <- args[1]
#fasta_file <- "longcovid.fasta"

dna<-ape::read.dna(fasta_file,format="fasta", as.matrix = TRUE)


dna.percentage.d = 100*(ape::dist.dna(dna, model='raw', pairwise.deletion = F, as.matrix=T))
dna.nt.d = ape::dist.dna(dna, model='N', pairwise.deletion = F, as.matrix=T)
dna.indel.d = ape::dist.dna(dna, model='indel', pairwise.deletion = F, as.matrix=T)

dna.d.melt <- reshape2::melt(dna.percentage.d)
dna.nt.melt<- reshape2::melt(dna.nt.d)
dna.indel.melt<- reshape2::melt(dna.indel.d)

dna.d.melt$valueCAT<- cut(dna.d.melt$value, breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25,0.3), right = FALSE)

pdf(paste(fasta_file, "pairwise.pdf", sep="_"))
ggplot(dna.nt.melt, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), color = 'white')+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, size = 7), axis.text.y = element_text(size = 7))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ scale_fill_gradient(low = "white", high = "#7f0bb0")+ggtitle("Number of Nucleotide differences")+
  theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position="bottom")+ 
  geom_text(aes(label = round(value, 1)))

ggplot(dna.indel.melt, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), color = 'white')+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, size = 7), axis.text.y = element_text(size = 7))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ scale_fill_gradient(low = "white", high = "red")+ggtitle("Number of indel differences")+
  theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position="bottom")+ 
  geom_text(aes(label = round(value, 1)))

ggplot(dna.d.melt, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = valueCAT), color = 'black')+theme_minimal()+
  scale_fill_manual(values = c('#f0f9e8',
                               '#a8ddb5',
                               '#238b45',
                               '#4eb3d3',
                               '#2b8cbe',
                               '#08589e',
                               '#c994c7'), name = 'Sequence\ndissimilarity, %')+
  theme(axis.text.x = element_text(angle = 90, size = 6), axis.text.y = element_text(size = 6))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+ggtitle("Percentage of differences")+
  theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position="bottom")+ 
  geom_text(aes(label = round(value, 3)))
dev.off()

colnames(dna.d.melt) <- c("Sequence 1", "Sequence2","Percentage of differences")
colnames(dna.nt.melt) <- c("Sequence 1", "Sequence2","Number of SNPs")
colnames(dna.indel.melt) <- c("Sequence 1", "Sequence2","Number of deletions")
write.table(dna.d.melt, file=paste(fasta_file, "percentage.tsv", sep="_"), quote=FALSE, sep='\t', row.names = F)
write.table(dna.nt.melt, file=paste(fasta_file, "snps.tsv", sep="_"), quote=FALSE, sep='\t', row.names = F)
write.table(dna.indel.melt, file=paste(fasta_file, "indels.tsv", sep="_"), quote=FALSE, sep='\t', row.names = F)
