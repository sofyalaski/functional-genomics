### library biomaRt for the genome annotation, we have to choose the corresponding version of the mm9 genome ##

library(biomaRt)

ensembl67 <- useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')

### now get all (protein coding) genes                     

#prot.gene <- getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id","gene_biotype","strand", "chromosome_name","start_position", "end_position"), mart = ensembl67, filters = "biotype", values = "protein_coding")
prot.gene <- getBM(attributes=c("ensembl_gene_id", "gene_biotype", "strand",  "chromosome_name", "start_position", "end_position"), mart = ensembl67, filters = "biotype", values = "protein_coding")

### now the data frame prot.gene has to be converted into GenomicRanges object. The chromosome names and the strand has to be changed before, as an example, two simple functions....

convert.strand <- function(strand.col){
  xx <-strand.col
  xx <- as.factor(xx)
  levels(xx)[levels(xx)=="-1"] <- "-"
  levels(xx)[levels(xx)=="1"] <- "+"
  levels(xx)[!levels(xx)%in%c("-","+")] <- "*" 
  return(xx)
}

###

correct.chr <- function(genes){
  genes <- genes[genes$chromosome_name%in%c(1:19,"X","Y"),]
  genes$chromosome_name <- paste0("chr",genes$chromosome_name)
  return(genes)
}


library(GenomicRanges)

prot.gene$strand <- convert.strand(prot.gene$strand)
prot.gene <- correct.chr(prot.gene)
prot.gene.gr <- makeGRangesFromDataFrame(df=prot.gene, start.field="start_position", end.field="end_position", keep.extra.columns=TRUE,starts.in.df.are.0based=FALSE) #ensembl:1-based


## get the promoters with promoters function

upstr = 2000
downstr = 500
prot.gene.prom <- promoters(prot.gene.gr, upstream = upstr, downstream = downstr)

setwd("/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped/")
library(bamsignals)
print("get the readcounts from our different histone modifications")
MEF_H3K4me3_counts <- bamCount("SRR5077625_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_H3K4me2_counts <- bamCount("SRR5077629_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_H3K4me1_counts <- bamCount("SRR5077633_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_H3K9ac_counts <- bamCount("SRR5077637_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_H3K27ac_counts <- bamCount("SRR5077641_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_H3K27me3_counts <- bamCount("SRR5077645_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_H3K79me2_counts <- bamCount("SRR5077649_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_H3K36me3_counts <- bamCount("SRR5077653_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_H3K9me3_counts <- bamCount("SRR5077657_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_MNase_counts <- bamCount("SRR5077669_rmDup_sorted.bam", prot.gene.prom, verbose=F)
MEF_WCE_counts <- bamCount("SRR5077673_rmDup_sorted.bam", prot.gene.prom, verbose=F)

gene_id <- prot.gene.prom$ensembl_gene_id

hm_count <- data.frame(H3K4me3=MEF_H3K4me3_counts,
                       H3K4me2=MEF_H3K4me2_counts,
                       H3K4me1=MEF_H3K4me1_counts,
                       H3K9ac=MEF_H3K9ac_counts,
                       H3K27ac=MEF_H3K27ac_counts,
                       H3K27me3=MEF_H3K27me3_counts,
                       H3K79me2=MEF_H3K79me2_counts,
                       H3K36me3=MEF_H3K36me3_counts,
                       H3K9me3=MEF_H3K9me3_counts,
                       MNase=MEF_MNase_counts,
                       WCE=MEF_WCE_counts)

print("get the readcounts from our different histone modifications from ESC")
setwd("/project/functional-genomics/2019/data/sra/MEF_G3/prefetched/mapped_ESC/")
ESC_H3K4me3_counts <- bamCount("SRR5077628_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_H3K4me2_counts <- bamCount("SRR5077632_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_H3K4me1_counts <- bamCount("SRR5077636_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_H3K9ac_counts <- bamCount("SRR5077640_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_H3K27ac_counts <- bamCount("SRR5077644_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_H3K27me3_counts <- bamCount("SRR5077648_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_H3K79me2_counts <- bamCount("SRR5077652_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_H3K36me3_counts <- bamCount("SRR5077656_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_H3K9me3_counts <- bamCount("SRR5077660_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_MNase_counts <- bamCount("SRR5077672_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)
ESC_WCE_counts <- bamCount("SRR5077675_1_rmDup_sorted.bam", prot.gene.prom, verbose=F)

ESC_hm_count <- data.frame(H3K4me3=ESC_H3K4me3_counts,
                       H3K4me2=ESC_H3K4me2_counts,
                       H3K4me1=ESC_H3K4me1_counts,
                       H3K9ac=ESC_H3K9ac_counts,
                       H3K27ac=ESC_H3K27ac_counts,
                       H3K27me3=ESC_H3K27me3_counts,
                       H3K79me2=ESC_H3K79me2_counts,
                       H3K36me3=ESC_H3K36me3_counts,
                       H3K9me3=ESC_H3K9me3_counts,
                       MNase=ESC_MNase_counts,
                       WCE=ESC_WCE_counts)

normalize_signal <- function(sample_count, control_count){
	S <- (sample_count + 1)
	C <- (control_count + 1)
	m <- median(S/C)
	Snorm <- S/C * 1/m
	Snorm_log2 <- log(Snorm)
	Snorm_log2_scaled <- scale(Snorm_log2-mean(Snorm_log2) + 1, center=F, scale=T)
	return(list("log"= Snorm_log2, "scaled"= Snorm_log2_scaled))
}

normalize_signal <- function(sample_count, control_count){
  Snorm_log2 <- log(sample_count+1)
  Snorm_log2_scaled <- scale(Snorm_log2-mean(Snorm_log2) + 1, center=F, scale=T)
  return(list("log"= Snorm_log2, "scaled"= Snorm_log2_scaled))
}
print("normalize the signal for our different histone modifications")
hm_count_norm_log <- data.frame(H3K4me3=normalize_signal(MEF_H3K4me3_counts,MEF_MNase_counts)$log,
                       H3K4me2=normalize_signal(MEF_H3K4me2_counts,MEF_MNase_counts)$log,
                       H3K4me1=normalize_signal(MEF_H3K4me1_counts,MEF_MNase_counts)$log,
                       H3K9ac=normalize_signal(MEF_H3K9ac_counts,MEF_MNase_counts)$log,
                       H3K27ac=normalize_signal(MEF_H3K27ac_counts,MEF_MNase_counts)$log,
                       H3K27me3=normalize_signal(MEF_H3K27me3_counts,MEF_MNase_counts)$log,
                       H3K79me2=normalize_signal(MEF_H3K79me2_counts,MEF_WCE_counts)$log,
                       H3K36me3=normalize_signal(MEF_H3K36me3_counts,MEF_MNase_counts)$log,
                       H3K9me3=normalize_signal(MEF_H3K9me3_counts,MEF_WCE_counts)$log)
rownames(hm_count_norm_log) <- gene_id
hm_count_norm_log_sorted <- hm_count_norm_log[ order(rownames(hm_count_norm_log)),]

hm_count_norm_scaled <- data.frame(H3K4me3=normalize_signal(MEF_H3K4me3_counts,MEF_MNase_counts)$scaled,
                                H3K4me2=normalize_signal(MEF_H3K4me2_counts,MEF_MNase_counts)$scaled,
                                H3K4me1=normalize_signal(MEF_H3K4me1_counts,MEF_MNase_counts)$scaled,
                                H3K9ac=normalize_signal(MEF_H3K9ac_counts,MEF_MNase_counts)$scaled,
                                H3K27ac=normalize_signal(MEF_H3K27ac_counts,MEF_MNase_counts)$scaled,
                                H3K27me3=normalize_signal(MEF_H3K27me3_counts,MEF_MNase_counts)$scaled,
                                H3K79me2=normalize_signal(MEF_H3K79me2_counts,MEF_WCE_counts)$scaled,
                                H3K36me3=normalize_signal(MEF_H3K36me3_counts,MEF_MNase_counts)$scaled,
                                H3K9me3=normalize_signal(MEF_H3K9me3_counts,MEF_WCE_counts)$scaled)
rownames(hm_count_norm_scaled) <- gene_id
hm_count_norm_scaled_sorted <- hm_count_norm_scaled[ order(rownames(hm_count_norm_scaled)),]



ESC_hm_count_norm_log <- data.frame(H3K4me3=normalize_signal(ESC_H3K4me3_counts,ESC_MNase_counts)$log,
                                H3K4me2=normalize_signal(ESC_H3K4me2_counts,ESC_MNase_counts)$log,
                                H3K4me1=normalize_signal(ESC_H3K4me1_counts,ESC_MNase_counts)$log,
                                H3K9ac=normalize_signal(ESC_H3K9ac_counts,ESC_MNase_counts)$log,
                                H3K27ac=normalize_signal(ESC_H3K27ac_counts,ESC_MNase_counts)$log,
                                H3K27me3=normalize_signal(ESC_H3K27me3_counts,ESC_MNase_counts)$log,
                                H3K79me2=normalize_signal(ESC_H3K79me2_counts,ESC_WCE_counts)$log,
                                H3K36me3=normalize_signal(ESC_H3K36me3_counts,ESC_MNase_counts)$log,
                                H3K9me3=normalize_signal(ESC_H3K9me3_counts,ESC_WCE_counts)$log)
rownames(ESC_hm_count_norm_log) <- gene_id
ESC_hm_count_norm_log_sorted <- ESC_hm_count_norm_log[ order(rownames(ESC_hm_count_norm_log)),]

ESC_hm_count_norm_scaled <- data.frame(H3K4me3=normalize_signal(ESC_H3K4me3_counts,ESC_MNase_counts)$scaled,
                                   H3K4me2=normalize_signal(ESC_H3K4me2_counts,ESC_MNase_counts)$scaled,
                                   H3K4me1=normalize_signal(ESC_H3K4me1_counts,ESC_MNase_counts)$scaled,
                                   H3K9ac=normalize_signal(ESC_H3K9ac_counts,ESC_MNase_counts)$scaled,
                                   H3K27ac=normalize_signal(ESC_H3K27ac_counts,ESC_MNase_counts)$scaled,
                                   H3K27me3=normalize_signal(ESC_H3K27me3_counts,ESC_MNase_counts)$scaled,
                                   H3K79me2=normalize_signal(ESC_H3K79me2_counts,ESC_WCE_counts)$scaled,
                                   H3K36me3=normalize_signal(ESC_H3K36me3_counts,ESC_MNase_counts)$scaled,
                                   H3K9me3=normalize_signal(ESC_H3K9me3_counts,ESC_WCE_counts)$scaled)
rownames(ESC_hm_count_norm_scaled) <- gene_id
ESC_hm_count_norm_scaled_sorted <- ESC_hm_count_norm_scaled[ order(rownames(ESC_hm_count_norm_scaled)),]

library(biomaRt)
library(GenomicFeatures)
# get all exons corresponding to all ensembl genes #
mm9.exons <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                    dataset="mmusculus_gene_ensembl",
                    transcript_ids=NULL,
                    circ_seqs=DEFAULT_CIRC_SEQS,
                    filter=NULL,
                    id_prefix="ensembl_",
                    host="may2012.archive.ensembl.org", #www.ensembl.org",
                    port=80,
                    taxonomyId=NA,
                    miRBaseBuild=NA)
# now get the exons per gene (list of genomic ranges)
exonic <- exonsBy(mm9.exons, by="gene")
# reduce the exons by the union (list of genomic ranges)
red.exonic <- reduce(exonic)
# lengts of all genes as a sum of exons
exon.lengths <- sum(width(red.exonic))


setwd("/project/functional-genomics/2019/group3/rna_mapping/tableReads/")
library("DESeq2")
MEF_repA<- read.table("SRR5077600_ReadsPerGene.out.tab", row.names=1, skip=4)
MEF_repB<- read.table("SRR5077601_ReadsPerGene.out.tab", row.names=1, skip=4)
MEF_repC<- read.table("SRR5077602_ReadsPerGene.out.tab", row.names=1, skip=4)
MEF_PE<- read.table("SRR5077621_ReadsPerGene.out.tab", row.names=1, skip=4)
ESC_PE <- read.table("SRR5077624_ReadsPerGene.out.tab", row.names=1, skip=4)
ESC_SE <- read.table("SRR5077609_ReadsPerGene.out.tab", row.names=1, skip=4)

cts <- data.frame(MEF_PE=MEF_PE[,3], MEF_repA=MEF_repA[,3], MEF_repB=MEF_repB[,3], MEF_repC=MEF_repC[,3], ESC_PE=ESC_PE[,3], ESC_SE=ESC_SE[,3])
rownames(cts) <- rownames(MEF_PE)
cts <- cts[ order(rownames(cts)),]

coldata <- matrix(c(rep("MEF",4),rep("ESC",2),"paired-end",rep("single-read",3),"paired-end","single-end" ), nrow=6)
colnames(coldata) <- c("tissue", "type")
rownames(coldata) <- colnames(cts)
print("construct dds")
dds <- DESeqDataSetFromMatrix(countData=cts, colData=DataFrame(coldata), design= ~tissue)
dds <- estimateSizeFactors(dds)
### you can add this information to the dds object of DESeq2 ##
rowRanges(dds) <- red.exonic
## and then just use the function to calculate FPKM values 
print("calculate fpkm")
fpkm.dds <- fpkm(dds)
print("filter for protein coding genss")
prot.indices <- rownames(fpkm.dds) %in% prot.gene$ensembl_gene_id
fpkm.prot.dds <- fpkm.dds[prot.indices,]
print("calculate medians")
Y <- as.matrix(rowMedians(fpkm.prot.dds[,1:4]))#uses mean because we have 4 columns
rownames(Y) <- rownames(fpkm.prot.dds)
fpkm.prot.dds_wo_PE <- cbind(fpkm.prot.dds[,2:4])
Y_wo_PE <- as.matrix(rowMedians(fpkm.prot.dds_wo_PE)) #"real" median but without paired end data
rownames(Y_wo_PE) <- rownames(fpkm.prot.dds)
print("use normalized counts instead of fpkm")
norm.counts <- counts(dds,normalized=T)
norm.counts.prot <- norm.counts[prot.indices,]
Y_norm.count <- as.matrix(rowMedians(log(norm.counts.prot[,1:4]+1)))
rownames(Y_norm.count) = rownames(norm.counts.prot)

norm.counts.prot_wo_PE <- norm.counts.prot[,2:4]
Y_norm.count_wo_PE <- as.matrix(rowMedians(log2(norm.counts.prot_wo_PE+1)))
rownames(Y_norm.count_wo_PE) = rownames(norm.counts.prot)

print("calculate medians for ESC")
ESC_Y <- as.matrix(rowMeans(fpkm.prot.dds[,5:6]))#uses mean because we have 4 columns
rownames(ESC_Y) <- rownames(fpkm.prot.dds)
print("use normalized counts instead of fpkm for ESC")
ESC_Y_norm.count <- as.matrix(rowMedians(log(norm.counts.prot[,5:6]+1)))
rownames(ESC_Y_norm.count) = rownames(norm.counts.prot)




set.seed(1234)
random_indices <- sample(length(Y_wo_PE)) #random permutation of indices
print("constructing training data set")
training_indices <- random_indices[1:(length(Y_wo_PE)/2)]
training_X_log <- hm_count_norm_log_sorted[training_indices,]
training_X_scaled <- hm_count_norm_scaled_sorted[training_indices,]

ESC_training_X_log <- ESC_hm_count_norm_log_sorted[training_indices,]
ESC_training_X_scaled <- ESC_hm_count_norm_scaled_sorted[training_indices,]
print("constructing validation data set")
validation_indices <- random_indices[(length(Y_wo_PE)/2):((length(Y_wo_PE)/2)+(length(Y_wo_PE)/4))]
validation_X_log <- hm_count_norm_log_sorted[validation_indices,]
validation_X_scaled <- hm_count_norm_scaled_sorted[validation_indices,]

ESC_validation_X_log <- ESC_hm_count_norm_log_sorted[validation_indices,]
ESC_validation_X_scaled <- ESC_hm_count_norm_scaled_sorted[validation_indices,]
print("constructing test data set")
test_indices <- random_indices[((length(Y_wo_PE)/2)+(length(Y_wo_PE)/4)):length(Y_wo_PE)]
test_X_log <- hm_count_norm_log_sorted[test_indices,]
test_X_scaled <- hm_count_norm_scaled_sorted[test_indices,]

ESC_test_X_log <- ESC_hm_count_norm_log_sorted[test_indices,]
ESC_test_X_scaled <- ESC_hm_count_norm_scaled_sorted[test_indices,]

training_Y <- Y_wo_PE[training_indices]
validation_Y <- Y_wo_PE[validation_indices]
test_Y <- Y_wo_PE[test_indices]

ESC_training_Y <- ESC_Y[training_indices]
ESC_validation_Y <- ESC_Y[validation_indices]
ESC_test_Y <- ESC_Y[test_indices]

training_Y_norm.count <- Y_norm.count_wo_PE[training_indices]
validation_Y_norm.count <- Y_norm.count_wo_PE[validation_indices]
test_Y_norm.count <- Y_norm.count_wo_PE[test_indices]

ESC_training_Y_norm.count <- ESC_Y_norm.count[training_indices]
ESC_validation_Y_norm.count <- ESC_Y_norm.count[validation_indices]
ESC_test_Y_norm.count <- ESC_Y_norm.count[test_indices]


print("build linear model for log2 data using fpkm as Y")
data_log2 <- as.data.frame(cbind(training_Y,training_X_log))
lm.hist_mod.log2 <- lm(training_Y ~ H3K4me3 + H3K4me2 + H3K4me1 + H3K9ac + H3K27ac + H3K27me3 + H3K79me2 + H3K36me3 + H3K9me3 , data = data_log2)
summary(lm.hist_mod.log2)

print("build linear model for log2 data with normalized read count as Y")
data_log2_norm.count <- as.data.frame(cbind(training_Y_norm.count,training_X_log))
lm.hist_mod.log2_norm.count <- lm(training_Y_norm.count ~ H3K4me3 + H3K4me2 + H3K4me1 + H3K9ac + H3K27ac + H3K27me3 + H3K79me2 + H3K36me3 + H3K9me3 , data = data_log2_norm.count)
summary(lm.hist_mod.log2_norm.count)

print("build linear model for log2 data with normalized read count as Y for ESC")
ESC_data_log2_norm.count <- as.data.frame(cbind(ESC_training_Y_norm.count,ESC_training_X_log))
ESC_lm.hist_mod.log2_norm.count <- lm(ESC_training_Y_norm.count ~ H3K4me3 + H3K4me2 + H3K4me1 + H3K9ac + H3K27ac + H3K27me3 + H3K79me2 + H3K36me3 + H3K9me3 , data = ESC_data_log2_norm.count)
summary(ESC_lm.hist_mod.log2_norm.count)

print("build linear model for scaled data using fpkm as Y")
data_scaled <- as.data.frame(cbind(training_Y,training_X_scaled))
lm.hist_mod.scaled <- lm(training_Y ~ H3K4me3 + H3K4me2 + H3K4me1 + H3K9ac + H3K27ac + H3K27me3 + H3K79me2 + H3K36me3 + H3K9me3 , data = data_scaled)
summary(lm.hist_mod.scaled)

print("build linear model for scaled data with normalized read count as Y")
data_scaled_norm.count <- as.data.frame(cbind(training_Y_norm.count,training_X_scaled))
lm.hist_mod.scaled_norm.count <- lm(training_Y_norm.count ~ H3K4me3 + H3K4me2 + H3K4me1 + H3K9ac + H3K27ac + H3K27me3 + H3K79me2 + H3K36me3 + H3K9me3 , data = data_scaled_norm.count)
summary(lm.hist_mod.scaled_norm.count)

print("build linear model for scaled data with normalized read count as Y for ESC")
ESC_data_scaled_norm.count <- as.data.frame(cbind(ESC_training_Y_norm.count,ESC_training_X_scaled))
ESC_lm.hist_mod.scaled_norm.count <- lm(ESC_training_Y_norm.count ~ H3K4me3 + H3K4me2 + H3K4me1 + H3K9ac + H3K27ac + H3K27me3 + H3K79me2 + H3K36me3 + H3K9me3 , data = ESC_data_scaled_norm.count)
summary(ESC_lm.hist_mod.scaled_norm.count)


pred_log <- predict.lm(lm.hist_mod.log2_norm.count,newdata = as.data.frame(test_X_log))
names(test_Y_norm.count) <- rownames(test_X_log)
pred_matrix <- cbind(test_Y_norm.count,pred_log)
colnames(pred_matrix) <- c("true","predicted")
cor(pred_matrix)

ESC_pred_log <- predict.lm(ESC_lm.hist_mod.log2_norm.count,newdata = as.data.frame(ESC_test_X_log))
names(test_Y_norm.count) <- rownames(ESC_test_X_log)
ESC_pred_matrix <- cbind(test_Y_norm.count,ESC_pred_log)
colnames(ESC_pred_matrix) <- c("true","predicted")
cor(ESC_pred_matrix)

hist(log2(training_Y), main = "RNA-Seq")
hist(hm_count_norm_log_sorted$H3K4me3, main = "H3K4me3")
hist(hm_count_norm_log_sorted$H3K4me2, main = "H3K4me2")
hist(hm_count_norm_log_sorted$H3K4me1, main = "H3K4me1")
hist(hm_count_norm_log_sorted$H3K9ac, main = "H3K9ac")
hist(hm_count_norm_log_sorted$H3K27ac, main = "H3K27ac")
hist(hm_count_norm_log_sorted$H3K27me3, main = "H3K27me3")
hist(hm_count_norm_log_sorted$H3K79me2, main = "H3K79me2")
hist(hm_count_norm_log_sorted$H3K36me3, main ="H3K36me3")
hist(hm_count_norm_log_sorted$H3K9me3, main = "H3K9me3")

library(glmnet)

glmnet.hist_mod.scaled <- glmnet(y=training_Y_norm.count ,x = as.matrix(training_X_scaled), family="gaussian")
glmnet.hist_mod.scaled.cv <- cv.glmnet(y=training_Y_norm.count ,x = as.matrix(training_X_scaled), family="gaussian")
plot(glmnet.hist_mod.scaled, label=T)
plot(glmnet.hist_mod.scaled.cv)
