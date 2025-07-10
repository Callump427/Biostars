library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(GenomicFeatures)
library(AnnotationDbi)
library(ggplot2)
library(plotly)
library(biomaRt)
library(tibble)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
samples <- read.csv("sample_metadata.csv",row.names = 1)
head(samples)

samples$Condition <- as.factor(samples$Condition)
levels(samples$Condition)

samples$Condition <- relevel(samples$Condition, ref = "spf")

# Build vector of quant.sf file paths
files <- file.path("salmom_quant", paste0(rownames(samples), "_quant"), "quant.sf")
names(files) <- rownames(samples)

gene_map = read.csv("gene_underscore_map.csv",col.names = c('ENMUST','ENMUSG'))

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = gene_map,
  ignoreTxVersion = TRUE
)


deseq_dataset<- DESeqDataSetFromTximport(txi, colData = samples, design=~Condition)
#size factor
deseq_dataset <- estimateSizeFactors(deseq_dataset)

# normalises the data to remove systematic effects
normalizationFactors(deseq_dataset)
#Generates the Counts, just to look at first 5 rows
counts(deseq_dataset,normalized = TRUE)[1:5,]

#stabilise the variance of data for downstream PCA
vst = varianceStabilizingTransformation(deseq_dataset)

#Plot the PCA, can treat it like ggplot 2 plot
plotPCA(vst,intgroup='Condition') + theme_bw()


#Hierarchical clustering, need a distance matrix
d = assay(vst)
#transposed data using t
d = t(d)
d = dist(d)

h = hclust(d)
plot(h)

#K means clustering, doesnt give a visual output like hclust, kmeans is supervised
k = kmeans(t(assay(vst)),centers = 2)
# The 3 steps to DESeq2 analysis
#1 Normalisation (size effect)

#2 estimate dispersion
deseq_dataset = estimateDispersions(deseq_dataset)
plotDispEsts(deseq_dataset)


#3 apply statistics (Wald)
deseq_dataset = nbinomWaldTest(deseq_dataset)
result_table = results(deseq_dataset)
summary(result_table)

#make a data.frame instead of DataFrame
result_df = as.data.frame
#can plot normalised count data for a gene, this highlights a filtered outlier
plotCounts(deseq_dataset,gene="ENSMUSG00000086922",intgroup = "Condition")

#remove filtered data (From mean counts)
filter_df1 = result_df[complete.cases(result_df),]

#filter results based on p adj and log2fold changed
filter_df1$padj < 0.05
filter_df2 = filter_df1[filter_df1$padj <0.05,]
#get absolute value, magnitude of 1 in either direction
abs(filter_df2$log2FoldChange) > 1
filter_df3 = filter_df2[abs(filter_df2$log2FoldChange) > 1,]

# Plotting DE genes using MA and Volcano plot
plotMA(result_table)

#volcano Plot
# filered df1 useful because it has all values but NA' removed
# padj is -log10 to make smaller p values higher on y axis

filter_df1$test = filter_df1$padj < 0.05 & abs(filter_df1$log2FoldChange) > 1
filter_df1 = rownames_to_column(filter_df1,var ='ensgene')


g = ggplot(filter_df1, aes(x= log2FoldChange, y=-log10(padj))) +
  geom_point(aes(colour = test), size=1,alpha=0.3) +
  scale_colour_manual(values = c('black','red')) +
  geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw() +
  theme(legend.position ='none')

# makes an ineractive version of the plot
#ggplotly(g)


#Will use BiomaRt to get useful gene names 
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

annotation = getBM(attributes = c('ensembl_gene_id','chromosome_name','start_position',
                                  'end_position','strand','gene_biotype','external_gene_name',
                                  'description'),
                   filters = c('ensembl_gene_id'),
                   values = filter_df1$ensgene,
                   mart = ensembl)
annotated_df = left_join(filter_df1,annotation,
                         by = c('ensgene' = 'ensembl_gene_id'))



#Heatmap, filter the now annoated df
anno_df2 = annotated_df[annotated_df$padj < 0.05,]
anno_df3 = anno_df2[abs(anno_df2$log2FoldChange) > 1,]

degs = anno_df3$ensgene
vst_mat = assay(vst)
data_hm = vst_mat[degs,]
rownames(data_hm) = anno_df3$external_gene_name

#one-liner to make heatmap using base heatmap
#heatmap(data_hm)

annotation <- data.frame(Condition = samples$Condition)
rownames(annotation) <- rownames(samples)

# nice cuts for up/down
pheatmap(data_hm,fontsize = 5, scale = 'row', cutree_cols = 2,cutree_rows = 2,annotation_col = annotation)

## GO Enrichment # first use BioMart

ent_gene = getBM(attributes = c('entrezgene_id'),
                   filters = c('ensembl_gene_id'),
                   values = anno_df3$ensgene,
                   mart = ensembl)
# make it a vector
ent_gene = (ent_gene$entrezgene_id)
ent_gene = as.character(ent_gene)

## for the universe, get all the genes that were tested for DE
ent_universe = getBM(attributes = c('entrezgene_id'),
                 filters = c('ensembl_gene_id'),
                 values = annotated_df$ensgene,
                 mart = ensembl)
ent_universe = (ent_universe$entrezgene_id)
ent_universe = as.character(ent_universe)

#BP is biological process
ego = enrichGO(gene = ent_gene,
               OrgDb = org.Mm.eg.db ,
               ont = "BP",
               universe = ent_universe
               )
barplot(ego, showCategory = 15)
dotplot(ego)
cnetplot(ego)

#Kegg enrichment

ekg = enrichKEGG(
  gene = ent_gene,
  universe = ent_universe,
  organism = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
  
write_tsv(anno_df3,"filtered_DEG.txt")

