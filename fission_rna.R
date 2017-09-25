#source of inspiration: https://www.bioconductor.org/help/course-materials/2015/LearnBioconductorFeb2015/B02.1.1_RNASeqLab.html


library(fission)
## use data command to load a prepared SummarizedExperiment that was generated  and aligned
data(fission)

fission_data<-as.data.frame(assay(fission))
col_data<-as.data.frame(colData(fission))
row_data<-as.data.frame(rowData(fission))
colData(fission)

###view the different groups WT vs. Mutated
fission$strain

#confirm matrix is count of reads##
assayNames(fission)


######### Construct DESeqDataset Object for downstream analysis ################
library("DESeq2")
#design is the subgroups or categories for differential expression (including interaction effect)
fission_dsd<-DESeqDataSet(fission, design = ~ strain + minute + strain:minute)
nrow(fission_dsd)

#the rlog-transformed data transforms data to homoskedastic#
rld <- rlog(fission_dsd)
# access sample similarity using dist() function to obtain Euclidean distance between samples
#dist function expects the different samples to be rows of the matrix,
similarity_matrix<-dist(t(assay(rld)))

#visualize similarity/distance using heat map
library(pheatmap)
library("RColorBrewer")

dist_matrix <- as.matrix( similarity_matrix)
#rename rows to show treatment group description
rownames(dist_matrix ) <- paste( rld$strain, rld$minute, sep="-" )
colnames(dist_matrix ) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)
pheatmap(dist_matrix,
         clustering_distance_rows=similarity_matrix,
         clustering_distance_cols=similarity_matrix,
         col=colors)

library(ggplot2)

# intgroup are the groups of interests for labeling the samples by colors
pcadata <- plotPCA(rld, intgroup = c( "strain", "minute"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))
ggplot(pcadata, aes(PC1, PC2, color=minute, shape=strain)) + geom_point(size=3) +xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


###diff.expression analysis###
diff_exp<-DESeq(fission_dsd)
report<-results(diff_exp)
report

######### Visualization of count for a specific gene ###################
#find the most differentially expressed gene candidate in time course based on adj. p-val

topGene<-rownames(report)[which.min(report$padj)]
#plotCounts function that takes DESeqDataSet as argument
topGeneCounts<-plotCounts(fission_dsd, gene=topGene, intgroup=c("strain", "minute"), returnData= TRUE)

ggplot(topGeneCounts, aes(x=minute, y=count, color=strain, group=strain)) +
  scale_y_log10() + geom_point(position=position_jitter(width=0.1,height=0), size=3)+geom_line()+
  ggtitle(sprintf("Count for top gene %s as function of time ", topGene))

###########  Cluster the top 10 genes (based on variance of rlogged counts) ###############
library("genefilter")
#rowVars calculates the row variance. In this case, the row represents the genes.
top10 <- head(order(rowVars(assay(rld)),decreasing=TRUE),10)
#mean center the genes to look at amount of deviation
select_10<-assay(rld)[top10,]
mc_top10<-select_10-rowMeans(select_10)
df <- as.data.frame(colData(rld)[,c("strain","minute")])
pheatmap(mc_top10, annotation = df)

######## time course analysis ############
#compare full model and reduced model via likelyhood ratio test. This tests for genes expressing strain-specific diff. over time course (strain:minute)
diff_timepoint<- DESeq(fission_dsd, test="LRT", reduced = ~ strain + minute)
#Use wald test to determine significance
timepoint_result<-results(diff_timepoint)
timepoint_result$symbol<-mcols(diff_timepoint)$symbol
head(timepoint_result)

head(timepoint_result[order(timepoint_result$padj),],4)

data <- plotCounts(diff_timepoint, which.min(timepoint_result$padj),
                   intgroup=c("minute","strain"), returnData=TRUE)
#name of the gene with the smallest padj val
genename<- rownames(timepoint_result[which.min(timepoint_result$padj),])

ggplot(data, aes(x=minute, y=count, color=strain, group=strain)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") + scale_y_log10()+ggtitle(sprintf("Count for top gene %s as function of time ", genename))

########## Obtain log2fold change at the individual timepoints.

resultsNames(diff_timepoint)
tp15 <- results(diff_timepoint, name="strainmut.minute15", test="Wald")[genename,]
tp30<-results(diff_timepoint, name="strainmut.minute30", test="Wald")[genename,]
tp60<-results(diff_timepoint, name="strainmut.minute60", test="Wald")[genename,]
tp120<-results(diff_timepoint, name="strainmut.minute120", test="Wald")[genename,]
tp180<-results(diff_timepoint, name="strainmut.minute180", test="Wald")[genename,]


log2_change<-c(tp15$log2FoldChange, tp30$log2FoldChange, tp60$log2FoldChange, tp120$log2FoldChange, tp180$log2FoldChange)
timepoint<-c(15,30,60, 120,180)
se<-c(tp15$lfcSE, tp30$lfcSE, tp60$lfcSE, tp120$lfcSE,tp180$lfcSE)
time_data<-data.frame(timepoint, log2_change,se)

ggplot(time_data,aes(x=timepoint, y=log2_change))+geom_point(color='#336666')+geom_errorbar(aes(ymin=time_course-se, ymax=time_course+se), width=.1,color='#336666')+geom_line(color='#336666')+
  ggtitle(sprintf("Log2 fold change for gene %s in mutant over time", genename))


############ Cluster significant Genes based on log2fold change values #######
weights<- coef(diff_timepoint)
colnames(weights)
#interested in significantly different genes
topGenes<-head(order(timepoint_result$padj),10)
#get the coefficients of the top genes
mat<-weights[topGenes, -c(1,2)]
thr <- 4
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),cluster_col=F)
