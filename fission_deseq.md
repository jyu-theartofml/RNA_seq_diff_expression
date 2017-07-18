The fission dataset from Bioconductor is collected from an experiment where 2 groups of fission yeast, WT and mutant (*atf21* deletion) undergo oxidative stress. RNA-seq was performed for the two groups at 6 time points to measure gene expression levels.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### Reference: <https://f1000researchdata.s3.amazonaws.com/manuscripts/9994/54602c98-5b2c-4d9d-9004-4e1f2c4e52c6_7035_-_michael_love_v2.pdf?doi=10.12688/f1000research.7035.2>

``` r
library(fission)
data(fission)
```

### After importing the SummarizedExperiment dataset fission, view the columns and row values.

``` r
colData(fission)
```

    ## DataFrame with 36 rows and 4 columns
    ##              strain   minute replicate          id
    ##            <factor> <factor>  <factor> <character>
    ## GSM1368273       wt        0        r1     wt_0_r1
    ## GSM1368274       wt        0        r2     wt_0_r2
    ## GSM1368275       wt        0        r3     wt_0_r3
    ## GSM1368276       wt       15        r1    wt_15_r1
    ## GSM1368277       wt       15        r2    wt_15_r2
    ## ...             ...      ...       ...         ...
    ## GSM1368304      mut      120        r2  mut_120_r2
    ## GSM1368305      mut      120        r3  mut_120_r3
    ## GSM1368306      mut      180        r1  mut_180_r1
    ## GSM1368307      mut      180        r2  mut_180_r2
    ## GSM1368308      mut      180        r3  mut_180_r3

``` r
rowData(fission)
```

    ## DataFrame with 7039 rows and 2 columns
    ##               symbol        biotype
    ##          <character>       <factor>
    ## 1               tlh1 protein_coding
    ## 2        SPAC212.09c     pseudogene
    ## 3         SPNCRNA.70          ncRNA
    ## 4         SPAC212.12 protein_coding
    ## 5        SPAC212.04c protein_coding
    ## ...              ...            ...
    ## 7035 SPMITTRNATYR.01           tRNA
    ## 7036 SPMITTRNAILE.02           tRNA
    ## 7037            atp9 protein_coding
    ## 7038 SPMITTRNAGLU.01           tRNA
    ## 7039            cox2 protein_coding

### To view different groups in this experiment.

``` r
fission$strain
```

    ##  [1] wt  wt  wt  wt  wt  wt  wt  wt  wt  wt  wt  wt  wt  wt  wt  wt  wt 
    ## [18] wt  mut mut mut mut mut mut mut mut mut mut mut mut mut mut mut mut
    ## [35] mut mut
    ## Levels: wt mut

### Construct DESeqDataset Object for downstream analysis in DESeq2 library.

``` r
library("genefilter")
library(DESeq2)

#design refers to variables for differential expression (including interaction effect)
fission_dsd<-DESeqDataSet(fission, design = ~ strain + minute + strain:minute)
```

### the rlog-transformed function transforms count data to homoskedastic distribution.

``` r
rld <- rlog(fission_dsd)
```

### Calculate sample similarity using dist() function to obtain Euclidean distance between samples.

#### Note: dist() expects the different samples to be rows of the matrix.

``` r
similarity_matrix<-dist(t(assay(rld)))
```

### Visualize similarity/distance using heat map.

``` r
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
```

![Figure 1.Similarity matrix of samples. Value of 0 indicates the most similar](unnamed-chunk-4-1.png) \#\#\# Use PCA to look at data in low dimension - plotPCA() produces a Principal Component Analysis (PCA) plot of the counts in object.

``` r
library(ggplot2)
#intgroup is variables of interest

pcadata <- plotPCA(rld,intgroup = c( "strain", "minute"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))
ggplot(pcadata, aes(PC1, PC2, color=minute, shape=strain)) + geom_point(size=3) +xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
```

![Figure 2.PCA plot of the counts in dataset with respect to strain and time point](unnamed-chunk-5-1.png)

``` r
library("pcaExplorer")
groups=colData(rld)$strain
cols <- scales::hue_pal()(2)[groups]
#don't use row names, too many genes
genespca(rld,ntop=100,
         choices = c(1,2),
         arrowColors=cols,groupNames=groups,
         alpha = 0.3,
         useRownamesAsLabels=F,
         varname.size = 4
        )
```

![Figure 3. PCA plot of the top 100 genes selected by highest variances](unnamed-chunk-6-1.png)

### In the two PCA plots, there does not seem to be a difference between the mutant and WT representations, but the time point data show distinct clustering. This is because the first PCA plot corresponds to aggregation of all genes in the data. To look at the difference of mutant vs. WT, we can perform differential expression analysis on the count data matrix to gain more insights.

### Use likelyhood ratio test (LRT) to test for genes expressing strain-specific diff. over time course (i.e., interaction term strain:minute).

``` r
diff_timepoint<- DESeq(fission_dsd, test="LRT", reduced = ~ strain + minute)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
class(diff_timepoint)
```

    ## [1] "DESeqDataSet"
    ## attr(,"package")
    ## [1] "DESeq2"

``` r
#to show the data table
timepoint_result<-results(diff_timepoint)
class(timepoint_result)
```

    ## [1] "DESeqResults"
    ## attr(,"package")
    ## [1] "DESeq2"

``` r
timepoint_result$symbol<-mcols(diff_timepoint)$symbol
head(timepoint_result[order(timepoint_result$padj),],4)
```

    ## log2 fold change (MLE): strainmut.minute180 
    ## LRT p-value: '~ strain + minute + strain:minute' vs '~ strain + minute' 
    ## DataFrame with 4 rows and 7 columns
    ##               baseMean log2FoldChange     lfcSE      stat       pvalue
    ##              <numeric>      <numeric> <numeric> <numeric>    <numeric>
    ## SPBC2F12.09c  174.6712    -2.65671953 0.7522613  97.28339 1.974151e-19
    ## SPAC1002.18   444.5050    -0.05093214 0.2042995  56.95360 5.169552e-11
    ## SPAC1002.19   336.3732    -0.39274898 0.5734940  43.53391 2.879804e-08
    ## SPAC1002.17c  261.7731    -1.13876477 0.6061288  39.31584 2.051371e-07
    ##                      padj      symbol
    ##                 <numeric> <character>
    ## SPBC2F12.09c 1.334526e-15       atf21
    ## SPAC1002.18  1.747308e-07        urg3
    ## SPAC1002.19  6.489157e-05        urg1
    ## SPAC1002.17c 3.466817e-04        urg2

``` r
#plotCounts take DESeqdataset object
data <- plotCounts(diff_timepoint, which.min(timepoint_result$padj),
                   intgroup=c("minute","strain"), returnData=TRUE)
#name of the gene with the smallest padj val
genename<- rownames(timepoint_result[which.min(timepoint_result$padj),])

ggplot(data, aes(x=minute, y=count, color=strain, group=strain)) +
  geom_point() + stat_smooth(se=FALSE,method="loess") + scale_y_log10()+ggtitle(sprintf("Count for top gene %s as function of time ", genename))
```

![Figure 4. Gene SPBC2F12.09c had the most significant difference in expression due to strain condition over time](unnamed-chunk-7-1.png) \#\# Obtain log2fold change at the individual timepoints.

``` r
resultsNames(diff_timepoint)
```

    ##  [1] "Intercept"           "strain_mut_vs_wt"    "minute_15_vs_0"     
    ##  [4] "minute_30_vs_0"      "minute_60_vs_0"      "minute_120_vs_0"    
    ##  [7] "minute_180_vs_0"     "strainmut.minute15"  "strainmut.minute30" 
    ## [10] "strainmut.minute60"  "strainmut.minute120" "strainmut.minute180"

``` r
tp15 <- results(diff_timepoint, name="strainmut.minute15", test="Wald")[genename,]
tp30<-results(diff_timepoint, name="strainmut.minute30", test="Wald")[genename,]
tp60<-results(diff_timepoint, name="strainmut.minute60", test="Wald")[genename,]
tp120<-results(diff_timepoint, name="strainmut.minute120", test="Wald")[genename,]
tp180<-results(diff_timepoint, name="strainmut.minute180", test="Wald")[genename,]

log2_change<-c(tp15$log2FoldChange, tp30$log2FoldChange, tp60$log2FoldChange, tp120$log2FoldChange, tp180$log2FoldChange)
timepoint<-c(15,30,60,120,180)
se<-c(tp15$lfcSE, tp30$lfcSE, tp60$lfcSE, tp120$lfcSE,tp180$lfcSE)
time_data<-data.frame(timepoint, log2_change,se)

ggplot(time_data,aes(x=timepoint, y=log2_change))+geom_point(color='#336666')+geom_errorbar(aes(ymin=time_course-se, ymax=time_course+se), width=.1,color='#336666')+geom_line(color='#336666')+
  ggtitle(sprintf("Log2 fold change for gene %s between mutant and WT over time", genename))
```

![Figure 5.The Log2 change of the mutant group vs. WT for top gene SPBC2F12.09C](unnamed-chunk-8-1.png)

Cluster top10 significant Genes based on log2fold change values.
----------------------------------------------------------------

``` r
weights<- coef(diff_timepoint)
colnames(weights)
```

    ##  [1] "Intercept"           "strain_mut_vs_wt"    "minute_15_vs_0"     
    ##  [4] "minute_30_vs_0"      "minute_60_vs_0"      "minute_120_vs_0"    
    ##  [7] "minute_180_vs_0"     "strainmut.minute15"  "strainmut.minute30" 
    ## [10] "strainmut.minute60"  "strainmut.minute120" "strainmut.minute180"

``` r
#interested in significantly different genes
topGenes<-head(order(timepoint_result$padj),10)
#get the coefficients of the top genes
mat<-weights[topGenes, -c(1,2)]
thr <- 4
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),cluster_col=F)
```

![Figure 6. Heat map of gene clustering based on log2 change of gene expression](unnamed-chunk-9-1.png)

``` r
sessionInfo()
```

    ## R version 3.4.0 (2017-04-21)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 7 x64 (build 7601) Service Pack 1
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] C
    ## 
    ## attached base packages:
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] pcaExplorer_2.2.0          factoextra_1.0.4.999      
    ##  [3] ggbiplot_0.55              scales_0.4.1              
    ##  [5] plyr_1.8.4                 devtools_1.13.2           
    ##  [7] rmarkdown_1.6              DESeq2_1.16.1             
    ##  [9] genefilter_1.58.1          BiocInstaller_1.26.0      
    ## [11] ggplot2_2.2.1              RColorBrewer_1.1-2        
    ## [13] pheatmap_1.0.8             fission_0.110.0           
    ## [15] SummarizedExperiment_1.6.3 DelayedArray_0.2.7        
    ## [17] matrixStats_0.52.2         Biobase_2.36.2            
    ## [19] GenomicRanges_1.28.3       GenomeInfoDb_1.12.2       
    ## [21] IRanges_2.10.2             S4Vectors_0.14.3          
    ## [23] BiocGenerics_0.22.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] colorspace_1.3-2        rprojroot_1.2          
    ##  [3] htmlTable_1.9           XVector_0.16.0         
    ##  [5] base64enc_0.1-3         d3heatmap_0.6.1.1      
    ##  [7] topGO_2.28.0            DT_0.2                 
    ##  [9] ggrepel_0.6.5           bit64_0.9-7            
    ## [11] AnnotationDbi_1.38.1    codetools_0.2-15       
    ## [13] splines_3.4.0           doParallel_1.0.10      
    ## [15] geneplotter_1.54.0      knitr_1.16             
    ## [17] jsonlite_1.5            Formula_1.2-2          
    ## [19] gridBase_0.4-7          annotate_1.54.0        
    ## [21] cluster_2.0.6           GO.db_3.4.1            
    ## [23] png_0.1-7               graph_1.54.0           
    ## [25] shinydashboard_0.6.1    shiny_1.0.3            
    ## [27] compiler_3.4.0          httr_1.2.1             
    ## [29] GOstats_2.42.0          backports_1.1.0        
    ## [31] assertthat_0.2.0        Matrix_1.2-9           
    ## [33] lazyeval_0.2.0          limma_3.32.3           
    ## [35] acepack_1.4.1           htmltools_0.3.6        
    ## [37] tools_3.4.0             Category_2.42.1        
    ## [39] gtable_0.2.0            GenomeInfoDbData_0.99.0
    ## [41] reshape2_1.4.2          Rcpp_0.12.12           
    ## [43] NMF_0.20.6              iterators_1.0.8        
    ## [45] stringr_1.2.0           mime_0.5               
    ## [47] rngtools_1.2.4          XML_3.98-1.9           
    ## [49] shinyAce_0.2.1          zlibbioc_1.22.0        
    ## [51] shinyBS_0.61            RBGL_1.52.0            
    ## [53] SparseM_1.77            yaml_2.1.14            
    ## [55] curl_2.7                memoise_1.1.0          
    ## [57] gridExtra_2.2.1         pkgmaker_0.22          
    ## [59] biomaRt_2.32.1          rpart_4.1-11           
    ## [61] latticeExtra_0.6-28     stringi_1.1.5          
    ## [63] RSQLite_2.0             highr_0.6              
    ## [65] foreach_1.4.3           checkmate_1.8.3        
    ## [67] BiocParallel_1.10.1     rlang_0.1.1            
    ## [69] pkgconfig_2.0.1         bitops_1.0-6           
    ## [71] evaluate_0.10.1         lattice_0.20-35        
    ## [73] htmlwidgets_0.9         labeling_0.3           
    ## [75] bit_1.1-12              AnnotationForge_1.18.0 
    ## [77] GSEABase_1.38.0         magrittr_1.5           
    ## [79] R6_2.2.2                Hmisc_4.0-3            
    ## [81] DBI_0.7                 foreign_0.8-67         
    ## [83] withr_1.0.2             survival_2.41-3        
    ## [85] RCurl_1.95-4.8          nnet_7.3-12            
    ## [87] tibble_1.3.3            locfit_1.5-9.1         
    ## [89] data.table_1.10.4       blob_1.1.0             
    ## [91] git2r_0.18.0            threejs_0.2.2          
    ## [93] digest_0.6.12           xtable_1.8-2           
    ## [95] tidyr_0.6.3             httpuv_1.3.5           
    ## [97] munsell_0.4.3           registry_0.3
