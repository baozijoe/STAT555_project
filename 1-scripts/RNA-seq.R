## 0-rawdata/RNA-seq/<cell_type>/<ENCODE_id>.tsv


## 设置package目录
.libPaths("/storage/group/dxl46/default/private/fanzhang/R_session/R/x86_64-pc-linux-gnu-library/4.2")


## 下载包
# BiocManager::install('org.Mm.eg.db')             # 小鼠数据库
# BiocManager::install("biomaRt",lib="/storage/group/dxl46/default/private/fanzhang/R_session/R/x86_64-pc-linux-gnu-library/4.2")          # 差异分析包，可以转换transcript_id到Gene Symbol



cellType_list<-c("HSC","CMP","Erythroblast")

path_RNA = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/RNA-seq"
path_ATAC = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq"
out_path="/storage/work/fjz5078/Courses/STAT555/project/2-out"

## 读取文件
original_RNA_HSC_1 = read.table(paste0(path_RNA,"/HSC/ENCFF247FEJ.tsv"), header=TRUE, sep="\t")
original_RNA_HSC_2 = read.table(paste0(path_RNA,"/HSC/ENCFF064MKY.tsv"), header=TRUE, sep="\t")
original_RNA_CMP_1 = read.table(paste0(path_RNA,"/CMP/ENCFF691MHW.tsv"), header=TRUE, sep="\t")
original_RNA_CMP_2 = read.table(paste0(path_RNA,"/CMP/ENCFF623OLU.tsv"), header=TRUE, sep="\t")
original_RNA_Erythroblast_1 = read.table(paste0(path_RNA,"/Erythroblast/ENCFF342WUL.tsv"), header=TRUE, sep="\t")
original_RNA_Erythroblast_2 = read.table(paste0(path_RNA,"/Erythroblast/ENCFF858JHF.tsv"), header=TRUE, sep="\t")


# 检查行名是否相同
all(rownames(original_RNA_HSC_1) == rownames(original_RNA_HSC_2))
all(rownames(original_RNA_CMP_1) == rownames(original_RNA_CMP_2))

# 合并read count
ReadCountMatrix<-data.frame(
  HSC_1=original_RNA_HSC_1$expected_count,
  HSC_2=original_RNA_HSC_2$expected_count,
  CMP_1=original_RNA_CMP_1$expected_count,
  CMP_2=original_RNA_CMP_2$expected_count,
  Erythroblast_1=original_RNA_Erythroblast_1$expected_count,
  Erythroblast_2=original_RNA_Erythroblast_2$expected_count
)
rownames(ReadCountMatrix)<-original_RNA_CMP_1$gene_id

# 过滤表达量全为0
filtered_ReadCountMatrix <- ReadCountMatrix[rowSums(ReadCountMatrix != 0) > 0, ]

# Normalize gene_id (removing .x from ENSMUSGxxxxx.x)
library(stringr)
gene_id<-rownames(filtered_ReadCountMatrix)
gene_id <- data.frame(stringr::str_replace(gene_id, "\\.[0-9]+$", ""))
colnames(gene_id)<-"gene_id"
rownames(filtered_ReadCountMatrix)<-gene_id$gene_id

# 转换gene_id, symbol相同则raw counts相加
library(org.Mm.eg.db)
library(biomaRt)
mart <- useMart("ensembl","mmusculus_gene_ensembl")
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = gene_id,
                      mart = mart)
id_to_symbol <- setNames(gene_symbols$external_gene_name, gene_symbols$ensembl_gene_id)
# missing_genes <- gene_id[!gene_id$gene_id %in% gene_symbols$ensembl_gene_id,]
filtered_ReadCountMatrix$Gene_Symbol <- ifelse(is.na(id_to_symbol[rownames(filtered_ReadCountMatrix)]),
                                               rownames(filtered_ReadCountMatrix), id_to_symbol[rownames(filtered_ReadCountMatrix)])
aggregate_df <- aggregate(. ~ Gene_Symbol, data = filtered_ReadCountMatrix, sum)

final_ReadCountMatrix<-aggregate_df[,2:ncol(aggregate_df)]
rownames(final_ReadCountMatrix)<-aggregate_df$Gene_Symbol

save(final_ReadCountMatrix,file=paste0(path_RNA,"/Annotated_RNA-seq_ReadCount.Rdata"))

## Sample_info (design_matrix)
load(file=paste0(path_RNA,"/Annotated_RNA-seq_ReadCount.Rdata"))
sample_info <- data.frame(
  cell_type = rep(c("HSC", "CMP", "Erythroblast"), each = 2),
  replicate = rep(1:2, 3),
  row.names = colnames(final_ReadCountMatrix)
)


## DESeq2
library(DESeq2)
integer_counts <- round(final_ReadCountMatrix) # need integer input

dds <- DESeqDataSetFromMatrix(countData = integer_counts,
                              colData = sample_info,
                              design = ~ cell_type)
dds <- DESeq(dds)
# 过滤取整后为0
dds <- dds[rowSums(counts(dds)) > 0,]


# QC
# First, apply the transformation necessary for PCA
rlg_data <- assay(rlog(dds))  # Using rlog transformation
t_data <- t(rlg_data)
# Perform PCA
pca_res <- prcomp(t_data, scale. = TRUE)

# Calculate the proportion of variance explained by the first few principal components
prop_var <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

# Extract scores for the first two principal components
pca_data <- data.frame(PC1 = pca_res$x[, "PC1"], PC2 = pca_res$x[, "PC2"])
pca_data$cell_type <- dds$cell_type  # Ensure 'cell_type' is a column in the metadata of dds

# Creating the PCA plot using ggplot2
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_type)) +
  geom_point(aes(shape = cell_type), size = 3) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(title = "PCA by Cell Type",
       x = paste("PC1 (", format(round(prop_var[1], 2), nsmall = 2), "% variance explained)", sep=""),
       y = paste("PC2 (", format(round(prop_var[2], 2), nsmall = 2), "% variance explained)", sep=""))


# Save the ggplot-enhanced PCA plot with specific dimensions
ggsave(paste0(out_path,"/RNA_seq/DESeq2_PCA_plot.png"), plot = pca_plot)

### Cluster

### Compute pairwise correlation values
rld_cor <- cor(rlg_data)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames

### Plot heatmap
library(pheatmap)
heatmap_plot <- pheatmap(rld_cor, silent = TRUE)

# Use ggsave to save the plot
ggsave(paste0(out_path,"/RNA_seq/DESeq2_cluster_plot.png"), plot = heatmap_plot$gtable, width = 8, height = 6)


### MA plot


### analysis
CMP_Ery_DESeq2 <- results(dds, contrast=c("cell_type", "CMP", "Erythroblast"))

HSC_CMP_DESeq2 <- results(dds, contrast=c("cell_type", "HSC", "CMP"))


library(EnhancedVolcano)
DGE_analysis<-function(dds,contrast_1,contrast_2,out_path){
  # 获取 1 vs 2的结果
  result<-results(dds,contrast=c("cell_type", contrast_1, contrast_2))
  
  contrast<-paste0(contrast_1,"_",contrast_2)
  DGE <-subset(result,padj < 0.05 & (log2FoldChange > 2 | log2FoldChange < -2))    #这里我们提取出2倍差异表达、校正p值<0.05的显著差异表达基因
  DGE <- DGE[order(DGE$log2FoldChange),]
  write.csv(DGE,paste0(out_path,"/RNA_seq/",contrast,"/",contrast,'_DGE.csv'),row.names = TRUE)
  
  ## 提取结果
  result <- as.data.frame(results(dds,contrast=c("cell_type", contrast_1, contrast_2))) # results()从DESeq分析中提取出结果表
  
  result <- na.omit(result)
  result1 <- result[c(2,6)]
  
  ## 火山图 (down for less expressed in contrast_1)
  plt_enhanced<-EnhancedVolcano(result1,
                       lab = rownames(result1),
                       title=paste0(contrast_1," vs ",contrast_2),
                       x = 'log2FoldChange',
                       y = 'padj')
  ggsave(paste0(out_path,"/RNA_seq/",contrast,"/",contrast,"_EnhancedVolcano.png"),width=10,height=10,plot=plt_enhanced)
  result1$change = ifelse(result1$padj < 0.05 & abs(result1$log2FoldChange) >= 2, ifelse(result1$log2FoldChange> 2 ,'Up','Down'),'Stable')
  library('ggplot2')
  ggplot(result1, aes(x = log2FoldChange, y = -log10(padj),colour=change)) +
    geom_point(alpha=0.3, size=3.5) +
    scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
    labs(title=paste0('Volcano plot for ',contrast_1," and ",contrast_2),x="log2(fold change)",y="-log10 (padj)")+ # 坐标轴# 坐标轴和图标题title="Volcano plot",
    theme_bw()+ #去除背景色
    #xlim(-20, 20)+ #设置坐标轴范围
    #ylim(0,40)+
    # 辅助线
    geom_vline(xintercept=c(-2,2),lty=3,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
    theme(panel.grid = element_blank())+ #去除网格线
    theme(plot.title = element_text(hjust = 0.5,size=24),
          legend.position="bottom",
          legend.title = element_blank(),
          legend.text=element_text(size=18),
          legend.key.size = unit(1, 'cm'),
          legend.background = element_rect(fill="gray90", linetype="solid",colour ="gray"),
          axis.title.x =element_text(size=18),
          axis.title.y=element_text(size=18),
          axis.text=element_text(size=14,face = "bold")
    )
  ggsave(paste0(out_path,"/RNA_seq/",contrast,"/",contrast,"_volcano.png"),width=8,height=6)# 保存为png格式
}

DGE_analysis(dds,"CMP", "Erythroblast",out_path)
DGE_analysis(dds,"HSC", "CMP",out_path)


## Limma-voom
library(limma)
library(edgeR)
#
dge <- DGEList(counts=integer_counts)
dge <- calcNormFactors(dge)
v <- voom(dge, design = model.matrix(~ cell_type, data = sample_info))

# model fitting
fit <- lmFit(v, design = model.matrix(~ cell_type, data = sample_info))

# CMP vs Erythroblast
cont.matrix_CMP_Ery <- makeContrasts(cell_typeErythroblast, levels = colnames(fit))
fit2 <- contrasts.fit(fit, cont.matrix)

# ebayes
fit2 <- eBayes(fit2)
results <- topTable(fit2, adjust="BH", number=Inf)

# Filter results based on log-fold change and adjusted p-value
# For example, to find genes with absolute logFC > 1 and adjusted p-value < 0.05
selected_genes <- results[abs(results$logFC) > 2 & results$adj.P.Val < 0.05, ]

write.csv(selected_genes,paste0(out_path,"/RNA_seq/CMP_Erythroblast/CMP_Erythroblast_DGE_limma.csv"),row.names = TRUE)

## Compare DESeq2 and Limma-voom

# 1. Venn
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Extract gene names for Venn diagram
genes_limma <- rownames(selected_genes)
rownames(CMP_Erythroblast_DGE)<-CMP_Erythroblast_DGE$...1
genes_deseq2 <- CMP_Erythroblast_DGE$...1

# Create Venn diagram
venn.diagram(
  x = list(
    Limma = genes_limma,
    DESeq2 = genes_deseq2
  ),
  filename = paste0(out_path,"/RNA_seq/CMP_Erythroblast/testcompare_venn.png"),
  resolution = 300,
  cex=2,
  cat.cex=2,
  disable.logging = T
)


# scatter p-value
# Identify overlapped genes
common_genes <- intersect(genes_limma, genes_deseq2)

# Data for plot
plot_data <- data.frame(
  padj_limma = selected_genes[common_genes, "adj.P.Val"],
  padj_deseq = CMP_Erythroblast_DGE[common_genes, "padj"]
)
colnames(plot_data)<- c("padj_limma","padj_DESeq2")
# Scatter plot
library(ggplot2)
plt_scatter<-ggplot(plot_data, aes(x = -log10(padj_DESeq2), y = -log10(padj_limma))) +
  geom_point(alpha = 0.6) +
  labs(x = "-log10 p-adj (DESeq2)", y = "-log10 p-adj (Limma)") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")  # Line y=x for reference

ggsave(paste0(out_path,"/RNA_seq/CMP_Erythroblast/compare_scatter_pvalue.png"),width=8,height=6,plot = plt_scatter)

## hierarchical tree
# which kind of normalization

## GO analysis
library(org.Mm.eg.db)
library(clusterProfiler)
BiocManager::install("pathview")
BiocManager::install("enrichplot")

function_plot <- function(pair,path_RNA){
  DGE <- read.csv(paste0(out_path,"/RNA_seq/",pair,"/",pair,"_DGE.csv"),header=T,row.names=1)
  
  # we want the log2 fold change 
  original_gene_list <- DGE$log2FoldChange
  
  # name the vector
  names(original_gene_list) <- rownames(DGE)
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  # GO
  require(DOSE)
  library(topGO)
  library(stringr)
  ego <- enrichGO(gene         = rownames(DGE),  # Adjust this if using different column names
                  OrgDb        = org.Mm.eg.db,
                  keyType      = "SYMBOL",  # Change if using a different identifier type like "ENTREZID"
                  ont          = "BP",  # BP for Biological Process, can also use "CC" or "MF"
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2)
  GO_plot<-dotplot(ego)+scale_x_discrete(labels=function(x) str_wrap(ego, width=10))
  ggsave(paste0(out_path,"/RNA_seq/",pair,"/",pair,"_GO.png"),width=9,height=8,plot = GO_plot)
  
  # GSEA
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mm.eg.db, 
               pAdjustMethod = "none")
  
  
  GSEA_plot<-dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)+scale_x_discrete(labels=function(x) str_wrap(gse, width=10))
  ggsave(paste0(out_path,"/RNA_seq/",pair,"/",pair,"_GSEA.png"),width=9,height=8,plot = GSEA_plot)
  
}


function_plot("CMP_Erythroblast",path_RNA)
function_plot("HSC_CMP",path_RNA)



##
# When choosing between limma-voom and DESeq2 for analyzing differential expression between two different cell types like CMP and Erythroblast, considering the biological variability is indeed very important. Here are some factors that might influence your decision:
#   
# 1. Understanding Biological Variability
# Biological variability refers to the natural differences observed between biological samples, which can stem from genetic, environmental, and technical factors. For distinct cell types like CMP and Erythroblast, you might expect substantial intrinsic biological variability due to their different roles and states in hematopoiesis.
# 
# 2. Assessment of DESeq2 for High Variability
# Model Flexibility: DESeq2 models count data using a negative binomial distribution, which is well-suited for data with high variance and overdispersion typical in RNA-seq experiments. This method might be particularly advantageous when comparing different cell types where the expression profile differences are expected to be large and variable.
# Robustness to Variability: DESeq2’s ability to model and account for dispersion at the individual gene level makes it robust in scenarios where gene expression variability is high, potentially providing more reliable p-value estimates under such conditions.
# 3. Limma-voom in Context of Lower Variability
# Assumptions of Lower Variance: While limma-voom is very effective, particularly with larger datasets, it assumes a somewhat consistent variance across expression levels after the voom transformation. If the actual variance significantly deviates from these assumptions (e.g., very high biological variability), this might affect the accuracy of its estimates.
# Empirical Bayes Moderation: Limma uses empirical Bayes to shrink estimates of standard errors, which is particularly effective in reducing type I error rates in datasets where true differences exist but are subtle. This method can also be powerful when the sample sizes are relatively large.
# 4. Practical Considerations
# Sample Size: If you have a large dataset, limma-voom might still perform well despite the high variability, due to its effective handling of variance through empirical Bayes methods.
# Data Integration: If your analysis will involve integrating these results with other datasets or types of analysis where limma's framework might be beneficial (e.g., microarray data, complex experimental designs), limma could be preferred.
# 5. Using Both Tools for Comprehensive Analysis
# Given the strengths of each tool, one approach is to use both DESeq2 and limma-voom to analyze your data. This can provide a more comprehensive view of the results, highlighting genes consistently identified across both methods and potentially capturing different aspects of the biological variability.
# 
# Conclusion
# In your specific case, if the primary concern is the robust handling of high biological variability between CMP and Erythroblast, DESeq2 might indeed be a more suitable choice due to its flexibility in modeling dispersion. However, using both tools where feasible and comparing their outputs can provide a more nuanced understanding of the differential expression results, helping to reinforce findings that are robust across analytical methods.

