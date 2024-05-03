## 0-rawdata/ATAC-seq/<cell_type>/<ENCODE_id>.bigBed

# diff analysis for "24" pair, i.e. CMP and Erythroblast. We will do analysis by DEseq2 and limma respectively.

# To facilitate latter process, we encode cell lines at first, and others can change the filename of input .tab file to 
# complete diff analysis for other pairs of data

path_ATAC = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq"
out_path="/storage/work/fjz5078/Courses/STAT555/project/2-out/ATAC-seq"

cell1_rep1 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/CMP_vs_ERY/CMP1.tab"
cell1_rep2 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/CMP_vs_ERY/CMP2.tab"
cell2_rep1 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/CMP_vs_ERY/ERY1.tab"
cell2_rep2 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/CMP_vs_ERY/ERY2.tab"

# 读取CSV文件
RNA_data <- read.csv("/storage/work/fjz5078/Courses/STAT555/project/2-out/RNA_seq/CMP_Erythroblast/CMP_Erythroblast_DGE.csv", header = TRUE)

# pkgs import
.libPaths("/storage/work/fjz5078/Courses/STAT555/project/4.2")

library(DESeq2)
library(limma)
library(dplyr)
library(tidyr)
library(edgeR)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnhancedVolcano)
library(ragg)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)


# Data preparations

files <- c(cell1_rep1, cell1_rep2, cell2_rep1, cell2_rep2)
samples <- c("Line1_Rep1", "Line1_Rep2", "Line2_Rep1", "Line2_Rep2")
column_names <- c("Region", "Size", "Covered", "Sum", "Mean0", "Mean")

data_list <- lapply(1:length(files), function(i) {
  dat <- read.table(files[i], header = FALSE, sep = "\t", check.names = FALSE,col.names = column_names)
  dat$sample <- samples[i]
  return(dat)
})

# combine all samples data

all_data <- bind_rows(data_list)


# 以Region为行，样本为列构建矩阵
count_data <- all_data %>%
  dplyr::select(Region, sample, Mean) %>%
  spread(key = sample, value = Mean) %>%
  replace(is.na(.), 0)  # 将NA替换为0，因为DESeq2需要完整的矩阵


region_names <- count_data$Region
# 对count_data矩阵中的Mean值进行取整

count_matrix <- as.matrix(count_data[, -1])
count_matrix <- round(count_matrix)  # Region列是第一列，对其他列进行取整
rownames(count_matrix) <- make.names(region_names)  

#filter > 15
rmean = rowMeans(count_data[,2:5])
keep = rmean > 15
count_matrix = count_matrix[keep,]

#kpt = rmean[keep]
#sum(keep)
#hist(kpt)
#hist(rmean)

sample_info <- data.frame(
  cell_type = rep(c("CMP", "Erythroblast"), each = 2),
  replicate = rep(1:2, 2),
  row.names = colnames(count_data)[-1]
)


# 使用DESeq2创建数据集
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ cell_type)

# 运行DESeq2进行标准化和差异表达分析
dds <- DESeq(dds)
# 过滤取整后为0 
dds <- dds[rowSums(counts(dds)) > 0,]
# 获取结果
results <- results(dds)


# 转换数据框并替换NA值
result_df <- as.data.frame(results)
result_df <- result_df[!is.na(result_df$padj),]  # 移除padj为NA的行
result_df$logP <- -log10(result_df$padj + 1e-300)  # 添加一个非常小的数避免log(0)

# 创建Volcano图
volcano_plot <- ggplot(result_df, aes(x=log2FoldChange, y=logP)) +
  geom_point(alpha=0.4) +
  theme_minimal() +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 Adjusted P-Value")

# 使用ggsave保存Volcano图
# 定义文件名
volcano_plot_file_name <- "Volcano_Plot_filtered.png"
# 创建一个基因与显著区域的位置数据框
gene_positions <- data.frame(
  Gene = nearest_gene_details$gene_name,
  Start = start(nearest_gene_details),
  End = end(nearest_gene_details),
  Chromosome = seqnames(nearest_gene_details)
)

significant_positions <- data.frame(
  Region = paste("Region", seq_along(start(significant_granges))),
  Start = start(significant_granges),
  End = end(significant_granges),
  Chromosome = seqnames(significant_granges)
)

# 创建图形对象
plot <- ggplot() +
  geom_segment(data = gene_positions, aes(x = Start, xend = End, y = Chromosome, yend = Chromosome), color = "blue", size=3) +
  geom_segment(data = significant_positions, aes(x = Start, xend = End, y = Chromosome, yend = Chromosome), color = "red", size =3) +
  theme_minimal() +
  labs(title = "Genes and Significant Regions on Chromosomes", x = "Genomic Position", y = "Chromosome")

# 保存图像，不需要显示
ggsave(file.path(out_path, "Gene_and_Significant_Region_Relationship_filtered.png"), plot = plot, width = 10, height = 8)
# 构建完整的文件保存路径
volcano_plot_path <- file.path(out_path, volcano_plot_file_name)


# 使用ggsave保存Volcano图
ggsave(volcano_plot_path, plot = volcano_plot, width = 10, height = 8, units = "in")

# 排除包含NA的行，然后筛选显著差异的区域 # Modify >2 
significant_regions <- results[!is.na(results$padj) & !is.na(results$log2FoldChange) & results$padj < 0.05 & abs(results$log2FoldChange) > 2, ]





# 设置GTF文件的路径
gtf_path <- paste0(path_ATAC, "/Mus_musculus.GRCm38.101.gtf")  # 文件路径

# 导入GTF文件
gtf <- rtracklayer::import(gtf_path)

# significant_regions的Region列包含位置信息如 "chr1_3048243_3048393"
# 分割出染色体和位置信息
region_info <- strsplit(as.character(rownames(significant_regions)), "_")
chromosomes <- sapply(region_info, function(x) x[1])
starts <- sapply(region_info, function(x) as.numeric(x[2]))
ends <- sapply(region_info, function(x) as.numeric(x[3]))

# 创建GRanges对象
significant_granges <- GRanges(
  seqnames = chromosomes,
  ranges = IRanges(start = starts, end = ends)
)

# region_filter = significant_regions@rownames
# count_data_filtered = count_data[count_data$Region %in% region_filter, ]
cts_anno = annotatePeak(significant_granges, TxDb= TxDb.Mmusculus.UCSC.mm10.knownGene,
                        tssRegion = c(-2000,2000), annoDb = 'org.Mm.eg.db')
labels = cts_anno@anno$SYMBOL



# EnhancedVolcano(significant_regions, x = 'log2FoldChange', y ='padj', lab = labels)

# 绘制饼图并赋值给变量
pie_chart <- EnhancedVolcano(significant_regions, x = 'log2FoldChange', y ='padj', lab = labels)

# 使用ggsave保存这个图表，前提是pie_chart是一个ggplot对象
ggsave(paste0(out_path,"EnhancedVolcano.png"), plot = pie_chart, width = 10, height = 8, dpi = 300)




# 将significant_granges中的染色体名称修改为不包含'chr'
seqlevelsStyle(significant_granges) <- "Ensembl"  # 假设significant_granges使用的是UCSC风格，转换为ENSEMBL风格

# 查找最近的基因
nearest_genes <- nearest(significant_granges, gtf)

# 获取这些基因的详细信息
nearest_gene_details <- gtf[nearest_genes]
# 写入CSV文件
write.csv(nearest_gene_details, paste0(out_path,"nearest_gene_details.csv"), row.names = FALSE)

# 提取基因ID
gene_ids <- sapply(elementMetadata(nearest_gene_details)$gene_id, function(x) strsplit(x, ";")[[1]][1])

# 去除重复的基因ID，因为可能有多个区域靠近同一个基因
gene_ids <- unique(gene_ids)
write.csv(gene_ids, paste0(out_path,"gene_ids.csv"), row.names = FALSE)


# library(Gviz)
# 
# # 遍历每个染色体，并逐个生成图像
# chromosomes <- levels(seqnames(significant_granges))
# 
# for (chrom in chromosomes) {
#   sig_tracks_chrom <- significant_granges[seqnames(significant_granges) == chrom]
#   if (length(sig_tracks_chrom) > 0) {
#     gene_tracks_chrom <- gtf[seqnames(gtf) == chrom]
# 
#     # 对于过于密集的区域，继续尝试仅绘制前10个显著区域
#     sig_tracks_chrom <- head(sig_tracks_chrom, 10)
# 
#     # 仅绘制显著区域中心附近的一小部分
#     center_positions <- (start(sig_tracks_chrom) + end(sig_tracks_chrom)) / 2
#     start_positions <- pmax(center_positions - 50000, 1)  # 向左扩展50 kb，确保不低于1
#     end_positions <- center_positions + 50000  # 向右扩展50 kb
# 
#     # 设置文件路径
#     file_path_svg <- file.path(out_path, paste0("Chromosome_", chrom, "_Visualization.svg"))
# 
#     # 开启SVG设备
#     svg(filename = file_path_svg, width = 20, height = 10)
# 
#     # 创建和设置图形轨道
#     genomeAxisTrack <- GenomeAxisTrack()
#     significantTrack <- AnnotationTrack(range = sig_tracks_chrom, name = "Significant Regions", col = "red")
#     geneTrack <- GeneRegionTrack(gene_tracks_chrom, name = "Genes", fill = "darkblue", transcriptAnnotation = "gene")
# 
#     # 为每个显著区域绘制独立图形
#     for (i in seq_along(sig_tracks_chrom)) {
#       plotTracks(list(genomeAxisTrack, significantTrack, geneTrack), chromosome = chrom, from = start_positions[i], to = end_positions[i])
#     }
# 
#     # 关闭SVG设备
#     dev.off()
#   }
# }


# 创建一个基因列表
# gene_list <- bitr(gene_ids, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene = labels,
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",
                ont = "BP",  # 可以改为"CC"或"MF"来探索其他类别
                pAdjustMethod = "BH",  # 调整p值的方法
                qvalueCutoff = 0.05,  # 设置显著性阈值
                readable = TRUE)  # 确保输出是易读的基因名称



# # 进行GO富集分析
# ego <- enrichGO(gene = gene_list$ENTREZID, 
#                 OrgDb = org.Mm.eg.db, 
#                 keyType = "ENTREZID",
#                 ont = "BP",
#                 pAdjustMethod = "BH",
#                 qvalueCutoff = 0.05,
#                 readable = TRUE)


# GO分析结果可视化
goplot <- barplot(ego, showCategory=10)

# 保存GO分析结果
ggsave(file.path(out_path, "GO_Enrichment_filtered.png"), goplot, width = 10, height = 8)


# 创建一个基因与显著区域的位置数据框
gene_positions <- data.frame(
  Gene = nearest_gene_details$gene_name,
  Start = start(nearest_gene_details),
  End = end(nearest_gene_details),
  Chromosome = seqnames(nearest_gene_details)
)

significant_positions <- data.frame(
  Region = paste("Region", seq_along(start(significant_granges))),
  Start = start(significant_granges),
  End = end(significant_granges),
  Chromosome = seqnames(significant_granges)
)

# 创建图形对象
plot <- ggplot() +
  geom_segment(data = gene_positions, aes(x = Start, xend = End, y = Chromosome, yend = Chromosome), color = "blue") +
  geom_segment(data = significant_positions, aes(x = Start, xend = End, y = Chromosome, yend = Chromosome), color = "red") +
  theme_minimal() +
  labs(title = "Genes and Significant Regions on Chromosomes", x = "Genomic Position", y = "Chromosome")

# 保存图像，不需要显示
ggsave(file.path(out_path, "Gene_and_Significant_Region_Relationship.png"), plot = plot, width = 10, height = 8)

significant_RNA_genes <- RNA_data[RNA_data$padj < 0.05, ]

significant_gene_expression <- RNA_data[RNA_data$Gene %in% gene_ids,]


significant_RNA_seq = significant_regions[nearest_gene_details$gene_name %in% RNA_data$X,]
intersection = nearest_gene_details$gene_name[nearest_gene_details$gene_name %in% RNA_data$X]#Benchmark

# 假设 RNA_data 的第一列名为 gene_name
# 假设 intersection 是你希望按顺序筛选和排序的基因名列表

RNA_data_ATAC <- RNA_data[RNA_data$X %in% intersection, ]

RNA_data_ATAC <- RNA_data_ATAC[match(intersection, RNA_data_ATAC$X), ]

correlation_coefficient_log2 <- cor(abs(RNA_data_ATAC$log2FoldChange), abs(significant_RNA_seq$log2FoldChange))
correlation_coefficient_p = cor(RNA_data_ATAC$padj, significant_RNA_seq$padj)

log2_frame = data.frame(Column1 = abs(RNA_data_ATAC$log2FoldChange), Column2 = abs(significant_RNA_seq$log2FoldChange))
p_frame = data.frame(Col1 = RNA_data_ATAC$padj, Col2 = significant_RNA_seq$padj )




p <- ggplot() +
  geom_point(data = log2_frame, aes(x = Column1, y = Column2), colour = "orange") +
  ggtitle(paste("Scatter Plot - Correlation: ", round(correlation_coefficient_log2, 2))) +
  xlab("Absolute value of RNA-seq-log2 Fold Change") +
  ylab("Absolute value of ATAC-seq log2 Fold Change ") +
  theme_minimal()

ggsave(paste0(out_path,"scatter_plot_log2.png"), plot = p, width = 8, height = 6, dpi = 300)

p_1 <- ggplot() +
  geom_point(data = p_frame, aes(x = log(Col1), y = log(Col2)), colour = "purple") +
  ggtitle(paste("Scatter Plot - Correlation: ", round(correlation_coefficient_p, 2))) +
  xlab("log value of RNA-seq adjusted p-value") +
  ylab("log value of ATAC-seq adjusted p-value") +
  theme_minimal()

ggsave(paste0(out_path,"scatter_plot_p1.png"), plot = p_1, width = 8, height = 6, dpi = 300)



