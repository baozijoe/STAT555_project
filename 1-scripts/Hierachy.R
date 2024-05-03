path_ATAC = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq"
out_path="/storage/work/fjz5078/Courses/STAT555/project/2-out/ATAC-seq/"

cell1_rep1 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/Hierachy/CMP1.tab"
cell1_rep2 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/Hierachy/CMP2.tab"
cell2_rep1 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/Hierachy/ERY1.tab"
cell2_rep2 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/Hierachy/ERY2.tab"
cell3_rep1 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/Hierachy/HSC1.tab"
cell3_rep2 = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/Hierachy/HSC2.tab"

# pkgs import
.libPaths("/storage/work/fjz5078/Courses/STAT555/project/4.2")

library(DESeq2)
library(limma)
library(dplyr)
library(tidyr)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)


# Data preparations

files <- c(cell1_rep1, cell1_rep2, cell2_rep1, cell2_rep2, cell3_rep1,cell3_rep2)
samples <- c("CMP_Rep1", "CMP_Rep2", "ERY_Rep1", "ERY_Rep2","HSC_Rep1","HSC_Rep2")
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

library(ggplot2)
library(ggdendro)

# 数据准备和标准化，同前
rmean = rowMeans(count_data[,2:7])
keep = rmean > 15
count_data = count_data[keep,]

sample_data <- t(count_data[-1])  # 删除Region列
row.names(sample_data) <- samples
dist_matrix <- dist(sample_data, method = "euclidean")

# 基于距离矩阵进行层次聚类
hc <- hclust(dist_matrix, method = "ward.D2")

# 创建热图来展示样本间的聚类关系
# 并在边缘显示样本间的聚类树状图
# 使用ComplexHeatmap包
library(ComplexHeatmap)
heatmap_data <- as.matrix(dist_matrix)

# 生成聚类热图和保存文件
pdf(file = paste0(out_path, "/sample_clustering_heatmap.pdf"))
Heatmap(heatmap_data,
        name = "Distance",
        clustering_distance_rows = dist_matrix,
        clustering_distance_columns = dist_matrix,
        cluster_rows = hc,
        cluster_columns = hc,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_dend_reorder = TRUE,
        column_dend_reorder = TRUE)
dev.off()

# 继续使用前面处理好的count_data


group_names <- sub("_Rep[12]", "", row.names(scores))  # 删除"_Rep1" 或 "_Rep2" 以得到细胞系的名称
samples_group <- factor(group_names)

# 将分组因子添加到`scores`数据框中
scores$Group <- samples_group

# 执行PCA绘图
p <- ggplot(scores, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 5) +  # 画点
  geom_text(aes(label = row.names(scores)), vjust = 3, color = "black") +  # 加黑色标签
  scale_color_manual(values = c("red", "green", "blue")) +  # 手动指定颜色
  theme_minimal() +
  ggtitle("PCA of ATAC-seq Data") +
  xlab(paste("PC1 - Variance:", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 1), "%")) +
  ylab(paste("PC2 - Variance:", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 1), "%"))

# 保存图像到文件
# 保存图像到文件
ggsave(filename = paste0(out_path, "/pca_atac_seq.png"), plot = p, width = 10, height = 8)

