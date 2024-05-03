## 0-rawdata/RNA-seq/<cell_type>/<ENCODE_id>.tsv

path_RNA = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/RNA-seq/"
path_ATAC = "/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/"

cellType_list<-c("HSC","CMP","Erythroblast")

original_RNA_HSC = read.table("/storage/work/fjz5078/Courses/STAT555/project/0-rawdata/RNA-seq/HSC/ENCFF064MKY.tsv", header=TRUE, sep="\t", row.names=1)

# 初始化一个向量，用于保存每列是否为整数的结果
columns_are_integers <- logical(ncol(data))

# 遍历数据框的每一列
for(i in seq_along(data)) {
  # 检查第 i 列的所有值是否为整数
  columns_are_integers[i] <- all(data[[i]] == floor(data[[i]]))
}

# 输出结果
column_names <- names(data)
integers_info <- paste(column_names, ":", ifelse(columns_are_integers, "Integers", "Not Integers"))
print(integers_info)
