options(future.globals.maxSize = 100 * 1024^3)
# 加载必要的R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(ggalluvial)

# 读取已注释的数据对象
hypo <- readRDS("E:/时空/安诺/脑/V5toV4/umap.celltype.V5toV4.rds")
tes <- readRDS("E:/时空/安诺/睾丸/2025/combind_sub/2000/cell_type/umap.celltype.V4.rds")
setwd("E:\\时空\\安诺\\hypo_testis_com")

# 为每个数据添加批次（来源）信息
hypo$dataset <- 'Hypothalamus'
tes$dataset <- 'Testis'

# 将两个数据对象合并
combined_obj <- merge(hypo, y = tes, add.cell.ids = c("Hypo", "Tes"))
# SCTransform 归一化和特征选择
combined_obj <- SCTransform(combined_obj, vars.to.regress = "percent.mt")

# PCA降维
combined_obj <- RunPCA(combined_obj, verbose = FALSE)

# 使用 Harmony 进行批次效应校正
combined_obj <- RunHarmony(combined_obj, group.by.vars = "dataset", assay.use = "SCT")

# UMAP降维和聚类分析
combined_obj <- RunUMAP(combined_obj, reduction = "harmony", dims = 1:30)
combined_obj <- FindNeighbors(combined_obj, reduction = "harmony", dims = 1:30)
combined_obj <- FindClusters(combined_obj, resolution = 0.5)
# 可视化整合和聚类结果
# 按数据集和聚类结果分别绘制UMAP图
p1 <- DimPlot(combined_obj, reduction = "umap", group.by = "dataset") + 
    ggtitle("UMAP by Dataset") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(combined_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
    ggtitle("UMAP by Clusters") + theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(combined_obj, reduction = "umap", group.by = "sample") + ggtitle("UMAP by sample") + theme(plot.title = element_text(hjust = 0.5))

# 将两张图并排显示
p1 + p2 + p3
plot <- p1 + p2 + p3
ggsave("UMAP by Dataset Clusters sample.png", plot, width = 12, height = 6)
# 保存整合后的 Seurat 对象，包含原始的 cell_type 注释
saveRDS(combined_obj, file = "integrated_hypo_tes_pre_annotated.rds")
# 可视化整合后的 UMAP，并按原始注释的细胞类型分组
p_celltype <- DimPlot(combined_obj, reduction = "umap", group.by = "cell_type", label = TRUE) + 
    ggtitle("UMAP by Original Cell Type") + 
    theme(plot.title = element_text(hjust = 0.5))
ggsave("UMAP_by_Original_CellType.png", p_celltype, width = 10, height = 8)

# 可视化按原始组织分组的 UMAP
p_dataset <- DimPlot(combined_obj, reduction = "umap", group.by = "dataset") + 
    ggtitle("UMAP by Original Dataset") + 
    theme(plot.title = element_text(hjust = 0.5))
ggsave("UMAP_by_Dataset.png", p_dataset, width = 8, height = 6)

# 将两张图并排显示
p_celltype | p_dataset
# 计算每个组织的细胞总数
cell_counts_by_dataset <- combined_obj@meta.data %>%
    count(dataset) %>%
    rename(total_count = n)

# 计算每个组织的每种细胞类型数量和百分比
cell_proportions <- combined_obj@meta.data %>%
    count(dataset, cell_type) %>%
    left_join(cell_counts_by_dataset, by = "dataset") %>%
    mutate(proportion = n / total_count * 100)

# 保存比例数据到CSV文件
write.csv(cell_proportions, "cell_proportions_by_tissue_pre_annotated.csv", row.names = FALSE)

# 打印比例表格
print(cell_proportions)

# 绘制堆叠柱状图
p_stacked <- ggplot(cell_proportions, aes(x = dataset, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "Cell Type Proportion in Hypothalamus and Testis",
        x = "Tissue",
        y = "Proportion (%)",
        fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("cell_proportion_stacked_barplot_pre_annotated.png", p_stacked, width = 8, height = 6)

# 绘制冲积图（Alluvial Plot）
p_alluvial <- ggplot(cell_proportions, 
                     aes(x = dataset, y = n, stratum = cell_type, 
                         alluvium = cell_type, fill = cell_type)) +
    geom_flow(stat = "alluvium", lode.guidance = "right", alpha = 0.5) +
    geom_stratum() +
    labs(
        title = "Cell Type Composition Flow between Tissues",
        x = "Tissue",
        y = "Number of Cells"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("cell_proportion_alluvial_plot_pre_annotated.png", p_alluvial, width = 8, height = 6)
