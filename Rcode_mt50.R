#setwd("C:/Users/Admin/Desktop/PhD/scRNAseq/Sherbrooke_scRNA seq/4mo_Hyp_SVZ_WT_vs_3xTg/35000_cells_Nov_2023")
library(dplyr)
library(Seurat)
library(patchwork)
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(scibetR)
library(scMCA)
library(tidyverse)
library(pheatmap)
#library(dyno)
library(data.table)
library(magrittr)
library(WriteXLS)
library(harmony)
# download the RDS object from GEO, series GSE269035
#scAD <- readRDS("C:/Users/Admin/Desktop/PhD/scRNAseq/Sherbrooke_scRNA seq/4mo_Hyp_SVZ_WT_vs_3xTg/35000_cells_Nov_2023")
scAD
#32039 features across 35000 samples
#filter cells that have at least more than 200 genes expressed
#filter genes that have at least expressed in more than 3 cells
selected_c <- WhichCells(scAD, expression = nFeature_RNA > 200)
selected_f <- rownames(scAD)[Matrix::rowSums(scAD) > 3]
scAD <- subset(scAD, features = selected_f, cells = selected_c)
scAD
#30204 features across 34646 samples
scAD[["percent.mt"]] <- PercentageFeatureSet(scAD, pattern="^mt") #calculate percent mt
#12x10 "Violin_plot_first_QC"
VlnPlot(scAD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.1)
scAD <- subset(scAD, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 50)
scAD
#30204 features across 23351 samples
scAD <- subset(scAD, subset = Sample_Name !="Multiplet")
scAD
#30204 features across 20923
scAD <- NormalizeData(scAD) #logNorm
scAD <- FindVariableFeatures(scAD, selection.method="vst", nfeatures=2000) #find Variable features
top10 <- head(VariableFeatures(scAD), 10) #top 10 features
plot1 <- VariableFeaturePlot(scAD)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T) #save this one
plot1 + plot2
all.genes <- rownames(scAD)
scAD <- ScaleData(scAD, features= all.genes) #scale data
scAD <- RunPCA(scAD, features=VariableFeatures(object=scAD)) #Run pca and export elbow plot: use n dimensions for further calculations

print(scAD[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scAD, dims = 1:2, reduction = "pca")
DimPlot(scAD, reduction = "pca")

DimHeatmap(scAD, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(scAD, dims = 1:20, cells = 500, balanced = TRUE)
scAD <- JackStraw(scAD, num.replicate = 100)
scAD <- ScoreJackStraw(scAD, dims = 1:20)
JackStrawPlot(scAD)
ElbowPlot(scAD)
scAD <- FindNeighbors(scAD, dims = 1:8)
scAD <- FindClusters(scAD, resolution = 0.5)
scAD <- RunUMAP(scAD, dims = 1:8)
scAD <- RunTSNE(scAD, dims = 1:8)
DimPlot(scAD, reduction = "umap", label=T)
DimPlot(scAD, reduction = "tsne")
pro.core <- function(scibet.core){
  cell.type <- unname(unlist(scibet.core[,1]))
  scibet.core <- as.data.frame(t(scibet.core[,-1]))
  colnames(scibet.core) <- cell.type
  return(as.matrix(scibet.core))
}
model <- readr::read_csv("/media/mariano/Elements1/Fernandes/Pratesi_single_cell/GSE109774_scibet_core.csv") 
model <- pro.core(model)
genes <- rownames(model[,-1])
querry <- GetAssayData(scAD, slot = "data") #sparse matrix
querryMatrix <- as.matrix(querry) #convert to dense matrix
querryMatrix_transposed <- as.data.frame(t(querryMatrix)) #convert to transposed data frame
prd <- LoadModel_R(model)
label <- prd(querryMatrix_transposed) #generate a vector of cell IDs
scAD[["scibet_20TM"]] <- label #attribute the IDs in the seurat object metadata
mca_result <- scMCA(scdata = querryMatrix, numbers_plot = 3)
#scMCA_vis(mca_result)
scAD[["MCA"]] <- mca_result$scMCA
Idents(scAD) <- "seurat_clusters"
scAD[["cell.ID"]] <- Idents(scAD)
Idents(scAD) <- "Sample_Name"
unique(scAD@meta.data[["Sample_Name"]])
scAD <- RenameIdents(scAD, "SampleTag13_flex" = "WT_Hyp", "SampleTag14_flex" = "3xTg_Hyp", "SampleTag15_flex" = "WT_Svz", "SampleTag16_flex" = "3xTg_Svz")
scAD[["Sample_Name"]] <- Idents(scAD)
Idents(scAD) <- "Sample_Name"
Det.cells <- CellsByIdentities(scAD, idents = c("WT_Hyp", "3xTg_Hyp", "WT_Svz", "3xTg_Svz"))
#10 (width) 6 (height)
DimPlot(subset(scAD, subset = Sample_Name %in% c("WT_Svz", "3xTg_Svz") ), group.by = "cell.ID", split.by = "Sample_Name")+NoLegend()
DimPlot(subset(scAD, subset = Sample_Name %in% c("WT_Hyp", "3xTg_Hyp") ), group.by = "cell.ID", split.by = "Sample_Name")+NoLegend()
#15 (width) 6 (height)
DimPlot(subset(scAD, subset = Sample_Name %in% c("WT_Hyp","3xTg_Hyp","WT_Svz","3xTg_Svz","Undetermined") ), group.by = "cell.ID", split.by = "Sample_Name")
scAD <- RenameIdents(scAD, "WT_Hyp" = "Hyp",
                     "WT_Svz" =  "Svz",
                     "3xTg_Hyp" = "Hyp",
                     "3xTg_Svz" = "Svz")
scAD[["zone"]] <- Idents(scAD)
Idents(scAD) <- "Sample_Name"
scAD <- RenameIdents(scAD, "WT_Hyp" = "WT",
                     "WT_Svz" =  "WT",
                     "3xTg_Hyp" = "3xTg",
                     "3xTg_Svz" = "3xTg")
scAD[["Genotype"]] <- Idents(scAD)
write.csv(table(scAD$Genotype, scAD$zone), "cell_distribution_mt50.csv")

plot <- DimPlot(scAD, group.by = "seurat_clusters")
HoverLocator(plot = plot, information = FetchData(scAD, vars = c("cell.ID","seurat_clusters", "scibet_20TM", "MCA")))
md <- scAD@meta.data %>% as.data.table
cl_0 <-md %>% filter(seurat_clusters==0) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_1 <-md %>% filter(seurat_clusters==1) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_2 <-md %>% filter(seurat_clusters==2) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_3 <-md %>% filter(seurat_clusters==3) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_4 <-md %>% filter(seurat_clusters==4) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_5 <-md %>% filter(seurat_clusters==5) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_6 <-md %>% filter(seurat_clusters==6) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_7 <-md %>% filter(seurat_clusters==7) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_8 <-md %>% filter(seurat_clusters==8) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_9 <-md %>% filter(seurat_clusters==9) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_10 <-md %>% filter(seurat_clusters==10) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_11 <-md %>% filter(seurat_clusters==11) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_12 <-md %>% filter(seurat_clusters==12) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_13 <-md %>% filter(seurat_clusters==13) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_14 <-md %>% filter(seurat_clusters==14) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_15 <-md %>% filter(seurat_clusters==15) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_16 <-md %>% filter(seurat_clusters==16) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
cl_17 <-md %>% filter(seurat_clusters==17) %>% select(scibet_20TM,MCA) %>% group_by(scibet_20TM,MCA) %>% summarise(total_count=n(), .groups = 'drop') %>%arrange(desc(total_count))
Idents(scAD) <- "seurat_clusters"
scAD <- RenameIdents(scAD, "0" = "Astrocyte_1",
                     "1" = "Astrocyte_2",
                     "2" = "Astrocyte_3",
                     "3" = "Microglia_1",
                     "4" = "Oligodendrocyte_1",
                     "5" = "Oligodendrocyte_2",
                     "6" = "Astro_Micro",
                     "7" = "GABA_nueron",
                     "8" = "Astrocyte_4",
                     "9" = "Astrocyte_5",
                     "10" = "Endothelial",
                     "11" = "Astro/Endo",
                     "12" = "Tanycyte",
                     "13" = "Astrocyte_6",
                     "14" = "Ependymocyte",
                     "15" = "Vascular_Leptomeningeal",
                     "16" = "Microglia_2",
                     "17" = "T_cell")
scAD[["cell.ID"]] <- Idents(scAD)
#16 (width) 6 (height) "All_conditions_UMAP_annotated"
DimPlot(subset(scAD, subset = Sample_Name %in% c("WT_Hyp","3xTg_Hyp","WT_Svz","3xTg_Svz","Undetermined") ), group.by = "cell.ID", split.by = "Sample_Name",label = F,label.size = 1)
#9 (width) 6 (height) "UMAP_with_annotated_clusters"
DimPlot(scAD, reduction = "umap", label=T)

Idents(scAD) <- "cell.ID"

#master table
table <- as.data.frame.matrix(table(scAD@active.ident,scAD$Sample_Name))
table <- table %>% mutate("Tot cells" = rowSums(across(where(is.numeric))))
table$"Tot cells tagged" <- table$`Tot cells`-table$Undetermined
table$"% tagged" <- (table$`Tot cells tagged`/table$`Tot cells`)*100
colnames(table)<-c("WT_Hyp","Hyp_3xTg","WT_Svz","Svz_3xTg","Undetermined","Tot cells","Tot cells tagged","% tagged")
table$"WT_Hyp % on tagged" <- (table$WT_Hyp/table$`Tot cells tagged`)*100
table$"Hyp_3xTg % on tagged" <- (table$Hyp_3xTg/table$`Tot cells tagged`)*100
table$"WT_Svz % on tagged" <- (table$WT_Svz/table$`Tot cells tagged`)*100
table$"Svz_3xTg % on tagged" <- (table$Svz_3xTg/table$`Tot cells tagged`)*100
table = table[,c("Tot cells", "Tot cells tagged", "Undetermined", "% tagged", "WT_Hyp", "WT_Hyp % on tagged","Hyp_3xTg", "Hyp_3xTg % on tagged", "WT_Svz","WT_Svz % on tagged", "Svz_3xTg", "Svz_3xTg % on tagged")]
table <- table[order(row.names(table)), ]
names(table)[names(table) == 'Undetermined'] <- 'Untagged'
WriteXLS(table,ExcelFileName="master_table_tagged_cell_types_4months_mt50.xls", col.names = T, row.names = T)

save.image(file="Data.RData")











