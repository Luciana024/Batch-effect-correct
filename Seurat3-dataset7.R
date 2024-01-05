########################
library(devtools)
BiocManager::install("Seurat")
library(Seurat)
packageVersion('Seurat')
library(magrittr)
library(cowplot)
library(SingleCellExperiment)
rm(list=ls())
#settings

normData = T
Datascaling = T
regressUMI = F
min_cells = 0
min_genes =0
norm_method = "LogNormalize"
scale_factor = 10000
numVG = 300
nhvg = 2000
npcs = 20
visualize = T
plotout_dir="/Users/fengjiao.li/Desktop/Seurat3/"
outfile_prefix = "Dataset7"
save_obj = T
saveout_dir="/Users/fengjiao.li/Desktop/Seurat3/"

# load evaluation & preprocess utils files #
this.dir <- '/Users/fengjiao.li/Desktop/demo_Seurat3/'
setwd(this.dir)
read_dir ="/Users/fengjiao.li/Desktop/dataset7/"
utils_dir <- '/Users/fengjiao.li/Desktop/Seurat3/'
source(paste0(utils_dir,'call_seurat_3.R'))
b1_exprs_filename = "b1_exprs.txt"
b2_exprs_filename = "b2_exprs.txt"
b1_celltype_filename = "b1_celltype.txt"
b2_celltype_filename = "b2_celltype.txt"

batch_label = "batch"
celltype_label = "celltype"

########################
# read data from flat text

b1_exprs <- read.table(file = paste0(read_dir,b1_exprs_filename),sep="\t",header=T,row.names=1,check.names = F)
b2_exprs <- read.table(file = paste0(read_dir,b2_exprs_filename),sep="\t",header=T,row.names=1,check.names = F)
b1_celltype <- read.table(file = paste0(read_dir,b1_celltype_filename),sep="\t",header=T,row.names=1,check.names = F)
b2_celltype <- read.table(file = paste0(read_dir,b2_celltype_filename),sep="\t",header=T,row.names=1,check.names = F)

b1_celltype$cell <- rownames(b1_celltype)
b1_celltype <- b1_celltype[colnames(b1_exprs),]
b2_celltype$cell <- rownames(b2_celltype)
b2_celltype <- b2_celltype[colnames(b2_exprs),]
b1_metadata <- as.data.frame(b1_celltype)
b2_metadata <- as.data.frame(b2_celltype)
b1_metadata$batch <- 1
b2_metadata$batch <- 2
b1_metadata$batchlb <- 'Batch_1'
b2_metadata$batchlb <- 'Batch_2'

expr= cbind(b1_exprs, b2_exprs)
metadata = rbind(b1_metadata, b2_metadata)
colnames(metadata)[colnames(metadata) == "CellType"] <- "celltype"
expr<- expr[, rownames(metadata)]
set.seed(123)
data7<-sample(expr,size=2000,replace=FALSE)
source(paste0(utils_dir,'ComBat_seq_functions.R'))
myFilteredData <- filter_data_mtx1(data7, base_name, is_filter_cells=TRUE,min_genes=300, 
                                  is_filter_genes=TRUE, min_cells=20)
lengths_list <- lengths(myFilteredData)
all_equal <- all(lengths_list == lengths_list[1])
matrix_data <- do.call(rbind, myFilteredData)
row_max <- matrixStats::rowMaxs(matrix_data)
row_min = matrixStats::rowMins(matrix_data)
if(sum(row_max==row_min)>0){
  myFilteredData = myFilteredData[row_max != row_min,]
}
cells_use <- colnames(myFilteredData)
metadata<- metadata[cells_use,]
myFilteredData_matrix=as.matrix(myFilteredData)
sce_obj <- SingleCellExperiment(assays=list(counts=myFilteredData_matrix),colData=metadata)
cellind=colData(sce_obj)$celltype
batchind=colData(sce_obj)$batchlb
ssq1<- function(x,batch,cell){
  mod <- lm(x~as.factor(batch)+as.factor(cell))
  av <- anova(mod)
  result <- av$`Sum Sq`
  return(result)
}
expvar=t(apply(assay(sce_obj,"counts"),1,ssq1,batch=batchind,cell=cellind))
colnames(expvar)=c("batch","cell","residual")
ss=apply(expvar,2,sum) ## This is the sum of all SS across all the genes
ss
r_raw_all <- ss[2]/ss[1] ## This is the overall ratio between cell and batch 
r_raw_all
########################
# run pipeline
# 检查数据对象中的细胞标识
source(paste0(utils_dir,'call_seurat_3.R'))
seurat_obj <- CreateSeuratObject(counts =myFilteredData,meta.data = metadata)
expr_mat<-seurat_obj@assays$RNA@data
expr_mat
batch_list = seurat3_preprocess(
  expr_mat, metadata, 
  normData = normData, Datascaling = Datascaling, regressUMI = regressUMI, 
  min_cells = min_cells, min_genes = min_genes, 
  norm_method = norm_method, scale_factor = scale_factor, 
  numVG = numVG, nhvg = nhvg, 
  batch_label = batch_label, celltype_label = celltype_label)

batches = call_seurat3(batch_list, batch_label, celltype_label, npcs, plotout_dir = this.dir, saveout_dir = this.dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
batches
umapplot_filename <- "_UMAP_Plot"  # 这只是一个例子，您可以根据自己的需要进行修改

p21 <- DimPlot(object = batches, reduction= 'umap', group.by = batch_label)
p22<- DimPlot(object = batches, reduction= 'umap', group.by = celltype_label)

png(paste0(plotout_dir,outfile_prefix,umapplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
print(plot_grid(p21, p22))
dev.off()
##############
seurat3_res <- as.data.frame(batches@reductions$pca@cell.embeddings)
cells_use <- rownames(seurat3_res)
seurat3_res$batch <- batches@meta.data[, batch_label]
seurat3_res$celltype <- batches@meta.data[, celltype_label]
write.csv(seurat3_res, file = paste0(saveout_dir,"dataset7—Seurat3_pca.csv"))
## We need new functions for calculate the coefficient beta
sce_obj <- SingleCellExperiment(assays=list(counts=expr_mat),colData=metadata)
cellind=colData(sce_obj)$celltype
batchind=colData(sce_obj)$batchlb
ssq <- function(x,batchlb,cell){
  mod <- lm(x~as.factor(batchlb)+as.factor(cell))
  av <- anova(mod)
  result <- av$`Sum Sq`
  return(result)
}
expvar=t(apply(assay(sce_obj,"counts"),1,ssq,batch=batchind,cell=cellind))
colnames(expvar)=c("batch","cell","residual")
ss=apply(expvar,2,sum) ## This is the sum of all SS across all the genes
r_raw_all <- ss[2]/ss[1] ## This is the overall ratio between cell and batch 
r_raw_all 
dim(assay(sce_obj, "counts"))
length(batchind)
length(cellind)
## Tables: cell type versus batch
table(colData(sce_obj)$batchlb, colData(sce_obj)$celltype)
table(colData(sce_obj)$CellType)
#########################
batches_data=batches@assays$integrated@data ###使用的是integrated的data而不是counts，counts矩阵为0
batches_data

# 现在，new_sce_obj 具有所需的行数和与原始数据相同的内容

cells_use <- colnames(batches_data)
colData<-sce_obj@colData
colData
# 创建一个新的colData对象，仅包含与样本数匹配的行
new_colData <- colData[cells_use, ]
new_colData
# 将新的colData对象分配给sce_obj
cellind=new_colData$celltype
batchind=new_colData$batchlb
batch_subset <- colData(sce_obj)$batchlb
sce_obj_corr<- SingleCellExperiment(assays=list(counts=batches_data),colData=new_colData)
expvar=t(apply(assay(sce_obj_corr,"counts"),1,ssq,batch=batchind,cell=cellind))
colnames(expvar)=c("batch","cell","residual")
ss=apply(expvar,2,sum) ## This is the sum of all SS across all the genes
ss
r_seurat3_all <- ss[2]/ss[1] ## This is the overall ratio between cell and batch 
r_seurat3_all

## Calculate the beta coefficients for seurat3

beta_seurat3_all <- r_seurat3_all/r_raw_all
beta_seurat3_all

