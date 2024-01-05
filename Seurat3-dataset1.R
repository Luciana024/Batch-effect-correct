#load packages
install.packages("devtools")
library(devtools)
BiocManager::install("Seurat")
remotes::install_version("Seurat", version = "3.0.1")
library(Seurat)
packageVersion('Seurat')
library(magrittr)
library(cowplot)
library(SingleCellExperiment)
rm(list=ls())
########################
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
outfile_prefix = "Dataset1"
save_obj = T
saveout_dir="/Users/fengjiao.li/Desktop/Seurat3/"
# load evaluation & preprocess utils files #
this.dir <- '/Users/fengjiao.li/Desktop/demo_Seurat3/'
setwd(this.dir)
read_dir <- '/Users/fengjiao.li/Desktop/dataset1/'
utils_dir <- '/Users/fengjiao.li/Desktop/Seurat3/'
source(paste0(utils_dir,'call_seurat_3.R'))
expr_filename = 'dataset1_sm_uc3.txt'
metadata_filename = 'sample_sm_uc3.txt'
batch_label = "batch"
celltype_label = "celltype"

########################
# read data 

expr<- read.table(file = paste0(read_dir,expr_filename),sep="\t",header=T,row.names=1,check.names = F)
metadata <- read.table(file = paste0(read_dir,metadata_filename),sep="\t",header=T,row.names=1,check.names = F)

expr<- expr[, rownames(metadata)]
expr_mat=as.matrix(expr)
class(expr_mat)
########################第二种尝试######
base_name <- 'Seurat3_dataset1/'
source(paste0(utils_dir,'ComBat_functions.R'))
###################我自己的过滤#########
myFilteredData<- filter_data_mtx(expr, base_name, is_filter_cells=TRUE,min_genes=300, 
                                  is_filter_genes=TRUE, min_cells=10)
lengths_list <- lengths(myFilteredData)
all_equal <- all(lengths_list == lengths_list[1])
matrix_data <- do.call(rbind, myFilteredData)
row_max <- matrixStats::rowMaxs(matrix_data)
row_min = matrixStats::rowMins(matrix_data)
if(sum(row_max==row_min)>0){
  myFilteredData = myFilteredData[row_max != row_min,]
}
cells_use <- colnames(myFilteredData)
############
metadata <- metadata[cells_use,]
seurat_obj <- CreateSeuratObject(counts =myFilteredData,meta.data = metadata)
expr_mat<-seurat_obj@assays$RNA@data
batch_list = seurat3_preprocess(
  expr_mat, metadata, 
  normData = normData, Datascaling = Datascaling, regressUMI = regressUMI, 
  min_cells = min_cells, min_genes = min_genes, 
  norm_method = norm_method, scale_factor = scale_factor, 
  numVG = numVG, nhvg = nhvg, 
  batch_label = batch_label, celltype_label = celltype_label)
######################################
# run pipeline
batches = call_seurat3(batch_list, batch_label, celltype_label, npcs, plotout_dir = this.dir, saveout_dir = this.dir, outfilename_prefix = outfile_prefix, visualize = visualize, save_obj = save_obj)
batches
########UMAP plot########
umapplot_filename <- "_UMAP_Plot"  

p21 <- DimPlot(object = batches, reduction= 'umap', group.by = batch_label)
p22 <- DimPlot(object = batches, reduction= 'umap', group.by = celltype_label)

png(paste0(plotout_dir,outfile_prefix,umapplot_filename,".png"),width = 2*1000, height = 800, res = 2*72)
print(plot_grid(p21, p22))
dev.off()
seurat3_res <- as.data.frame(batches@reductions$pca@cell.embeddings)
cells_use <- rownames(seurat3_res)
seurat3_res$batch <- batches@meta.data[, batch_label]
seurat3_res$celltype <- batches@meta.data[, celltype_label]
write.csv(seurat3_res, file = paste0(saveout_dir,"dataset1_seurat3_pca.csv"))

## We need new functions for calculate the coefficient beta
sce_obj <- SingleCellExperiment(assays=list(counts=expr_mat),colData=metadata)
cellind=colData(sce_obj)$celltype
batchind=colData(sce_obj)$batch
ssq <- function(x,batch,cell){
  means <- tapply(x,list(cell,batch),mean)
  gmean <- matrix(rep(tapply(x,cell,mean),times=2),ncol=2,byrow=FALSE)
  sizes <- table(cell,batch)
  result <- (means-gmean)^2*sizes
  result <- na.omit(result)
  return(sum(result))
}
ssq1<- function(x,batch,cell){
  mod <- lm(x~as.factor(batch)+as.factor(cell))
  av <- anova(mod)
  result <- av$`Sum Sq`
  return(result)
}
expvar=t(apply(assay(sce_obj,"counts"),1,ssq1,batch=batchind,cell=cellind))
colnames(expvar)=c("batch","cell","residual")
ss=apply(expvar,2,sum) ## This is the sum of all SS across all the genes
r_raw_all <- ss[2]/ss[1] ## This is the overall ratio between cell and batch 
r_raw_all 
ss[2]
ss[1]
expvar_raw=apply(assay(sce_obj,"counts"),1,ssq,batch=batchind,cell=cellind)
## expvar is a vector of SSQ between 2 batches for all different genes conditional on the cell types
## We can sum it to get an overall estimate across genes
sum(expvar_raw)
## Tables: cell type versus batch
table(colData(sce_obj)$batch, colData(sce_obj)$celltype)
table(colData(sce_obj)$celltype)
#########################
batches_data=batches@assays$integrated@data ###使用的是integrated的data而不是counts，counts矩阵为0
sample_subset <- colData(sce_obj)
batch_subset <- colData(sce_obj)$batch
rownames(batches_data)
rownames(sample_subset)
sce_obj_corr<- SingleCellExperiment(assays=list(counts=batches_data),colData=sample_subset)
## Calculate the ratio between cell and batch
#######第一种指标######
expvar=t(apply(assay(sce_obj_corr,"counts"),1,ssq1,batch=batchind,cell=cellind))
colnames(expvar)=c("batch","cell","residual")
ss=apply(expvar,2,sum) ## This is the sum of all SS across all the genes
ss
r_seurat3_all <- ss[2]/ss[1] ## This is the overall ratio between cell and batch 
r_seurat3_all
ss[2]
ss[1]
## Calculate the beta coefficients for seurat3
beta_seurat3_all <- r_seurat3_all/r_raw_all
beta_seurat3_all
