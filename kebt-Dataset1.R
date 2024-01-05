library(ggplot2)
library(kBET)
library(ggpubr)

rm(list=ls())
this_dir <-'/Users/fengjiao.li/Desktop/dataset1/'
setwd(this_dir)
eval_metric <- 'kBET/'
dataset_use <- 'dataset1'

utils_dir <- '/Users/fengjiao.li/Desktop/ComBat/'
source(paste0(utils_dir,'evaluation_utils.R'))
source(paste0(utils_dir,'kbet_utils.R'))

data_dir <- '/Users/fengjiao.li/Desktop/demo_ComBat/combat_dataset1/'

method_dir <- c('raw','combat','combatcell','combatsva','seurat3','combatseq','combatcellseq','combatsvaseq')
methods_use <- c('Raw','combat','combatcell','combatsva','Seurat_3','combatseq','combatcellseq','combatsvaseq')
# fn <- 'dataset2_resnet_pca.csv'
# myResnet <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
# head(myResnet)
# rownames(myResnet) <- gsub('-[0-9]$','',rownames(myResnet))
# write.csv(myResnet, file = paste0(data_dir, fn), row.names = T, quote = F)


dir.create(paste0(this_dir, eval_metric), showWarnings = F)
fn <- paste0(dataset_use, '_raw_pca.csv')
meta_ls <- get_celltype_common_kbet_V2(data_dir, fn)

size = length(meta_ls$cells_extract)
print('******Sample size is: ********')
print(size) 
pct = c(seq(from = 5, to = 25, by = 5))
kns = floor(pct*size/100)
print(kns)

res_dir = paste0(this_dir, eval_metric,'result/')
dir.create(res_dir, showWarnings = F)

for(i in 1:length(kns)){
  
  kn = kns[i]
  summary_KBET_v2(meta_ls, data_dir, methods_use, method_dir, dataset_use, eval_metric, this_dir, kn)
}

generate_median_kbet(pct, kns, res_dir, dataset_use, this_dir, kbet_output=T)
# 加载所需的库
install.packages("ggplot2")
# 输入数据
# 加载所需的库
library(ggplot2)
library(tidyr)

# 输入数据
data <- data.frame(
  Methods = c("ComBat", "ComBat-Cell", "ComBat-SVA", "Srurat3", "ComBat-Cell-Seq", "ComBat-SVA-Seq", "Before adjusted", "ComBat-Seq"),
  `5%` = c(0.632, 0.632, 0.395, 0.658, 0.237, 0.105, 0.105, 0.079),
  `10%` = c(0.632, 0.579, 0.316, 0.395, 0.211, 0.053, 0.079, 0),
  `15%` = c(0.658, 0.579, 0.342, 0.316, 0.158, 0.053, 0.053, 0),
  `20%` = c(0.684, 0.605, 0.342, 0.263, 0.237, 0.079, 0.053, 0),
  `25%` = c(0.724, 0.605, 0.395, 0.211, 0.276, 0.079, 0.053, 0.026)
)

names(data)[2] <- "5%"
names(data)[3] <- "10%"
names(data)[4] <- "15%"
names(data)[5] <- "20%"
names(data)[6] <- "25%"
# 数据整理为长格式
long_data <- gather(data, key = "Percentage", value = "Value", -Methods)
long_data$Percentage <- factor(long_data$Percentage, levels = c("5%", "10%", "15%", "20%", "25%"))

# 为Percentage列设置一个有序
# 使用已定义的数据和已经转换的长格式数据
# 制作折线图
ggplot(long_data, aes(x = Percentage, y = Value, color = Methods, group = Methods)) +
  geom_line() +
  geom_point() +
  labs(
       x = "% Sample Size",
       y = "kBET (acceptance rate)") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black")  # Add black axis lines
  )