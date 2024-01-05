library(devtools)
install_github('theislab/kBET')
install.packages("ggpubr")
library(ggplot2)
library(kBET)
library(ggpubr)

rm(list=ls())
this_dir <-'/Users/fengjiao.li/Desktop/demo_ComBat/'
setwd(this_dir)
eval_metric <- 'kBET/'
dataset_use <- 'dataset7'

utils_dir <- '/Users/fengjiao.li/Desktop/ComBat/'
source(paste0(utils_dir,'evaluation_utils.R'))
source(paste0(utils_dir,'kbet_utils.R'))

data_dir <- '/Users/fengjiao.li/Desktop/demo_ComBat/combat_dataset/'

method_dir <- c('raw')
methods_use <- c('Raw')
# fn <- 'dataset2_resnet_pca.csv'
# myResnet <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
# head(myResnet)
# rownames(myResnet) <- gsub('-[0-9]$','',rownames(myResnet))
# write.csv(myResnet, file = paste0(data_dir, fn), row.names = T, quote = F)


dir.create(paste0(this_dir, eval_metric), showWarnings = F)
fn <- paste0(dataset_use, '_raw_pca.csv')
meta_ls <- get_celltype_common_kbet_ds(data_dir, fn, percent_ext=0.2)


#### iteration 
size = length(meta_ls$cells_extract)  # 16600
print(paste0("Sample size", size))
pct = c(seq(from = 5, to = 25, by = 5))
kns = floor(pct*size/100)
print(kns)
res_dir <- paste0(this_dir, eval_metric,"result/")

install.packages("parallel")

library(parallel)
mclapply(1:length(kns), function(i){
  kn = kns[i]
  print(paste0('____k0 value use is: ', kn))
  summary_KBET_v2(meta_ls, data_dir, methods_use, method_dir, dataset_use, eval_metric, this_dir, kn)
  
}, mc.cores = 10)



generate_median_kbet(pct, kns, res_dir, dataset_use, this_dir, kbet_output=T)



############################
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
library(ggplot2)
library(tidyr)
library(scales)

# 输入数据
data <- data.frame(
  Methods = c("Seurat3", "ComBat-Cell", "Before adjusted", "ComBat-Cell-Seq", "ComBat- SVA-Seq", "ComBat-Seq", "ComBat-SVA", "ComBat"),
  `5%` = c(0.565, 0.5, 0.37, 0.196, 0.109, 0.109, 0.109, 0.109),
  `10%` = c(0.457, 0.522, 0.359, 0.359, 0.152, 0.13, 0.109, 0.065),
  `15%` = c(0.522, 0.522, 0.457, 0.413, 0.304, 0.283, 0.174, 0.152),
  `20%` = c(0.587, 0.435, 0.457, 0.391, 0.391, 0.391, 0.326, 0.283),
  `25%` = c(0.609, 0.435, 0.413, 0.337, 0.304, 0.304, 0.348, 0.293)
)
names(data)[2] <- "5%"
names(data)[3] <- "10%"
names(data)[4] <- "15%"
names(data)[5] <- "20%"
names(data)[6] <- "25%"

long_data <- gather(data, key = "Percentage", value = "Value", -Methods)
long_data$Percentage <- factor(long_data$Percentage, levels = c("5%", "10%", "15%", "20%", "25%"))
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