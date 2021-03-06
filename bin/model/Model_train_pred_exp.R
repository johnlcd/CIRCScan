# please make sure required packages already installed


# library packages
require(mlbench)
require(caret)
require(DT)
library(doParallel)
library(PerfMeas)
library(gplots)
library(ROCR)
library(randomForest)
library(Metrics)

# set arguments
args2 <- commandArgs(T)
cell <- args2[1]
PM <- args2[2]

# load data of feature selection
load(paste(cell, PM, 'pred_exp_all.RData', sep = '_'))

cores <- as.numeric(args2[3])
#seed <- 111

# load data_train
#file <- paste(cell, 'exp_train', sep = '_')
#data_train_all <- read.table(file, head = T)
#data_train_all <- data.frame(data_train_all)
#data_train_raw <- data_train_all
#data_train_all_mat <- data_train_all[,5:(FN+4)]
#data_train_all$SRPBM <- log2(data_train_all$SRPBM)
#data_train <- data_train_all
#data_train_mat <- data_train[,5:(FN+4)]
#all_num <- dim(data_train)[1]
#if (PM == 'glm') {
#	data_train$SRPBM <- data_train$SRPBM/10 # value = log2(SRPBM)/10
#}
#summary(data_train)
cat('>>> Dimension of data matrix: \n')
dim(data_train)
n <- dim(data_train)[1]
FN <- dim(data_train)[2] - 5
col_name <- names(data_train)
fea_all <- col_name[5:(FN+4)]
fea_sel <- fea_all
data_train_mat <- data_train[,5:(FN+4)]


# show which libraries were loaded  
cat('>>> Session Info: \n')
sessionInfo()

# register parallel front-end
#cl <- makeCluster(cores); registerDoParallel(cl)

cat('>>> [1] Summary of trained model for circRNAs expression (SRPBM) with all training data ... ... \n')
#ctrl <- trainControl(method = 'cv', number = 10 , savePredictions = "all", returnData = T, returnResamp = 'all', 
#					 verboseIter = T, allowParallel = T)
#Model_all <- train(y = data_train$SRPBM, 
#			   x = data_train_mat, 
#			   method = PM, trControl = ctrl, tuneGrid = tune_best, metric = "RMSE", importance = T, preProc = c("center", "scale"))

# model summary
#cat('    Summary of model: \n')
print(Model_all)

# stop cluster and register sepuntial front end
#stopCluster(cl); registerDoSEQ()

# load pred data
#load(paste('CELL', 'train_pred_exp.RData', sep = '_'))
cat('>>> [2] Predict circRNAs expression (SRPBM) with final model ... ... \n')
data_pred <- read.table(paste(cell, 'exp_pred', sep = '_'), head = T)
data_pred.bed <- data_pred[1:4]
all_pred_mat <- data_pred[-(1:3)]
#sel_pred_mat <- all_pred_mat[,c('Alu','H3K36me3','H3K79me2')]
#sel_pred_mat <- all_pred_mat[,c('Alu','DNaseI_HS','H3K27ac','H3K36me3','H3K4me1','H3K4me2','H3K4me3','H3K79me2')]
sel_pred_mat <- all_pred_mat[,fea_sel]
rownames(sel_pred_mat) <- data_pred$Intron_pair
cat('>>> Number of intron pairs and all feature: \n')
cat(paste(dim(all_pred_mat)[1], dim(all_pred_mat)[2]-2, "\n", sep = ';'))

# prediction by trained model
circ_pred <- predict(Model_all, sel_pred_mat)
circ_pred <- 2^circ_pred
circ_pred_df <- data.frame(circ_pred)
circ_pred_bed <- data.frame(cbind(data_pred.bed, circ_pred_df))
colnames(circ_pred_bed) <- c('Chr','Start','End','Intron_pair','SRPBM')
#circ_pred_bed$SRPBM <- 2^(circ_pred_bed$SRPBM)
#circ_pred_all <- predict(Model_all_all, all_pred_mat[fea_all])
#circ_pred_all <- 2^circ_pred_all
#circ_pred_all_df <- data.frame(circ_pred_all)
#circ_pred_all_bed <- data.frame(cbind(data_pred.bed, circ_pred_all_df))
#colnames(circ_pred_all_bed) <- c('Chr','Start','End','Intron_pair','SRPBM')
max_exp <- max(circ_pred_df)
min_exp <- min(circ_pred_df)
max_exp_bed <- circ_pred_bed[c(which(circ_pred_df == max_exp)),]
min_exp_bed <- circ_pred_bed[c(which(circ_pred_df == min_exp)),]
#max_exp_df <- data.frame(cbind(max_exp_bed,data.frame(circ_pred[c(which(circ_pred_df == max_exp))])))
#colnames(max_exp_df) <- c('Chr','Start','End','Intron_pair','SRPBM')
#min_exp_df <- data.frame(cbind(min_exp_bed,data.frame(circ_pred[c(which(circ_pred_df == min_exp))])))
#colnames(min_exp_df) <- c('Chr','Start','End','Intron_pair','SRPBM')
cat('    The max expressed circRNA: \n')
print(max_exp_bed)
cat('    The min expressed circRNA: \n')
print(min_exp_bed)

write.table(circ_pred_bed, paste(cell, PM, 'pred_exp.bed', sep = '_'), row.names = F, col.names = F, quote = F, sep = '\t')
#write.table(cbind(data_pred.bed[,1:4], circ_pred_all_df), paste(cell, PM, 'pred_exp_all.bed', sep = '_'), row.names = F, col.names = F, quote = F, append = T, sep = '\t')


# save R data
rm(list = 'args')
save(list = objects(), file=paste(cell, PM, 'train_pred_exp.RData', sep = '_'))


### END
