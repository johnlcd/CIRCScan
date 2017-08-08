# Used for prediction


# library packages
require(mlbench)
require(caret)
require(DT)
library(doParallel)
library(doMC)
library(PerfMeas)
library(gplots)
library(ROCR)

# set arguments
args <- commandArgs(T)
cell <- args[1]
PM <- args[2]
#FN_sel <- as.numeric(args[3])

# load trained models
load(paste(cell, PM, 'FS.RData', sep = '_'))

# show R objects
print('>>> Objects:')
objects()

# register parallel front-end
cores <- as.numeric(args[3])
cl <- makeCluster(cores); registerDoParallel(cl)

# re-train model with all data
data_train_mat <- data_train[fea_best]
data_train$Type <- as.factor(data_train$Type)
ctrl <- trainControl(method = 'cv', number = 10, allowParallel = T)
Model_final <- train(y = data_train$Type, x = data_train_mat, method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale"))
print('>>> Model summary:')
summary(Model_final)

# stop cluster and register sequntial front end
stopCluster(cl); registerDoSEQ();

# load pred data
data_pred <- read.table(paste(cell, 'pred', sep = '_'), head = T)
data_pred.bed <- data_pred[1:4]
all_pred_mat <- data_pred[-(1:3)]
sel_pred_mat <- all_pred_mat[fea_best]
rownames(sel_pred_mat) <- data_pred$Intron_pair
print('>>> Number of intron pairs and all feature:')
print(paste(dim(all_pred_mat)[1], dim(all_pred_mat)[2]-2, sep = ';'))
print('>>> Selected features:')
print(fea_best)
print('>>> Selected feature number:')
print(length(fea_best))
print('>>> Head of prediction matrix:')
head(sel_pred_mat)
#model_pred <- Model_best

# get threshold
roc_num <- length(Perf.roc_best@x.values[[1]])
youden <- list()
i = 0
while (i < roc_num){
i = i + 1
youden_i <- Perf.roc_best@y.values[[1]][i]-Perf.roc_best@x.values[[1]][i]
youden <- c(youden,youden_i)
}
youden_all <- do.call(rbind, youden)
line <- which(youden_all==youden_all[which.max(youden_all)], arr.ind = T)[1,1]
threshold <- Perf.roc_best@alpha.values[[1]][line]
print('>>> Threshold for prediction:')
print(threshold)

# prediction by trained model
circ_pred_prob <- predict(Model_final, sel_pred_mat, type = 'prob')
true_prob <- circ_pred_prob[,2]
true_num <- length(which(true_prob >= threshold))
print('>>> Number of predicted circRNAs:')
print(true_num)

true_row_num <- which(true_prob >= threshold)
write.table(data_pred.bed[true_row_num,], paste(cell, PM, 'pred_true.bed', sep = '_'), row.names = F, col.names = F, quote = F, append = TRUE, sep = '\t')


# summary of prediction
print('>>> Sunmmary of prediction:')
summary(circ_pred_prob)

# write into file and save R data
save(list = objects(), file = paste(cell, PM, 'pred.RData', sep = '_'))


### END
