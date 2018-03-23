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
cores <- as.numeric(args[3])

# load training data
file <- paste(cell, 'train', sep = '_')
data_train <- read.table(file, head = T)
summary(data_train)
fea <- c('Alu', 'H3K36me3', 'H3K79me2')
data_train_mat <- data_train[fea]
data_train$Type <- as.factor(data_train$Type)

# register parallel front-end
cl <- makeCluster(cores); registerDoParallel(cl)

# model training
ctrl <- trainControl(method = 'cv', number = 10, allowParallel = T)
Model <- train(y = data_train$Type, x = data_train_mat, method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale"))
print('>>> Model summary:')
summary(Model)

Pred <- predict(Model, data_train_mat)
Prob <- predict(Model, data_train_mat, type = 'prob')
Pred_prob <- Prob[, 2]
Prediction <- prediction(predictions = Pred_prob, labels = data_train$Type)
Perf.roc <- performance(Prediction, measure = 'tpr', x.measure = 'fpr')
Perf.auc <- performance(Prediction, measure = 'auc')
results <- confusionMatrix(Pred, data_train$Type, positive = 'TRUE')
print('>>> Confusion matrix:')
print(results)
ACC <- as.numeric(results$byClass['Balanced Accuracy'])
Precision <- as.numeric(results$byClass['Precision'])
Recall <- as.numeric(results$byClass['Recall']) # Sensitivity
Spe <- as.numeric(results$byClass['Specificity'])
AUC <- unlist(Perf.auc@y.values)
F1 <- as.numeric(results$byClass['F1'])
print('>>> Precision is:')
print(Precision)
print('>>> Sensitivity is:')
print(Recall)
print('>>> Specificity is:')
print(Spe)
print('>>> ACC is:')
print(ACC)
print('>>> AUC is:')
print(AUC)
print('>>> F1 score is:')
print(F1)


# stop cluster and register sequntial front end
stopCluster(cl); registerDoSEQ();

# load pred data
data_pred <- read.table(paste(cell, 'pred', sep = '_'), head = T)
data_pred.bed <- data_pred[1:4]
all_pred_mat <- data_pred[-(1:3)]
sel_pred_mat <- all_pred_mat[fea]
rownames(sel_pred_mat) <- data_pred$Intron_pair
print('>>> Number of intron pairs to predict and features:')
print(paste(dim(all_pred_mat)[1], length(fea), sep = ';'))
print('>>> Features:')
print(fea)
print('>>> Head of prediction matrix:')
head(sel_pred_mat)

# get threshold
roc_num <- length(Perf.roc@x.values[[1]])
youden <- list()
i = 0
while (i < roc_num){
i = i + 1
youden_i <- Perf.roc@y.values[[1]][i]-Perf.roc@x.values[[1]][i]
youden <- c(youden,youden_i)
}
youden_all <- do.call(rbind, youden)
line <- which(youden_all==youden_all[which.max(youden_all)], arr.ind = T)[1,1]
threshold <- Perf.roc@alpha.values[[1]][line]
print('>>> Threshold for prediction:')
print(threshold)

# prediction by trained model
circ_pred_prob <- predict(Model, sel_pred_mat, type = 'prob')
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
save(list = objects(), file = paste(cell, PM, 'fast_model.RData', sep = '_'))


### END
