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
cat('>>> Objects: \n')
objects()

# register parallel front-end
cores <- as.numeric(args[3])
cl <- makeCluster(cores); registerDoParallel(cl)

# re-train model with all data
cat('>>> Re-train model with best features: \n')

cat('+++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

cat('>>> Best features are: \n')
print(fea_best)
cat('>>> Feature number: \n')
print(length(fea_best))
data_train_mat <- data.frame(data_train[fea_best])
data_train$Type <- as.factor(data_train$Type)
ctrl <- trainControl(method = 'cv', number = 10, savePredictions = "all", returnData = T, returnResamp = 'all', 
					 verboseIter = T, allowParallel = T)
Model_final <- train(y = data_train$Type, x = data_train_mat, method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale"))
cat('>>> Summary of final model: \n')
summary(Model_final)

Pred <- predict(Model_final, data_train_mat)
if (Pred_type == 'prob') {
	Prob <- predict(Model_final, data_train_mat, type = 'prob')
	Pred_prob <- Prob[, 2]
	Prediction <- prediction(predictions = Pred_prob, labels = data_train$Type)
	Perf.roc <- performance(Prediction, measure = 'tpr', x.measure = 'fpr')
	Perf.auc <- performance(Prediction, measure = 'auc')
}
results <- confusionMatrix(Pred, data_train$Type, positive = 'TRUE')
cat('>>> Confusion matrix: \n')
print(results)
ACC <- as.numeric(results$byClass['Balanced Accuracy'])
Precision <- as.numeric(results$byClass['Precision'])
Recall <- as.numeric(results$byClass['Recall']) # Sensitivity
Spe <- as.numeric(results$byClass['Specificity'])
if (Pred_type == 'prob') {
	AUC <- unlist(Perf.auc@y.values)
}
F1 <- as.numeric(results$byClass['F1'])
cat('>>> Precision is: \n')
print(Recall)
cat('>>> Specificity is: \n')
print(Spe)
cat('>>> ACC is: \n')
print(ACC)
if (Pred_type == 'prob') {
	cat('>>> AUC is:\n')
	print(AUC)
}
cat('>>> F1 score is: \n')
print(F1)

cat('+++++++++++++++++++++++++++++++++++++++++++++++++++++\n')


# stop cluster and register sequntial front end
stopCluster(cl); registerDoSEQ();

# load pred data
cat('>>> Predict circRNAs with final model ... ... \n')

data_pred <- read.table(paste(cell, 'pred', sep = '_'), head = T)
data_pred.bed <- data_pred[1:4]
all_pred_mat <- data_pred[-(1:3)]
sel_pred_mat <- all_pred_mat[fea_best]
rownames(sel_pred_mat) <- data_pred$Intron_pair
cat('>>> Number of intron pairs and all feature: \n')
cat(paste(dim(all_pred_mat)[1], dim(all_pred_mat)[2]-2, "\n", sep = ';'))
cat('>>> Selected features: \n')
print(fea_best)
cat('>>> Selected feature number: \n')
print(length(fea_best))
cat('>>> Head of prediction matrix: \n')
head(sel_pred_mat)

# get threshold
cat('>>> [1] Get the threshold -->> \n')
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
cat('    Threshold for prediction: \n')
print(threshold)

# prediction by trained model
cat('>>> [2] Prediction and classify by threshold -->> \n')
circ_pred_prob <- predict(Model_final, sel_pred_mat, type = 'prob')
true_prob <- circ_pred_prob[,2]
true_num <- length(which(true_prob >= threshold))
cat('    Number of predicted circRNAs: \n')
print(true_num)

true_row_num <- which(true_prob >= threshold)
write.table(data_pred.bed[true_row_num,], paste(cell, PM, 'pred_true.bed', sep = '_'), row.names = F, col.names = F, quote = F, append = TRUE, sep = '\t')


# summary of prediction
cat('>>> Sunmmary of prediction: \n')
summary(circ_pred_prob)

# write into file and save R data
obj_save <- c('PM', 'cell', 'fea_best', 'Model_final', 'data_train', 'data_train_mat', 'ctrl', 'data_pred', 'data_pred.bed', 'sel_pred_mat', 'roc_num', 'youden', 'threshold', 'circ_pred_prob', 'true_prob', 'true_num', 'true_row_num')
save(list = obj_save, file = paste(cell, PM, 'pred.RData', sep = '_'))


### END
