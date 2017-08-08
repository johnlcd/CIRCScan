## Feature selection


# library packages
require(mlbench)
require(caret)
require(DT)
library(doParallel)
library(PerfMeas)
library(gplots)
library(ROCR)

# set arguments
args <- commandArgs(T)
cell <- args[1]
PM <- args[2]
if (args[4] != 'all'){
	FN_list <- as.numeric(strsplit(args[4], ',')[[1]])
} else{
	FN_list <- seq(1,FN)
}
print('>>> Feature number list:')
print(FN_list)

## load data
load(paste(cell, PM, 'train.RData', sep = '_'))

# set current time
t1 = Sys.time()

# show objects and which libraries were loaded
print('>>> Objects are:')
ls()
print('>>> Session info:')
sessionInfo()

# register parallel front-end
cores <- as.numeric(args[3])
cl <- makeCluster(cores); registerDoParallel(cl)

# model traning with different feature number
ctrl_cv <- trainControl(method = 'cv', number = 10, allowParallel = T)

F1_best <- 0.5
gp_best <- 0
for (gp in 1:5)
{

	print('############################################')
	print(paste('>>> Group', gp, ':'))
	train <- get(paste('train', gp, sep = ''))
	test <- get(paste('test', gp, sep = ''))
	test$Type <- as.factor(test$Type)
	X = train[,5:(FN+4)]
	str(X)
	Y = as.factor(train$Type)

# select features by importance
#	sort_imp <- get(paste('sort_imp_gp', gp, sep = ''))
	assign(paste('sort_fea_gp', gp, '_all', sep = ''), c(rownames(get(paste('sort_imp_gp', gp, sep = '')))))
	assign(paste('old_sort_fea_gp', gp, sep = ''), c(rownames(get(paste('sort_imp_gp', gp, sep = '')))))
	
	assign(paste('F1_best_gp', gp, sep = ''), 0.5)
	for (i in 1:length(FN_list))
	{
		print('=====================================================')
		n <- rev(sort(FN_list))[i]
		print('>>> Selcted top feature number:')
		print(n)
		assign(paste('tmp_sort_fea_gp', gp, sep = ''), get(paste('old_sort_fea_gp', gp, sep = ''))[1:n])
		print(paste('>>> Top', n, 'Feature(s) is(are):', sep = ' '))
		print(get(paste('tmp_sort_fea_gp', gp, sep = '')))
		assign(paste('Model_top', n, sep = '_'), train(y = Y, x = train[get(paste('tmp_sort_fea_gp', gp, sep = ''))], method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale")))
		Model_tmp <- get(paste('Model_top', n, sep = '_'))
		test_mat_tmp <- test[get(paste('tmp_sort_fea_gp', gp, sep = ''))]
		Pred_tmp <- predict(Model_tmp, test_mat_tmp)
		Prob_tmp <- predict(Model_tmp, test_mat_tmp, type = 'prob')
		Pred_prob_tmp <- Prob_tmp[, 2]
		print('>>> Length of predicted prob:')
		print(length(Pred_prob_tmp))
		print('>>> Length of labels:')
		print(length(test$Type))
		Prediction_tmp <- prediction(predictions = Pred_prob_tmp, labels = test$Type)
		Perf.roc_tmp <- performance(Prediction_tmp, measure = 'tpr', x.measure = 'fpr')
		Perf.auc_tmp <- performance(Prediction_tmp, measure = 'auc')
		print(confusionMatrix(Pred_tmp, test$Type, positive = 'TRUE'))
		
		Precision_tmp <- posPredValue(Pred_tmp, test$Type, positive = 'TRUE')
		Recall_tmp <- sensitivity(Pred_tmp, test$Type, positive = 'TRUE') # Sensitivity
		Spe_tmp <- specificity(Pred_tmp, test$Type, negative = 'FALSE')
		AUC_tmp <- unlist(Perf.auc_tmp@y.values)
		F1_tmp <- (2 * Precision_tmp * Recall_tmp) / (Precision_tmp + Recall_tmp) 
		print('>>> Precision is:')
		print(Precision_tmp)
		print('>>> Sensitivity is:')
		print(Recall_tmp)
		print('>>> Specificity is:')
		print(Spe_tmp)
		print('>>> AUC is:')
		print(AUC_tmp)
		print('>>> F1 score is:')
		print(F1_tmp)

		if (F1_tmp > get(paste('F1_best_gp', gp, sep = ''))){
			assign(paste('Pred_best_gp', gp, sep = ''), Pred_tmp)
			assign(paste('Prob_best_gp', gp, sep = ''), Prob_tmp)
			assign(paste('Pred_prob_best_gp', gp, sep = ''), Pred_prob_tmp)
			assign(paste('Prediction_best_gp', gp, sep = ''), Prediction_tmp)
			assign(paste('Perf.roc_best_gp', gp, sep = ''), Perf.roc_tmp)
			assign(paste('Perf.auc_best_gp', gp, sep = ''), Perf.auc_tmp)
			assign(paste('Model_best_gp', gp, sep = ''), Model_tmp)
			assign(paste('Precision_best_gp', gp, sep = ''), Precision_tmp)
			assign(paste('Recall_best_gp', gp, sep = ''), Recall_tmp)
			assign(paste('Spe_best_gp', gp, sep = ''), Spe_tmp)
			assign(paste('AUC_best_gp', gp, sep = ''), AUC_tmp)
			assign(paste('F1_best_gp', gp, sep = ''), F1_tmp)
			assign(paste('fn_best_gp', gp, sep = ''), n)
			assign(paste('fea_best_gp', gp, sep = ''), get(paste('tmp_sort_fea_gp', gp, sep = '')))
		}

# estimate new importance
		assign(paste('importance_gp', gp, '_fea', n, sep = ''), varImp(Model_tmp, scale=FALSE))
		assign(paste('sort_imp_gp', j, '_fea', n, sep = ''), as.data.frame(sortImp(get(paste('importance_gp', gp, '_fea', n, sep = '')), n)))
		print('>>> New ranking of importance:')
		assign(paste('old_sort_fea_gp', gp, sep = ''), c(rownames(get(paste('sort_imp_gp', j, '_fea', n, sep = '')))))
		print(get(paste('old_sort_fea_gp', gp, sep = '')))

# running time
		print(difftime(Sys.time(), t1, units = 'sec'))
		print('=====================================================')
	}

	print('>>> Best model (feature number):')
	print(paste('Top', get(paste('fn_best_gp', gp, sep = ''))))
	print('>>> Best features are:')
	print(get(paste('fea_best_gp', gp, sep = '')))
	print('>>> Confusion matrix:')
	print(confusionMatrix(get(paste('Pred_best_gp', gp, sep = '')), test$Type, positive = 'TRUE'))
	print('>>> Precision is:')
	print(get(paste('Precision_best_gp', gp, sep = '')))
	print('>>> Sensitivity is:')
	print(get(paste('Recall_best_gp', gp, sep = '')))
	print('>>> Specificity is:')
	print(get(paste('Spe_best_gp', gp, sep = '')))
	print('>>> AUC is:')
	print(get(paste('AUC_best_gp', gp, sep = '')))
	print('>>> F1 score is:')
	print(get(paste('F1_best_gp', gp, sep = '')))
	print(paste('>>> Feature selection of group', gp, 'finished.'))
	print('############################################')
	
	if (get(paste('F1_best_gp', gp, sep = '')) > F1_best){
		gp_best <- gp
		Pred_best <- get(paste('Pred_best_gp', gp, sep = ''))
		Prob_best <- get(paste('Prob_best_gp', gp, sep = ''))
		Pred_prob_best <- get(paste('Pred_prob_best_gp', gp, sep = ''))
		Prediction_best <- get(paste('Prediction_best_gp', gp, sep = ''))
		Perf.roc_best <- get(paste('Perf.roc_best_gp', gp, sep = ''))
		Perf.auc_best <- get(paste('Perf.auc_best_gp', gp, sep = ''))
		Model_best <- get(paste('Model_best_gp', gp, sep = ''))
		Precision_best <- get(paste('Precision_best_gp', gp, sep = ''))
		Recall_best <- get(paste('Recall_best_gp', gp, sep = ''))
		Spe_best <- get(paste('Spe_best_gp', gp, sep = ''))
		AUC_best <- get(paste('AUC_best_gp', gp, sep = ''))
		F1_best <- get(paste('F1_best_gp', gp, sep = ''))
		fn_best <- get(paste('fn_best_gp', gp, sep = ''))
		fea_best <- get(paste('fea_best_gp', gp, sep = ''))
	}

}

# best model
pdf(file = paste(cell, PM, 'ROC_curve.pdf', sep = '_'))

plot(Perf.roc_best, main ='ROC curve', col = 'blue', lwd = 3)
abline(a = 0, b = 1, col = 'red', lwd = 2.5, lty = 2)

dev.off()

print('+++++++++++++++++++++++++++++++++++++++++++++++++++++')
if (exists('Model_best')){
	print('>>> Best group:')
	print(paste('Group', gp_best))
	print('>>> Best model:')
	print(paste('Top', fn_best))
	print('>>> Best features are:')
	print(fea_best)
#	print('>>> Confusion matrix:')
#	print(confusionMatrix(Pred_best, test$Type, positive = 'TRUE'))
	print('>>> Precision is:')
	print(Precision_best)
	print('>>> Sensitivity is:')
	print(Recall_best)
	print('>>> Specificity is:')
	print(Spe_best)
	print('>>> AUC is:')
	print(AUC_best)
	print('>>> F1 score is:')
	print(F1_best)
	print('#END')
} else{
	print('>>> No best model (F1 score greater than 0.5).')
}
print('+++++++++++++++++++++++++++++++++++++++++++++++++++++')


# stop cluster and register sequntial front end
stopCluster(cl); registerDoSEQ();

# save R data
save(list = objects(), file = paste(cell, PM, 'FS.RData', sep = '_'))


### END
