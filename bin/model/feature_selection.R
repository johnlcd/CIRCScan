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

## load data
load(paste(cell, PM, 'train.RData', sep = '_'))
FN_list <- c(as.numeric(seq(1,FN)))

args <- commandArgs(T)
if (args[4] != 'all'){
	FN_list <- as.numeric(strsplit(args[4], ',')[[1]])
}
# sort and reverse FN_list (decrease)
FN_list <- rev(sort(FN_list))
print('>>> Feature number list:')
print(FN_list) 
len_FN_list <- length(FN_list)
max_FN <- FN_list[1]
max_fea <- sort_fea_final[1:max_FN]

ref_index_raw <- args[5]
if (ref_index_raw == '--auc') {
   ref_index <<- 'AUC'
}
if (ref_index_raw == '--f1') {
   ref_index <<- 'F1'
} 

print('>>> The reference index is: ')
print(ref_index)
Pred_type <- args[6]

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

print('############################################')
print(paste('>>> Top', max_FN, 'Features are:', sep = ' '))
print(max_fea)

for (gp in 1:5)
{
	print('=====================================================')
	print(paste('>>> Group', gp, ':'))
	train <- get(paste('train', gp, sep = ''))
	test <- get(paste('test', gp, sep = ''))
	test$Type <- factor(as.character(test$Type))
	print('>>> Structrue and dimension of training matrix: ')
	X <- data.frame(train[max_fea])
	str(X)
	print(dim(X))
	Y <- factor(as.character(train$Type))
	print('>>> Lenght of training set labels: ')
	print(length(Y))

	assign(paste('Model_fea', max_FN, 'gp', gp, sep = '_'), train(y = Y, x = X, method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale")))
	Model_tmp <- get(paste('Model_fea', max_FN, 'gp', gp, sep = '_'))
	test_mat_tmp <- test[max_fea]
	Pred_tmp <- predict(Model_tmp, test_mat_tmp)
	assign(paste('Pred_fea', max_FN, 'gp', gp, sep = '_'), Pred_tmp)
	if (Pred_type == 'prob') {
		Prob_tmp <- predict(Model_tmp, test_mat_tmp, type = 'prob')
		Pred_prob_tmp <- Prob_tmp[, 2]
		print('>>> Length of predicted prob:')
		print(length(Pred_prob_tmp))
		print('>>> Length of testing set labels:')
		print(length(test$Type))
		Prediction_tmp <- prediction(predictions = Pred_prob_tmp, labels = test$Type)
		Perf.roc_tmp <- performance(Prediction_tmp, measure = 'tpr', x.measure = 'fpr')
		Perf.auc_tmp <- performance(Prediction_tmp, measure = 'auc')
	}
# confusion matrix
	results_tmp <- confusionMatrix(Pred_tmp, test$Type, positive = 'TRUE')
	print('>>> Confusion matrix: ')
	print(results_tmp)

	ACC_tmp <- as.numeric(results_tmp$byClass['Balanced Accuracy'])
	Precision_tmp <- as.numeric(results_tmp$byClass['Precision'])
	Recall_tmp <- as.numeric(results_tmp$byClass['Recall']) # Sensitivity
	Spe_tmp <- as.numeric(results_tmp$byClass['Specificity'])
	if (Pred_type == 'prob') {
		AUC_tmp <- unlist(Perf.auc_tmp@y.values)
	}
	F1_tmp <- as.numeric(results_tmp$byClass['F1'])
	print('>>> Precision is:')
	print(Precision_tmp)
	print('>>> Sensitivity is:')
	print(Recall_tmp)
	print('>>> Specificity is:')
	print(Spe_tmp)
	print('>>> ACC is:')
	print(ACC_tmp)
	if (Pred_type == 'prob') {
		print('>>> AUC is:')
		print(AUC_tmp)
	}
	print('>>> F1 score is:')
	print(F1_tmp)

	if (Pred_type == 'prob') {
		if (gp == 1) {
			AUC_list_tmp <- c(AUC_tmp)
		} else {
			AUC_list_tmp <- c(AUC_list_tmp, AUC_tmp)
		}
	}
	if (gp == 1) {
		F1_list_tmp <- c(F1_tmp)
		ACC_list_tmp <- c(ACC_tmp)
		Precision_list_tmp <- c(Precision_tmp)
		Recall_list_tmp <- Recall_tmp
		Spe_list_tmp <- Spe_tmp
	} else {
		F1_list_tmp <- c(F1_list_tmp, F1_tmp)
		ACC_list_tmp <- c(ACC_list_tmp, ACC_tmp)
		Precision_list_tmp <- c(Precision_list_tmp, Precision_tmp)
		Recall_list_tmp <- c(Recall_list_tmp, Recall_tmp)
		Spe_list_tmp <- c(Spe_list_tmp, Spe_tmp)
	}
	
# running time
	print(difftime(Sys.time(), t1, units = 'sec'))
	print('=====================================================')

	if (Pred_type == 'prob') {
		assign(paste('AUC_list_fea', max_FN, sep = ''), AUC_list_tmp)
	}
	assign(paste('F1_list_fea', max_FN, sep = ''), F1_list_tmp)
	assign(paste('ACC_list_fea', max_FN, sep = ''), ACC_list_tmp)
	assign(paste('Precision_list_fea', max_FN, sep = ''), Precision_list_tmp)
	assign(paste('Recall_list_fea', max_FN, sep = ''), Recall_list_tmp)
	assign(paste('Spe_list_fea', max_FN, sep = ''), Spe_list_tmp)

}


if (Pred_type == 'prob') {
	assign(paste('AUC_sum_fea', max_FN, sep = ''), sum(get(paste('AUC_list_fea', max_FN, sep = ''))))
	assign(paste('AUC_mean_fea', max_FN, sep = ''), mean(get(paste('AUC_list_fea', max_FN, sep = ''))))
}
assign(paste('F1_sum_fea', max_FN, sep = ''), sum(get(paste('F1_list_fea', max_FN, sep = ''))))
assign(paste('F1_mean_fea', max_FN, sep = ''), mean(get(paste('F1_list_fea', max_FN, sep = ''))))
assign(paste('ACC_sum_fea', max_FN, sep = ''), sum(get(paste('ACC_list_fea', max_FN, sep = ''))))
assign(paste('ACC_mean_fea', max_FN, sep = ''), mean(get(paste('ACC_list_fea', max_FN, sep = ''))))
assign(paste('Precision_sum_fea', max_FN, sep = ''), sum(get(paste('Precision_list_fea', max_FN, sep = ''))))
assign(paste('Precision_mean_fea', max_FN, sep = ''), mean(get(paste('Precision_list_fea', max_FN, sep = ''))))
assign(paste('Recall_sum_fea', max_FN, sep = ''), sum(get(paste('Recall_list_fea', max_FN, sep = ''))))
assign(paste('Recall_mean_fea', max_FN, sep = ''), mean(get(paste('Recall_list_fea', max_FN, sep = ''))))
assign(paste('Spe_sum_fea', max_FN, sep = ''), sum(get(paste('Spe_list_fea', max_FN, sep = ''))))
assign(paste('Spe_mean_fea', max_FN, sep = ''), mean(get(paste('Spe_list_fea', max_FN, sep = ''))))

if (Pred_type == 'prob') {
	print('>>> AUC of 5 groups: ')
	print(get(paste('AUC_list_fea', max_FN, sep = '')))
	print('>>> Mean value of AUC: ')
	print(get(paste('AUC_mean_fea', max_FN, sep = '')))
}
print('>>> F1 of 5 groups: ')
print(get(paste('F1_list_fea', max_FN, sep = '')))
print('>>> Mean value of F1: ')
print(get(paste('F1_mean_fea', max_FN, sep = '')))
print('>>> ACC of 5 groups: ')
print(get(paste('ACC_list_fea', max_FN, sep = '')))
print('>>> Mean value of ACC: ')
print(get(paste('ACC_mean_fea', max_FN, sep = '')))
print('>>> Precision of 5 groups: ')
print(get(paste('Precision_list_fea', max_FN, sep = '')))
print('>>> Mean value of precision: ')
print(get(paste('Precision_mean_fea', max_FN, sep = '')))
print('>>> Sensitivity (Recall) of 5 groups: ')
print(get(paste('Recall_list_fea', max_FN, sep = '')))
print('>>> Mean value of Sensitivity (Recall): ')
print(get(paste('Recall_mean_fea', max_FN, sep = '')))
print('>>> Specificity of 5 groups: ')
print(get(paste('Spe_list_fea', max_FN, sep = '')))
print('>>> Mean value of Specificity: ')
print(get(paste('Spe_mean_fea', max_FN, sep = '')))

print(paste('>>> Summary of ', max_FN, ' features finished.', sep = ''))
print('############################################')

FN_best <- max_FN
fea_best <- max_fea
if (Pred_type == 'prob') {
	AUC_mean_best <<- get(paste('AUC_mean_fea', max_FN, sep = ''))
}
F1_mean_best <- get(paste('F1_mean_fea', max_FN, sep = ''))
ACC_mean_best <- get(paste('ACC_mean_fea', max_FN, sep = ''))
Precision_mean_best <- get(paste('Precision_mean_fea', max_FN, sep = ''))
Recall_mean_best <- get(paste('Recall_mean_fea', max_FN, sep = ''))
Spe_mean_best <- get(paste('Spe_mean_fea', max_FN, sep = ''))


## Summary of feature groups with different numbers ( in FN_list) 
for (fn in FN_list[2:len_FN_list]) {

#	assign(paste('fea_', fn, sep = ''), sort_fea_final[1:fn])
	fea_tmp <- get(paste('sort_fea', fn, sep = ''))

	print('############################################')
	print(paste('>>> Top', fn, 'Features are:', sep = ' '))
	print(fea_tmp)

	for (gp in 1:5)
	{
		print('=====================================================')
		print(paste('>>> Group', gp, ':'))
		train <- get(paste('train', gp, sep = ''))
		test <- get(paste('test', gp, sep = ''))
		test$Type <- factor(as.character(test$Type))
		X <- data.frame(train[fea_tmp])
		print('>>> Structrue and dimension of training matrix: ')
		str(X)
		print(dim(X))
		Y <- factor(as.character(train$Type))
		print('>>> Lenght of training set labels: ')
		print(length(Y))
	
		assign(paste('Model_fea', fn, 'gp', gp, sep = '_'), train(y = Y, x = X, method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale")))
		Model_tmp <- get(paste('Model_fea', fn, 'gp', gp, sep = '_'))
		test_mat_tmp <- test[fea_tmp]
		Pred_tmp <- predict(Model_tmp, test_mat_tmp)
		assign(paste('Pred_fea', fn, 'gp', gp, sep = '_'), Pred_tmp)
		if (Pred_type == 'prob') {
			Prob_tmp <- predict(Model_tmp, test_mat_tmp, type = 'prob')
			Pred_prob_tmp <- Prob_tmp[, 2]
			print('>>> Length of predicted prob:')
			print(length(Pred_prob_tmp))
			print('>>> Length of testing set labels:')
			print(length(test$Type))
			Prediction_tmp <- prediction(predictions = Pred_prob_tmp, labels = test$Type)
			Perf.roc_tmp <- performance(Prediction_tmp, measure = 'tpr', x.measure = 'fpr')
			Perf.auc_tmp <- performance(Prediction_tmp, measure = 'auc')
		}
# confusion matrix
		results_tmp <- confusionMatrix(Pred_tmp, test$Type, positive = 'TRUE')
		print('>>> Confusion matrix: ')
		print(confusionMatrix(Pred_tmp, test$Type, positive = 'TRUE'))

		ACC_tmp <- as.numeric(results_tmp$byClass['Balanced Accuracy'])
		Precision_tmp <- as.numeric(results_tmp$byClass['Precision'])
		Recall_tmp <- as.numeric(results_tmp$byClass['Recall']) # Sensitivity
		Spe_tmp <- as.numeric(results_tmp$byClass['Specificity'])
		if (Pred_type == 'prob') {
			AUC_tmp <- unlist(Perf.auc_tmp@y.values)
		}
		F1_tmp <- as.numeric(results_tmp$byClass['F1'])
		print('>>> Precision is:')
		print(Precision_tmp)
		print('>>> Sensitivity is:')
		print(Recall_tmp)
		print('>>> Specificity is:')
		print(Spe_tmp)
		print('>>> ACC is:')
		print(ACC_tmp)
		if (Pred_type == 'prob') {
			print('>>> AUC is:')
			print(AUC_tmp)
		}
		print('>>> F1 score is:')
		print(F1_tmp)

		if (Pred_type == 'prob') {
			if (gp == 1) {
				AUC_list_tmp <- c(AUC_tmp)
			} else {
				AUC_list_tmp <- c(AUC_list_tmp, AUC_tmp)
			}
		}
		if (gp == 1) {
			F1_list_tmp <- c(F1_tmp)
			ACC_list_tmp <- c(ACC_tmp)
			Precision_list_tmp <- c(Precision_tmp)
			Recall_list_tmp <- Recall_tmp
			Spe_list_tmp <- Spe_tmp
		} else {
			F1_list_tmp <- c(F1_list_tmp, F1_tmp)
			ACC_list_tmp <- c(ACC_list_tmp, ACC_tmp)
			Precision_list_tmp <- c(Precision_list_tmp, Precision_tmp)
			Recall_list_tmp <- c(Recall_list_tmp, Recall_tmp)
			Spe_list_tmp <- c(Spe_list_tmp, Spe_tmp)
		}
# running time
		print(difftime(Sys.time(), t1, units = 'sec'))
		print('=====================================================')
		
		if (Pred_type == 'prob') {
			assign(paste('AUC_list_fea', fn, sep = ''), AUC_list_tmp)
		}
		assign(paste('F1_list_fea', fn, sep = ''), F1_list_tmp)
		assign(paste('ACC_list_fea', fn, sep = ''), ACC_list_tmp)
		assign(paste('Precision_list_fea', fn, sep = ''), Precision_list_tmp)
		assign(paste('Recall_list_fea', fn, sep = ''), Recall_list_tmp)
		assign(paste('Spe_list_fea', fn, sep = ''), Spe_list_tmp)
	
	}
	
	
	if (Pred_type == 'prob') {
		assign(paste('AUC_sum_fea', fn, sep = ''), sum(get(paste('AUC_list_fea', fn, sep = ''))))
		assign(paste('AUC_mean_fea', fn, sep = ''), mean(get(paste('AUC_list_fea', fn, sep = ''))))
	}
	assign(paste('F1_sum_fea', fn, sep = ''), sum(get(paste('F1_list_fea', fn, sep = ''))))
	assign(paste('F1_mean_fea', fn, sep = ''), mean(get(paste('F1_list_fea', fn, sep = ''))))
	assign(paste('ACC_sum_fea', fn, sep = ''), sum(get(paste('ACC_list_fea', fn, sep = ''))))
	assign(paste('ACC_mean_fea', fn, sep = ''), mean(get(paste('ACC_list_fea', fn, sep = ''))))
	assign(paste('Precision_sum_fea', fn, sep = ''), sum(get(paste('Precision_list_fea', fn, sep = ''))))
	assign(paste('Precision_mean_fea', fn, sep = ''), mean(get(paste('Precision_list_fea', fn, sep = ''))))
	assign(paste('Recall_sum_fea', fn, sep = ''), sum(get(paste('Recall_list_fea', fn, sep = ''))))
	assign(paste('Recall_mean_fea', fn, sep = ''), mean(get(paste('Recall_list_fea', fn, sep = ''))))
	assign(paste('Spe_sum_fea', fn, sep = ''), sum(get(paste('Spe_list_fea', fn, sep = ''))))
	assign(paste('Spe_mean_fea', fn, sep = ''), mean(get(paste('Spe_list_fea', fn, sep = ''))))
	
	if (Pred_type == 'prob') {
		print('>>> AUC of 5 groups: ')
		print(get(paste('AUC_list_fea', fn, sep = '')))
		print('>>> Mean value of AUC: ')
		print(get(paste('AUC_mean_fea', fn, sep = '')))
	}
	print('>>> F1 of 5 groups: ')
	print(get(paste('F1_list_fea', fn, sep = '')))
	print('>>> Mean value of F1: ')
	print(get(paste('F1_mean_fea', fn, sep = '')))
	print('>>> ACC of 5 groups: ')
	print(get(paste('ACC_list_fea', fn, sep = '')))
	print('>>> Mean value of ACC: ')
	print(get(paste('ACC_mean_fea', fn, sep = '')))
	print('>>> Precision of 5 groups: ')
	print(get(paste('Precision_list_fea', fn, sep = '')))
	print('>>> Mean value of precision: ')
	print(get(paste('Precision_mean_fea', fn, sep = '')))
	print('>>> Sensitivity (Recall) of 5 groups: ')
	print(get(paste('Recall_list_fea', fn, sep = '')))
	print('>>> Mean value of Sensitivity (Recall): ')
	print(get(paste('Recall_mean_fea', fn, sep = '')))
	print('>>> Specificity of 5 groups: ')
	print(get(paste('Spe_list_fea', fn, sep = '')))
	print('>>> Mean value of Specificity: ')
	print(get(paste('Spe_mean_fea', fn, sep = '')))
	
	print(paste('>>> Summary of ', fn, ' features finished.', sep = ''))
	print('############################################')

	print(paste('>>> Compare the performance of ', fn, ' features by index of ', ref_index, ' ... ... ', sep = ''))

	if (Pred_type == 'prob') {
		if (ref_index == 'AUC') {
			if (get(paste('AUC_mean_fea', fn, sep = '')) > AUC_mean_best) {
				FN_best <<- fn
				fea_best <<- fea_tmp
				AUC_mean_best <<- get(paste('AUC_mean_fea', fn, sep = ''))
				F1_mean_best <<- get(paste('F1_mean_fea', fn, sep = ''))
				ACC_mean_best <<- get(paste('ACC_mean_fea', fn, sep = ''))
				Precision_mean_best <<- get(paste('Precision_mean_fea', fn, sep = ''))
				Recall_mean_best <<- get(paste('Recall_mean_fea', fn, sep = ''))
				Spe_mean_best <<- get(paste('Spe_mean_fea', fn, sep = ''))
			}
		}
		if (ref_index == 'F1') {
			if (get(paste('F1_mean_fea', fn, sep = '')) > F1_mean_best) {
				FN_best <<- fn
				fea_best <<- fea_tmp
				AUC_mean_best <<- get(paste('AUC_mean_fea', fn, sep = ''))
				F1_mean_best <<- get(paste('F1_mean_fea', fn, sep = ''))
				ACC_mean_best <<- get(paste('ACC_mean_fea', fn, sep = ''))
				Precision_mean_best <<- get(paste('Precision_mean_fea', fn, sep = ''))
				Recall_mean_best <<- get(paste('Recall_mean_fea', fn, sep = ''))
				Spe_mean_best <<- get(paste('Spe_mean_fea', fn, sep = ''))
			}
		}
	}
	if (Pred_type == 'raw') {
		if (get(paste('F1_mean_fea', fn, sep = '')) > F1_mean_best) { 
			FN_best <<- fn
			fea_best <<- fea_tmp
			F1_mean_best <<- get(paste('F1_mean_fea', fn, sep = ''))
			ACC_mean_best <<- get(paste('ACC_mean_fea', fn, sep = ''))
			Precision_mean_best <<- get(paste('Precision_mean_fea', fn, sep = ''))
			Recall_mean_best <<- get(paste('Recall_mean_fea', fn, sep = ''))
			Spe_mean_best <<- get(paste('Spe_mean_fea', fn, sep = ''))
		}
	}

	print('>>> Now the best feature number is: ')
	print(FN_best)
	print('>>> Best features: ')
	print(fea_best)
	print(paste('>>> Finish the comparation of ', fn, ' features. ', sep = ''))

}


# best features

print('+++++++++++++++++++++++++++++++++++++++++++++++++++++')
print('>>> Best features are:')
print(fea_best)
print('>>> Best feature number is: ')
print(paste('    ', FN_best, sep = ''))
print('>>> Mean value of precision is:')
print(Precision_mean_best)
print('>>> Mean value of sensitivity is:')
print(Recall_mean_best)
print('>>> Mean value of specificity is:')
print(Spe_mean_best)
print('>>> Mean value of ACC is:')
print(ACC_mean_best)
if (Pred_type == 'prob') {
	print('>>> Mean value of AUC is:')
	print(AUC_mean_best)
}
print('>>> Mean value of F1 score is:')
print(F1_mean_best)
print('#END')
print('+++++++++++++++++++++++++++++++++++++++++++++++++++++')


# stop cluster and register sequntial front end
stopCluster(cl); registerDoSEQ();

# save R data
save(list = objects(), file = paste(cell, PM, 'FS.RData', sep = '_'))


### END
