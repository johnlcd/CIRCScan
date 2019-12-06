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
FN_list <- c(as.numeric(seq(2,FN-1)))
args <- commandArgs(T)
if (args[4] != 'all'){
	FN_list <- as.numeric(strsplit(args[4], ',')[[1]])
}
# sort and reverse FN_list (decrease)
FN_list <- rev(sort(FN_list))
cat('>>> Feature number list:\n')
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

cat('>>> The reference index is: ')
print(ref_index)
Pred_type <- args[6]

# set current time
t1 = Sys.time()

# show objects and which libraries were loaded
cat('>>> Objects are:\n')
ls()
cat('>>> Session info:\n')
sessionInfo()

# register parallel front-end
cores <- as.numeric(args[3])
cl <- makeCluster(cores); registerDoParallel(cl)

## function get_result
train_cv <- function(fd) {

	cat('=====================================================\n')
	cat(paste('>>> Fold ', fd, " : \n", sep = ''))
	if (fd %in% 1:9) {
		fd <- paste('Fold0', fd, sep = '')
	} else if (fd == 10) {
		fd <- paste('Fold', fd, sep = '')
	}
	pred_fd <<- pred_tune[pred_tune$Resample==fd,]
	ind_fd <<- pred_fd$rowIndex
	train_fd_mat <<- data_train_sel_mat[-ind_fd,]
	test_fd_mat <<- data_train_sel_mat[ind_fd,]
	ctrl_fd <- trainControl(method = 'none', savePredictions = "all", returnData = T, verboseIter = T, allowParallel = T)
	Model_fd <<- train(y = data_train$Type[-ind_fd], x = train_fd_mat, method = PM, tuneGrid = tunegrid, trControl = ctrl_fd, prob.model = TRUE, preProc = c("center", "scale"))
	Pred_fd <<- predict(Model_fd, test_fd_mat)
	if (Pred_type == 'prob') {
		Prob_fd <<- predict(Model_fd, test_fd_mat, type = 'prob')
		Prob_fd_test <<- predict(Model_fd, test_fd_mat, type = 'prob')
		pp_fd <<- Prob_fd[, 2]
		Prediction_fd <<- prediction(predictions = pp_fd, labels = data_train$Type[ind_fd])
		Perf.auc_fd <<- performance(Prediction_fd, measure = 'auc')
	}
	results_fd <<- confusionMatrix(Pred_fd, data_train$Type[ind_fd], positive = 'TRUE')
	cat('>>> Confusion matrix:\n')
	print(results_fd)
	ACC_fd <<- as.numeric(results_fd$byClass['Balanced Accuracy'])
	Pre_fd <<- as.numeric(results_fd$byClass['Precision'])
	Rec_fd <<- as.numeric(results_fd$byClass['Recall']) # Sensitivity
	Spe_fd <<- as.numeric(results_fd$byClass['Specificity'])
	if (Pred_type == 'prob') {
		AUC_fd <<- unlist(Perf.auc_fd@y.values)
		AUC_list <<- c(AUC_list, AUC_fd)
		if (fd == 'Fold01') {
			AUC_list <<- AUC_list[-1]
		}        
	}
	F1_fd <<- as.numeric(results_fd$byClass['F1'])
	ACC_list <<- c(ACC_list, ACC_fd)
	Pre_list <<- c(Pre_list, Pre_fd)
	Rec_list <<- c(Rec_list, Rec_fd)
	Spe_list <<- c(Spe_list, Spe_fd)
	F1_list <<- c(F1_list, F1_fd)
	if (fd == 'Fold01') {
		ACC_list <<- ACC_list[-1]
		Pre_list <<- Pre_list[-1]
		Rec_list <<- Rec_list[-1]
		Spe_list <<- Spe_list[-1]
		F1_list <<- F1_list[-1]
	}
	cat('>>> Precision is: \n')
	print(Pre_fd)
	cat('==> New Precision list: \n')
	print(Pre_list)
	cat('>>> Sensitivity is: \n')
	print(Rec_fd)
	cat('==> New Sensitivity list: \n')
	print(Rec_list)
	cat('>>> Specificity is: \n')
	print(Spe_fd)
	cat('==> New Specificity list: \n')
	print(Spe_list)
	cat('>>> ACC is: \n')
	print(ACC_fd)
	cat('==> New ACC list: \n')
	print(ACC_list)
	if (Pred_type == 'prob') {
		cat('>>> AUC is: \n')
		print(AUC_fd)
		cat('==> New AUC list: \n')
		print(AUC_list)
	}
	cat('>>> F1 score is: \n')
	print(F1_fd)
	cat('==> New F1 list: \n')
	print(F1_list)

# running time
	print(difftime(Sys.time(), t1, units = 'sec'))
	cat('=====================================================\n')

}

# model traning with different feature number
ctrl <- trainControl(method = "none", savePredictions = "all", returnData = T, verboseIter = T, allowParallel = T)
ctrl_cv <- trainControl(method = 'cv', number = 10, savePredictions = "all", returnData = T, returnResamp = 'all', 
						verboseIter = T, allowParallel = T)

cat(paste('>>> Start the feature selection (cross validition) of all ', FN, ' features ... ... \n', sep = ''))
cat(paste('>>> All', FN, 'Features are: \n', sep = ' '))
print(fea_all)
cat('############################################\n')
cat(paste('>>> Top ', FN, ' features: \n', sep = ''))
print(fea_all)
assign(paste('Model_FN', FN, sep = ''), Model_fea_all)
cat('>>> Model summary: \n')
print(Model_fea_all)

# Get results of resample
all_best_tune <- Model_fea_all$bestTune
tunegrid <- expand.grid(all_best_tune)
tune_met_all <- colnames(all_best_tune)
cat('>>> Best tune: \n')
print(all_best_tune)
# tuned pred 
pred_tune <- Model_fea_all$pred
for (i in 1:length(tune_met_all)) {
	met <- tune_met_all[i]
	met_val <- all_best_tune[,met]
	pred_tune <- pred_tune[pred_tune[,met]==met_val,]
}
cat('>>> Head of tune predict results: \n')
print(head(pred_tune))
# initialize the performance list
ACC_list <- c(0)
Pre_list <- c(0)
Rec_list <- c(0)
Spe_list <- c(0)
F1_list <- c(0)
ACC_test <- c()
Pre_test <- c()
Rec_test <- c()
Spe_test <- c()
FPR_test <- c()
F1_test <- c()
if (Pred_type == 'prob') {
	AUC_list <- c(0)
	AUC_test <- c()
}

data_train_all_mat <- data_train_mat
data_train_sel_mat <- data_train_all_mat
for (fold in 1:10) {
	Model_tmp <- Model_fea_all
	train_cv(fold)
}

pred_test_all <- predict(Model_fea_all, data_test_mat)
if (Pred_type == 'prob') {
	Prob_test_all <<- predict(Model_fea_all, data_test_mat, type = 'prob')
	pp_test_all <<- Prob_test_all[, 2]
	Prediction_test_all <<- prediction(predictions = pp_test_all, labels = data_test$Type)
	Perf.auc_test_all <<- performance(Prediction_test_all, measure = 'auc')
}
results_test_all <<- confusionMatrix(pred_test_all, data_test$Type, positive = 'TRUE')
cat('>>> Confusion matrix:\n')
print(results_test_all)
ACC_test_all <<- as.numeric(results_test_all$byClass['Balanced Accuracy'])
Pre_test_all <<- as.numeric(results_test_all$byClass['Precision'])
Rec_test_all <<- as.numeric(results_test_all$byClass['Recall']) # Sensitivity
Spe_test_all <<- as.numeric(results_test_all$byClass['Specificity'])
FPR_test_all <<- 1 - Spe_test_all
if (Pred_type == 'prob') {
	AUC_test_all <<- unlist(Perf.auc_test_all@y.values)
	AUC_test <<- c(AUC_test, AUC_test_all)
}
F1_test_all <<- as.numeric(results_test_all$byClass['F1'])
ACC_test <<- c(ACC_test, ACC_test_all)
Pre_test <<- c(Pre_test, Pre_test_all)
Rec_test <<- c(Rec_test, Rec_test_all)
Spe_test <<- c(Spe_test, Spe_test_all)
FPR_test <<- c(FPR_test, FPR_test_all)
F1_test <<- c(F1_test, F1_test_all)

fold_list <- paste('Fold0', seq(1,9), sep = '')
fold_list <- c(fold_list, 'Fold10')
cell_list <- rep(cell, 10)
model_list <- rep(PM, 10)
fn_list <- rep(FN, 10)
if (Pred_type == 'prob') {
	AUC_result_df <- data.frame(cbind(AUC_list, cell_list, model_list, fn_list, fold_list))
	colnames(AUC_result_df) <- c('AUC', 'Cell_type', 'Model', 'Feature_num', 'Fold')
	write.table(AUC_result_df, paste(cell, PM, 'AUC_cv.txt', sep = '_'), row.names = F, col.names = T, quote = F, sep = '\t')
}
F1_result_df <- data.frame(cbind(F1_list, cell_list, model_list, fn_list, fold_list))
colnames(F1_result_df) <- c('F1', 'Cell_type', 'Model', 'Feature_num', 'Fold')
write.table(F1_result_df, paste(cell, PM, 'F1_cv.txt', sep = '_'), row.names = F, col.names = T, quote = F, sep = '\t')

ACC_list_best <- ACC_list
Pre_list_best <- Pre_list
Rec_list_best <- Rec_list
Spe_list_best <- Spe_list
F1_list_best <- F1_list
ACC_mean <- mean(ACC_list)
Pre_mean <- mean(Pre_list)
Rec_mean <- mean(Rec_list)
Spe_mean <- mean(Spe_list)
F1_mean <- mean(F1_list)
F1_mean_best <- F1_mean
FN_best <- FN
fea_best<- fea_all

assign(paste('ACC_list_fea', FN, sep = '_'), ACC_list)
assign(paste('Pre_list_fea', FN, sep = '_'), Pre_list)
assign(paste('Rec_list_fea', FN, sep = '_'), Rec_list)
assign(paste('Spe_list_fea', FN, sep = '_'), Spe_list)
assign(paste('F1_list_fea', FN, sep = '_'), F1_list)
assign(paste('ACC_mean_fea', FN, sep = '_'), ACC_mean)
assign(paste('Pre_mean_fea', FN, sep = '_'), Pre_mean)
assign(paste('Rec_mean_fea', FN, sep = '_'), Rec_mean)
assign(paste('Spe_mean_fea', FN, sep = '_'), Spe_mean)
assign(paste('F1_mean_fea', FN, sep = '_'), F1_mean)
if (Pred_type == 'prob') {
	AUC_list_best <- AUC_list
	AUC_mean <- mean(AUC_list)
	AUC_mean_best <- AUC_mean
	assign(paste('AUC_list_fea', FN, sep = '_'), AUC_list)
	assign(paste('AUC_mean_fea', FN, sep = '_'), AUC_mean)
}


## Summary of feature groups with different numbers ( in FN_list) 
for (fn in FN_list[1:len_FN_list]) {

	fea_tmp <- get(paste('sort_fea', fn, sep = ''))

	cat('############################################\n')
	cat(paste('>>> Top', fn, 'Features are:\n', sep = ' '))
	print(fea_tmp)
	cat('>>> Model summary: \n')
	Model_tmp <- get(paste('Model_FN', fn, sep = ''))
	print(Model_tmp)

	# Get results of resample
	best_tune <- Model_tmp$bestTune
	tunegrid <- expand.grid(best_tune)
	tune_met <- colnames(best_tune)
	cat('>>> Best tune: \n')
	print(best_tune)

	# tuned pred 
	pred_tune <- Model_tmp$pred
	for (i in 1:length(tune_met)) {
			met <- tune_met[i]
		met_val <- best_tune[,met]
			pred_tune <- pred_tune[pred_tune[,met]==met_val,]
	}
	cat('>>> Head of tune predict results: \n')
	print(head(pred_tune))
	# initialize the performance list
	ACC_list <- c(0)
	Pre_list <- c(0)
	Rec_list <- c(0)
	Spe_list <- c(0)
	F1_list <- c(0)
	if (Pred_type == 'prob') {
		AUC_list <- c(0)
	}

	data_train_mat_tmp <- data_train[, fea_tmp]
	data_train_sel_mat <- data_train_mat_tmp
	for (fold in 1:10) {
		train_cv(fold)
	}
	fn_list <- rep(fn, 10)
	if (Pred_type == 'prob') {
		AUC_result_df <- data.frame(cbind(AUC_list, cell_list, model_list, fn_list, fold_list))
		write.table(AUC_result_df, paste(cell, PM, 'AUC_cv.txt', sep = '_'), row.names = F, col.names = F, append = T, quote = F, sep = '\t')
	}
	F1_result_df <- data.frame(cbind(F1_list, cell_list, model_list, fn_list, fold_list))
	write.table(F1_result_df, paste(cell, PM, 'F1_cv.txt', sep = '_'), row.names = F, col.names = F, append = T, quote = F, sep = '\t')

	ACC_mean <- mean(ACC_list)
	Pre_mean <- mean(Pre_list)
	Rec_mean <- mean(Rec_list)
	Spe_mean <- mean(Spe_list)
	F1_mean <- mean(F1_list)
	
	assign(paste('ACC_list_fea', fn, sep = '_'), ACC_list)
	assign(paste('Pre_list_fea', fn, sep = '_'), Pre_list)
	assign(paste('Rec_list_fea', fn, sep = '_'), Rec_list)
	assign(paste('Spe_list_fea', fn, sep = '_'), Spe_list)
	assign(paste('F1_list_fea', fn, sep = '_'), F1_list)
	assign(paste('ACC_mean_fea', fn, sep = '_'), ACC_mean)
	assign(paste('Pre_mean_fea', fn, sep = '_'), Pre_mean)
	assign(paste('Rec_mean_fea', fn, sep = '_'), Rec_mean)
	assign(paste('Spe_mean_fea', fn, sep = '_'), Spe_mean)
	assign(paste('F1_mean_fea', fn, sep = '_'), F1_mean)
	if (Pred_type == 'prob') {
		AUC_mean <- mean(AUC_list)
		assign(paste('AUC_list_fea', fn, sep = '_'), AUC_list)
		assign(paste('AUC_mean_fea', fn, sep = '_'), AUC_mean)
	}

	if (Pred_type == 'prob') {
		if (ref_index == 'AUC') {
			if (AUC_mean > AUC_mean_best) {
				FN_best <<- fn
				fea_best <<- fea_tmp
				AUC_mean_best <<- AUC_mean
				F1_mean_best <<- F1_mean
			}
		} else if (ref_index == 'F1') {
			if (F1_mean > F1_mean_best) {
				FN_best <<- fn
				fea_best <<- fea_tmp
				AUC_mean_best <<- AUC_mean
				F1_mean_best <<- F1_mean
			}
		}
	} else if (Pred_type == 'raw') {
		if (F1_mean > F1_mean_best) { 
			FN_best <<- fn
			fea_best <<- fea_tmp
			F1_mean_best <<- F1_mean
		}
	}

	cat('>>> Now the best feature number is: \n')
	print(FN_best)
	cat('>>> Best features: \n')
	print(fea_best)
	Model_best <<- get(paste('Model_FN', FN_best, sep = ''))
	cat('>>> Best model: \n')
	print(Model_best)
	cat(paste('>>> Finish the comparation of ', fn, ' features. \n', sep = ''))

	# evaluation of model performance
	pred_test_tmp <- predict(Model_tmp, data_test_mat)
	if (Pred_type == 'prob') {
		Prob_test_tmp <<- predict(Model_tmp, data_test_mat, type = 'prob')
		pp_test_tmp <<- Prob_test_tmp[, 2]
		Prediction_test_tmp <<- prediction(predictions = pp_test_tmp, labels = data_test$Type)
		Perf.auc_test_tmp <<- performance(Prediction_test_tmp, measure = 'auc')
	}
	results_test_tmp <<- confusionMatrix(pred_test_tmp, data_test$Type, positive = 'TRUE')
	cat('>>> Confusion matrix:\n')
	print(results_test_tmp)
	ACC_test_tmp <<- as.numeric(results_test_tmp$byClass['Balanced Accuracy'])
	Pre_test_tmp <<- as.numeric(results_test_tmp$byClass['Precision'])
	Rec_test_tmp <<- as.numeric(results_test_tmp$byClass['Recall']) # Sensitivity
	Spe_test_tmp <<- as.numeric(results_test_tmp$byClass['Specificity'])
	FPR_test_tmp <<- 1 - Spe_test_tmp
	if (Pred_type == 'prob') {
		AUC_test_tmp <<- unlist(Perf.auc_test_tmp@y.values)
		AUC_test <<- c(AUC_test, AUC_test_tmp)
	}
	F1_test_tmp <<- as.numeric(results_test_tmp$byClass['F1'])
	ACC_test <<- c(ACC_test, ACC_test_tmp)
	Pre_test <<- c(Pre_test, Pre_test_tmp)
	Rec_test <<- c(Rec_test, Rec_test_tmp)
	Spe_test <<- c(Spe_test, Spe_test_tmp)
	FPR_test <<- c(FPR_test, FPR_test_tmp)
	F1_test <<- c(F1_test, F1_test_tmp)


}


# best features

cat('+++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
cat('>>> Best features are: \n')
print(fea_best)
cat('>>> Best feature number is: \n')
cat(paste('    ', FN_best, "\n", sep = ''))
cat('>>> Mean value of precision is:\n')
Pre_mean_best <- mean(Pre_list_best)
print(Pre_mean_best)
cat('>>> Mean value of sensitivity is:\n')
Rec_mean_best <- mean(Rec_list_best)
print(Rec_mean_best)
cat('>>> Mean value of specificity is:\n')
Spe_mean_best <- mean(Spe_list_best)
print(Spe_mean_best)
cat('>>> Mean value of FPR is: \n')
FPR_mean_best <- 1 - Spe_mean_best
print(FPR_mean_best)
cat('>>> Mean value of ACC is:\n')
ACC_mean_best <- mean(ACC_list_best)
print(ACC_mean_best)
if (Pred_type == 'prob') {
	cat('>>> Mean value of AUC is:\n')
	AUC_mean_best <- mean(AUC_list_best)
	print(AUC_mean_best)
}
cat('>>> Mean value of F1 score is:\n')
F1_mean_best <- mean(F1_list_best)
print(F1_mean_best)
cat('>>> Summary of best model: \n')
print(Model_best)

# Model performance
FN_test <- c(FN, FN_list)
cell_test <- rep(cell, length(FN_test))
model_test <- rep(PM, length(FN_test))
if (Pred_type == 'prob') {
	perf_test_df <- data.frame(cbind(AUC_test, FPR_test, F1_test, cell_test, model_test, FN_test))
	colnames(perf_test_df) <- c('AUC', 'FPR', 'F1', 'Cell_type', 'Model', 'Feature_num')
} else {
	perf_test_df <- data.frame(cbind(FPR_test, F1_test, cell_test, model_test, FN_test))
	colnames(perf_test_df) <- c('FPR', 'F1', 'Cell_type', 'Model', 'Feature_num')
}
write.table(perf_test_df, paste(cell, PM, 'perf_test.txt', sep = '_'), row.names = F, col.names = T, quote = F, sep = '\t')
cat(">>> Summary of model performance in testing set:\n")
print(perf_test_df)

cat('#END\n')
cat('+++++++++++++++++++++++++++++++++++++++++++++++++++++\n')


# stop cluster and register sequntial front end
stopCluster(cl); registerDoSEQ();

# save R data
rm(list = 'args')
save(list = objects(), file = paste(cell, PM, 'FS.RData', sep = '_'))


### END

