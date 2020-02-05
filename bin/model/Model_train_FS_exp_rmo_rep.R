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
library(data.table)

# set arguments
args <- commandArgs(T)
cell <- args[1]
PM <- args[2]
cores <- as.numeric(args[3])
part <- args[4]

# load data_train
file <- paste(cell, 'exp_train', sep = '_')
data_train_all <- read.table(file, head = T)
data_train_all <- data.frame(data_train_all)
data_train_raw <- data_train_all
data_train_all$SRPBM <- log2(data_train_all$SRPBM)
summary(data_train_all)
n_all <- dim(data_train_all)[1]
all_num <- n_all
FN <- dim(data_train_all)[2] - 5
col_name <- names(data_train_all)
fea_all <- col_name[5:(FN+4)]
cat(">>> Remove outlier circRNAs data points ...\n")
order_srpbm <- order(data_train_all[,'SRPBM'], decreasing = F) # remove outlier circRNAs (highest and lowest 30)
rmn <- 30
rm_list <- as.data.table(data_train_all)[,Intron_pair][c((1:rmn),((all_num-rmn+1):all_num))]
write.table(as.data.table(matrix(rm_list,length(rm_list))), "RM_IP.list", col.name=F, quote=F, sep='\t')
data_train_all <- data_train_all[order_srpbm[-c((1:rmn),((all_num-rmn+1):all_num))],]
set.seed(123)
inTraining_rep <- createDataPartition(data_train_all$SRPBM, p = .9, list = FALSE)
data_train_rep <- data_train_all[inTraining_rep,]
data_test <- data_train_all[-inTraining_rep,]
inTraining1 <- createDataPartition(data_train_rep$SRPBM, p = 1/5, list = FALSE)
data_train1 <- data_train_rep[inTraining1,]
out_train1 <- data_train_rep[-inTraining1,]
inTraining2 <- createDataPartition(out_train1$SRPBM, p = 1/4, list = FALSE)
data_train2 <- out_train1[inTraining2,]
out_train2 <- out_train1[-inTraining2,]
inTraining3 <- createDataPartition(out_train2$SRPBM, p = 1/3, list = FALSE)
data_train3 <- out_train2[inTraining3,]
out_train3 <- out_train2[-inTraining3,]
inTraining4 <- createDataPartition(out_train3$SRPBM, p = 1/2, list = FALSE)
data_train4 <- out_train3[inTraining4,]
data_train5 <- out_train3[-inTraining4,]
#set.seed(seed)
data_train <- get(paste("data_train", part, sep=""))
if (PM == 'glm') {
	data_train$SRPBM <- data_train$SRPBM/10 # value = log2(SRPBM)/10
}
summary(data_train)
cat('>>> Dimension of training data matrix: ( replication', part, ') \n')
dim(data_train)
#all_num <- dim(data_train)[1]
n <- dim(data_train)[1]
data_train_mat <- data_train[,5:(FN+4)]
data_test_mat <- data_test[,5:(FN+4)]
y_scale_test <- mean(data_test$SRPBM)

# show which libraries were loaded  
cat('>>> Session Info: \n')
sessionInfo()

# register parallel front-end
cl <- makeCluster(cores); registerDoParallel(cl)

ctrl <- trainControl(method = 'cv', number = 10 , savePredictions = "all", returnData = T, returnResamp = 'all', 
					 verboseIter = T, allowParallel = T)

# Normalize RMSE (mean of y)
scale_fd <- function(model) {

	y_scale <- c()
	num <- c(paste("0",1:9,sep=''),"10")
	for (i in num) {
	fd <- paste("Fold",i,sep='')
	pred_tmp <- model$pred
	ind_tmp <- unique(pred_tmp$rowIndex[which(pred_tmp$Resample==fd)])
	mean_tmp <- mean(model$trainingData$.outcome[ind_tmp])
	y_scale <- c(y_scale,mean_tmp)
	}
	return(y_scale)

}


# Train best function
feature_sel <- function(fn) {
	('=====================================================\n')
	cat('>>> Sorted features: \n')
	print(sort_fea)
	new_fea_list <<- sort_fea[1: (length(sort_fea)-1)]
	cat('>>> Selected new features: \n')
	print(new_fea_list)
	cat('>>> Feature number input: \n')
	print(fn)
	if (fn == length(new_fea_list)) {
		cat('>>> Temp feature number EQUAL to sorted feature number, PASS ==>> \n')
#		print(paste('>>> Top ', fn+1, ' features ==> ', sep = ''))
#		print(sort_fea)
	} else {
		print('>>> Incoordinate temp feature number and sorted feature number ! ! ! \n')
	}
	cat('>>> Best features: \n')
	print(fea_best)
	
	data_train_mat_tmp <- data_train_mat[, new_fea_list]
	## train model
	Model_tmp <- train(y = data_train$SRPBM, x = data_train_mat_tmp, 
				   method = PM, trControl = ctrl, metric = "RMSE", importance = T, preProc = c("center", "scale"))
	
	# Get results of resample
	cat('>>> summary of resample: \n')
	resample_tmp <- Model_tmp$resample
	print(resample_tmp)
	best_tune <- Model_tmp$bestTune
	tunegrid <- expand.grid(best_tune)
	tune_met <- colnames(best_tune)
	cat('>>> Best tune: \n')
	print(best_tune)

	# tuned resample
	resample_tune <- resample_tmp
	for (i in 1:length(tune_met)) {
		met <- tune_met[i]
		met_val <- best_tune[,met]
		resample_tune <- resample_tune[resample_tune[,met]==met_val,]
	}
	cat('>>> Resamples of best tune: \n')
	print(resample_tune)

	## Cross-validation resample
	cv_perf_df_tmp <- cbind(resample_tune, rep(fn, 10), rep(cell, 10))
	colnames(cv_perf_df_tmp) <- c(colnames(resample_tune), 'Fea_number', 'Cell_type')
	cv_perf_df_tmp$scale_y <- scale_fd(Model_tmp)
	cv_perf_df_tmp$RMSE_norm <- cv_perf_df_tmp$RMSE/cv_perf_df_tmp$scale_y
	cv_perf_df <<- rbind(cv_perf_df, cv_perf_df_tmp)
	cat('>>> New data frame of resample performance: \n')
	print(cv_perf_df)
	write.table(cv_perf_df_tmp, paste(cell, PM, 'cv_perf', sep = '_'), row.names = F, col.names = F, 
				append = T, quote = F, sep = '\t')
	
	## model performance (RMSE, R2)
	if (PM == 'rf') {
		mse_tmp <- mse(data_train$SRPBM, Model_tmp$finalModel$predicted)
		rmse_tmp <- sqrt(mse_tmp)
		nrmse_tmp <- rmse_tmp/mean(data_train$SRPBM)
		SSE_tmp <- sum((Model_tmp$finalModel$predicted - data_train$SRPBM)^2)
		SST_tmp <- sum((data_train$SRPBM - mean(data_train$SRPBM)) ^ 2)
		R2_tmp <- 1- SSE_tmp/SST_tmp
		PCC_tmp <- cor.test(Model_tmp$finalModel$predicted, data_train$SRPBM,method = "pearson")
		cat('>>> Pearson\'s r (PCC): \n')
		print(PCC_tmp)
		PCC_tmp <- PCC_tmp$estimate[[1]]
	} else {
		tmp_pred <- Model_tmp$pred
		obs_tmp <- Model_tmp$pred$obs
		for (i in 1:length(tune_met)) {
			met <- tune_met[i]
			met_val <- best_tune[,met]
			tmp_pred <- tmp_pred[tmp_pred[,met]==met_val,]
			obs_tmp <- tmp_pred$obs
			pred_tmp <- tmp_pred$pred
		}
		mse_tmp <- mse(obs_tmp, pred_tmp)
		rmse_tmp <- sqrt(mse_tmp)
		nrmse_tmp <- rmse_tmp/mean(data_train$SRPBM)
		SSE_tmp <- sum((pred_tmp - obs_tmp) ^ 2)
		SST_tmp <- sum((data_train$SRPBM - mean(data_train$SRPBM)) ^ 2)
		R2_tmp <- 1 - SSE_tmp/SST_tmp
	}
	cat('>>> Model performance ==> \n')
	cat(paste('>>> Total RMSE: ', rmse_tmp, "\n",sep = ''))
	cat(paste('>>> Normalized RMSE: ', nrmse_tmp, "\n", sep = ''))
	cat(paste('>>> Total R2: ', R2_tmp, "\n", sep = ''))
	
	## feature importance
	imp_tmp <- varImp(Model_tmp)
	imp_tmp <- round(imp_tmp$importance, 2)
	imp_tmp <- data.frame(imp_tmp)
	order_imp_tmp <- order(imp_tmp[,'Overall'], decreasing = T)
	sort_fea_tmp <- rownames(imp_tmp)[order_imp_tmp]
	sort_imp_tmp <- data.frame(imp_tmp[,'Overall'][order_imp_tmp])
	rownames(sort_imp_tmp) <- sort_fea_tmp
	colnames(sort_imp_tmp) <- 'Importance'
	cat('>>> Sorted importance of MSE ==> \n')
	print(sort_imp_tmp)
	cat('>>> Summary of model: \n')
	print(Model_tmp)

	## compare to previous best model
	if (nrmse_tmp < nrmse_best) {
		Model_best <<- Model_tmp
		rmse_best <<- rmse_tmp
		nrmse_best <<- nrmse_tmp
		SSE_best <<- SSE_tmp
		SST_best <<- SST_tmp
		R2_best <<- R2_tmp
		if (PM == 'rf') {
			PCC_best <<- PCC_tmp
		}
		fea_best <<- new_fea_list
		fn_best <<- fn
		tune_best <<- best_tune
	}

	sort_fea <<- rownames(sort_imp_tmp)
	cat('   New rank of features: \n')
	print(sort_fea)

	# evaluation of model performance
	pred_exp_test_tmp <- predict(Model_tmp, data_test_mat)
	pred_exp_test_tmp <- as.vector(pred_exp_test_tmp)
	mmse_test_tmp <- mse(data_test$SRPBM,pred_exp_test_tmp)
	rmse_test_tmp <- sqrt(mmse_test_tmp)
	nrmse_test_tmp <- rmse_test_tmp/y_scale_test
	SSE_test_tmp <- sum((pred_exp_test_tmp - data_test$SRPBM) ^ 2)
	SST_test_tmp <- sum((data_test$SRPBM - y_scale_test) ^ 2)
	R2_test_tmp <- 1 - SSE_test_tmp/SST_test_tmp
	rmse_test <<- c(rmse_test, rmse_test_tmp)
	nrmse_test <<- c(nrmse_test, nrmse_test_tmp)
	R2_test <<- c(R2_test, R2_test_tmp)
	if (PM == 'rf') {
		PCC_test_tmp <- cor.test(data_test$SRPBM,pred_exp_test_tmp,method = "pearson")
		PCC_test_tmp <- PCC_test_tmp$estimate[[1]]
		PCC_test <<- c(PCC_test, PCC_test_tmp)
	}


	if (fn == 3) {
		save(list = objects(), file=paste(cell, PM, 'FS_exp.RData', sep = '_'))
	}

	cat('=====================================================\n')

}


cat('>>> [1] Train model with all features ==> \n')
Model_all <- train(y = data_train$SRPBM, x = data_train_mat, 
			   method = PM, trControl = ctrl, metric = "RMSE", importance = T, preProc = c("center", "scale"))
Model_best <- Model_all

# Get results of resample
cat('>>> Summary of resample (all features): \n')
resample_all <- Model_all$resample
print(resample_all)
best_tune_all <- Model_all$bestTune
tune_best <- best_tune_all
tunegrid <- expand.grid(best_tune_all)
tune_met_all <- colnames(best_tune_all)
cat('>>> Best tune of all features: \n')
print(best_tune_all)

# tuned resample
resample_tune_all <- resample_all
for (i in 1:length(tune_met_all)) {
	met_all <- tune_met_all[i]
	met_val_all <- best_tune_all[,met_all]
	resample_tune_all <- resample_tune_all[resample_tune_all[,met_all]==met_val_all,]
}
cat('>>> Resamples of best tune: \n')
print(resample_tune_all)

## Cross-validation resample
cv_perf_df <- cbind(resample_tune_all, rep(FN, 10), rep(cell, 10))
colnames(cv_perf_df) <- c(colnames(resample_tune_all), 'Fea_number', 'Cell_type')
cv_perf_df$scale_y <- scale_fd(Model_all)
cv_perf_df$RMSE_norm <- cv_perf_df$RMSE/cv_perf_df$scale_y
write.table(cv_perf_df, paste(cell, PM, 'cv_perf', sep = '_'), row.names = F, col.names = T, quote = F, sep = '\t')

## model performance of all feature (RMSE, R2)
if (PM == 'rf') {
	mmse <- mse(data_train$SRPBM, Model_all$finalModel$predicted)
	rmse <- sqrt(mmse)
	nrmse <- rmse/mean(data_train$SRPBM) # normalize RMSE by mean of y
	SSE <- sum((Model_all$finalModel$predicted - data_train$SRPBM)^2)
	SST <- sum((data_train$SRPBM - mean(data_train$SRPBM)) ^ 2)
	R2 <- 1- SSE/SST
	PCC <- cor.test(Model_all$finalModel$predicted, data_train$SRPBM,method = "pearson")
	PCC_best <- PCC
	cat('>>> Pearson\'s r (PCC): \n')
	print(PCC)
} else {
	all_pred <- Model_all$pred
	for (i in 1:length(tune_met_all)) {
		met_all <- tune_met_all[i]
		met_val_all <- best_tune_all[,met_all]
		all_pred <- all_pred[all_pred[,met_all]==met_val_all,]
		obs_all <- all_pred$obs
		pred_all <- all_pred$pred
	}
	mmse <- mse(obs_all, pred_all)
	rmse <- sqrt(mmse)
	nrmse <- rmse/mean(data_train$SRPBM)
	SSE <- sum((pred_all - obs_all) ^ 2)
	SST <- sum((data_train$SRPBM - mean(data_train$SRPBM)) ^ 2)
	R2 <- 1 - SSE/SST
}
rmse_best <- rmse
nrmse_best <- nrmse
SSE_best <- SSE
SST_best <- SST
R2_best <- R2

cat('>>> Model performance ==> \n')
cat(paste('>>> Total RMSE: ', rmse, "\n", sep = ''))
cat(paste('>>> Normalized RMSE: ', nrmse, "\n", sep = ''))
cat(paste('>>> Total R2: ', R2, "\n", sep = ''))

## feature importance
imp <- varImp(Model_all)
imp <- round(imp$importance, 2)
imp <- data.frame(imp)
order_imp <- order(imp[,'Overall'], decreasing = T)
sort_fea_all <- rownames(imp)[order_imp]
sort_imp <- data.frame(imp[,'Overall'][order_imp])
rownames(sort_imp) <- sort_fea_all
colnames(sort_imp) <- 'Importance'
sort_imp <- data.frame(sort_imp)
colnames(sort_imp) <- c('Importance')
#max_fea <- sort_fea_all[1:max_fn]
fea_best <- sort_fea_all
fn_best <- length(fea_best)
sort_fea <- sort_fea_all
write.table(sort_imp, paste(cell, PM, 'sort_Imp', sep = '_'), row.names = T, col.names = F, quote = F, sep = '\t')
cat('>>> Sorted importance of MSE ==> \n')
print(sort_imp)
if (PM == 'rf') {
	pdf(file = paste(cell, PM, "Imp.pdf", sep = "_"))
	varImpPlot(Model_all$finalModel, type = 1, main = 'Feature Importance')
	dev.off()
}
cat('>>> Summary of model: \n')
print(Model_all)


# initialize the performance list
rmse_test <- c()
nrmse_test <- c()
R2_test <- c()
PCC_test <- c()
# evaluation of performance in testing set
pred_exp_test_all <- predict(Model_all, data_test_mat)
pred_exp_test_all <- as.vector(pred_exp_test_all)
mmse_test_all <- mse(data_test$SRPBM,pred_exp_test_all)
rmse_test_all <- sqrt(mmse_test_all)
nrmse_test_all <- rmse_test_all/y_scale_test
SSE_test_all <- sum((pred_exp_test_all - data_test$SRPBM) ^ 2)
SST_test_all <- sum((data_test$SRPBM - y_scale_test) ^ 2)
R2_test_all <- 1 - SSE_test_all/SST_test_all
rmse_test <- c(rmse_test, rmse_test_all)
nrmse_test <- c(nrmse_test, nrmse_test_all)
R2_test <- c(R2_test, R2_test_all)
if (PM == 'rf') {
	PCC_test_all <- cor.test(data_test$SRPBM,pred_exp_test_all,method = "pearson")
	PCC_test_all <- PCC_test_all$estimate[[1]]
	PCC_test <- c(PCC_test, PCC_test_all)
}

cat('>>> [2] Feature selsction ==> \n')
for (fn in seq(length(sort_fea_all)-1,2)) {
	feature_sel(fn)
}
cat('>>> Feature selsction finished. \n')


cat('>>> [3] Summary of feature selection ==> \n')
cat('>>> Best features: \n')
print(fea_best)
cat('>>> Number of best feature is: \n')
print(fn_best)
cat('>>> Model summary: \n')
print(Model_best)
cat('>>> Performance of best model ==> \n')
cat(paste('>>> Total RMSE: ', rmse_best, "\n", sep = ''))
cat(paste('>>> Normalized RMSE: ', nrmse_best, "\n", sep = ''))
cat(paste('>>> Total R2: ', R2_best, "\t", sep = ''))
if (PM == 'rf') {
	cat('>>> Pearson\'s r (PCC): \n')
	print(PCC_best)
}
tune_met_best <- colnames(tune_best)
cat('>>> Best tune of final model: \n')
print(tune_best)

# Model performance
FN_test <- c(FN,seq(length(sort_fea_all)-1,2))
cell_test <- rep(cell, length(FN_test))
model_test <- rep(PM, length(FN_test))
if (PM == 'rf') {
	perf_test_df <- data.frame(cbind(rmse_test, nrmse_test, R2_test, PCC_test, cell_test, model_test, FN_test))
	colnames(perf_test_df) <- c ("RMSE","RMSE_norm","R2","PCC","Cell_type","Model","Feature_num")
} else {
	perf_test_df <- data.frame(cbind(rmse_test, nrmse_test, R2_test, cell_test, model_test, FN_test))
	colnames(perf_test_df) <- c ("RMSE","RMSE_norm","R2","Cell_type","Model","Feature_num")
}
write.table(perf_test_df, paste(cell, PM, 'perf_test_reg.txt', sep = '_'), row.names = F, col.names = T, quote = F, sep = '\t')


## write the observed and predicted expression to file
best_pred_test <- predict(Model_best, data_test_mat)
best_pred_test <- as.vector(best_pred_test)
data_test <- as.data.table(data_test)
obs_pred_exp <- data.frame(cbind(data_test[,"Intron_pair"], data_test$SRPBM, best_pred_test))
colnames(obs_pred_exp) <- c('Intron_pair','true_exp', 'pred_exp')
write.table(obs_pred_exp, paste(cell, PM, 'train_pred_exp', sep = '_'), row.names = F, col.names = T, quote = F, sep = '\t')

# stop cluster and register sepuntial front end
stopCluster(cl); registerDoSEQ()

# save R data
rm(list = 'args')
save(list = objects(), file=paste(cell, PM, 'FS_exp.RData', sep = '_'))
cat('>>> Task DONE !!!\n')


### END
