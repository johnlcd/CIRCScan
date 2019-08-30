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
args <- commandArgs(T)
cell <- args[1]
PM <- args[2]
cores <- as.numeric(args[3])
seed <- as.numeric(args[4])

# load data_train
file <- paste(cell, 'exp_train', sep = '_')
data_train_all <- read.table(file, head = T)
data_train_all <- data.frame(data_train_all)
data_train_raw <- data_train_all
data_train_all$SRPBM <- log2(data_train_all$SRPBM)
data_train <- data_train_all
all_num <- dim(data_train)[1]
cat(">>> Remove outlier circRNAs data points ...\n")
order_srpbm <- order(data_train[,'SRPBM'], decreasing = F)
data_train[order_srpbm[-c((1:20),(all_num-20:all_num))],] # remove outlier circRNAs (highest and lowest 20) 
if (PM == 'glm') {
	data_train$SRPBM <- data_train$SRPBM/10 # value = log2(SRPBM)/10
}
summary(data_train)
set.seed(seed)
cat('>>> Dimension of data matrix: ')
dim(data_train)
n <- dim(data_train)[1]
FN <- dim(data_train)[2] - 5
col_name <- names(data_train)
fea_all <- col_name[5:(FN+4)]
data_train_mat <- data_train[,5:(FN+4)]

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
	cat('=====================================================\n')
	cat('>>> Sorted features: \n')
	print(sort_fea)
	new_fea_list <<- sort_fea[1: (length(sort_fea)-1)]
	cat('>>> Selected new features: \n')
	print(new_fea_list)
	cat('>>> Feature number input: \n')
	print(fn)
	if (fn == length(new_fea_list)) {
		cat('Temp feature number EQUAL to sorted feature number, PASS ==>> \n')
	} else {
		cat('>>> Incoordinate temp feature number and sorted feature number ! ! ! \n')
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
		SSE_tmp <- sum((Model_tmp$finalModel$predictedi - data_train$SRPBM)^2)
		SST_tmp <- sum((data_train$SRPBM - mean(data_train$SRPBM)) ^ 2)
		R2_tmp <- 1- SSE_tmp/SST_tmp
		PCC_tmp <- cor.test(Model_tmp$finalModel$predicted, data_train$SRPBM,method = "pearson")
		cat('>>> Pearson\'s r (PCC): \n')
		print(PCC_tmp)
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
	cat(paste('>>> Total RMSE: ', rmse_tmp, "\n", sep = ''))
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
	if (rmse_tmp < rmse_best) {
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


cat('>>> [2] Feature selsction ==> \n')
for (fn in seq(length(sort_fea_all)-1,2)) {feature_sel(fn)}
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
cat(paste('>>> Total R2: ', R2_best, "\n", sep = ''))
if (PM == 'rf') {
	cat('>>> Pearson\'s r (PCC): \n')
	print(PCC_best)
}
tune_met_best <- colnames(tune_best)
cat('>>> Best tune of final model: \n')
print(tune_best)


## write the observed and predicted expression to file
if (PM == 'rf') {
	obs_pred_exp <- data.frame(cbind(levels(data_train$Intron_pair), data_train$SRPBM, Model_best$finalModel$predicted))
	colnames(obs_pred_exp) <- c('Intron_pair','true_exp', 'pred_exp')
} else {
	best_pred <- Model_best$pred
	for (i in 1:length(tune_met_best)) {
		met_best <- tune_met_best[i]
		met_val_best <- tune_best[,met_best]
		best_pred <- best_pred[best_pred[,met_best]==met_val_best,]
		pred_best <- best_pred$pred
		obs_best <- best_pred$obs
		index_best <- best_pred$rowIndex
	}
	obs_pred_exp <- data.frame(cbind(levels(data_train[index_best,]$Intron_pair), obs_best, pred_best))
	colnames(obs_pred_exp) <- c('Intron_pair', 'true_exp', 'pred_exp')
}
write.table(obs_pred_exp, paste(cell, PM, 'train_pred_exp', sep = '_'), row.names = F, col.names = T, quote = F, sep = '\t')

# stop cluster and register sepuntial front end
stopCluster(cl); registerDoSEQ()

# save R data
rm(list = 'args')
save(list = objects(), file=paste(cell, PM, 'FS_exp.RData', sep = '_'))
cat('>>> Task DONE !!! \n')


### END
