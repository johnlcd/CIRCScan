# please make sure required packages already installed


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
cores <- as.numeric(args[3])
seed <- as.numeric(args[4])

# load train data
data_train <- read.table(paste(cell, 'train', sep = '_'), head = T)
summary(data_train)
feature.names=names(data_train)
data_train$Type <- as.factor(data_train$Type)
print('>>> Dimension of data matrix: ')
dim(data_train)
n <- dim(data_train)[1]
FN <- dim(data_train)[2] - 5
sel_out <- 1:n
for (i in 1:4)
{
	assign(paste('indextest', i, sep = ''), sort(sample(sel_out, n*0.2)))
	assign(paste('indextrain', i, sep = ''), setdiff(1:n, get(paste('indextest', i, sep = ''))))
	sel_out <- setdiff(sel_out, get(paste('indextest', i, sep = '')))
	assign(paste('train', i, sep = ''), data_train[get(paste('indextrain', i, sep = '')), ])
	assign(paste('test', i, sep = ''), data_train[get(paste('indextest', i, sep = '')), ])
}
indextest5 <- sel_out
indextrain5 <- setdiff(1:n, indextest5)
train5 <- data_train[indextrain5, ]
test5 <- data_train[indextest5, ]

# show which libraries were loaded  
sessionInfo()

# register parallel front-end
cl <- makeCluster(cores); registerDoParallel(cl)

# train model
ctrl <- trainControl(method = 'cv', number = 10, allowParallel = T)

F1_best <- 0.5
gp_best <- 0
for (gp in 1:5)
{

	print('############################################')
	print(paste('>>> Group', gp, ':'))

# load X and Y (this will be transferred to to train function)
	train <- get(paste('train', gp, sep = ''))
	test <- get(paste('test', gp, sep = ''))
	X = get(paste('train', gp, sep = ''))[,5:(FN+4)]
	str(X)
	Y = as.factor(get(paste('train', gp, sep = ''))$Type)
	assign(paste('test_mat_gp', gp, sep = ''), get(paste('test', gp, sep = ''))[,5:(FN+4)])
#	assign(get(paste('test', gp, sep = ''))$Type, as.factor(get(paste('test', gp, sep = ''))$Type))
				
# train model
	assign(paste('Model_gp', gp, sep = ''), train(y = Y, x = X, method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale")))
	get(paste('Model_gp', gp, sep = ''))

## estimate variable importance  
	assign(paste('importance_gp', gp, sep = ''), varImp(get(paste('Model_gp', gp, sep = '')), scale = FALSE))

## summarize importance  
	print('>>> Importance:')
	print(get(paste('importance_gp', gp, sep = '')))

	assign(paste('Pred_gp', gp, sep = ''), predict(get(paste('Model_gp', gp, sep = '')), get(paste('test_mat_gp', gp, sep = ''))))
	assign(paste('Prob_gp', gp, sep = ''), predict(get(paste('Model_gp', gp, sep = '')), get(paste('test_mat_gp', gp, sep = '')), type = 'prob'))
	assign(paste('Pred_prob_gp', gp, sep = ''), get(paste('Prob_gp', gp, sep = ''))[, 2])
	assign(paste('Prediction_gp', gp, sep = ''), prediction(predictions = get(paste('Pred_prob_gp', gp, sep = '')), labels = as.factor(get(paste('test', gp, sep = ''))$Type)))
	assign(paste('Perf.roc_gp', gp, sep = ''), performance(get(paste('Prediction_gp', gp, sep = '')), measure = 'tpr', x.measure = 'fpr'))
	assign(paste('Perf.auc_gp', gp, sep = ''), performance(get(paste('Prediction_gp', gp, sep = '')), measure = 'auc'))
	if (cell != 'GM12878'){
		print('>>> Confusion matrix:')
		print(confusionMatrix(get(paste('Pred_gp', gp, sep = '')), test$Type, positive = 'TRUE'))
	}
#	if (PM == 'svmRadial'){
#		if (cell %in% c('K562', "H1-hESC", "HeLa-S3", "K562")){
#			print('>>> Confusion matrix:')
#			print(confusionMatrix(get(paste('Pred_gp', gp, sep = '')), test$Type, positive = 'TRUE'))
#		}
#	}

	assign(paste('Precision_gp', gp, sep = ''), posPredValue(get(paste('Pred_gp', gp, sep = '')), test$Type, positive = 'TRUE'))
	assign(paste('Recall_gp', gp, sep = ''), sensitivity(get(paste('Pred_gp', gp, sep = '')), test$Type, positive = 'TRUE'))
	assign(paste('Spe_gp', gp, sep = ''), specificity(get(paste('Pred_gp', gp, sep = '')), test$Type, negative = 'FALSE'))
	assign(paste('AUC_gp', gp, sep = ''), unlist(get(paste('Perf.auc_gp', gp, sep = ''))@y.values))
	assign(paste('F1_gp', gp, sep = ''), (2 * get(paste('Precision_gp', gp, sep = '')) * get(paste('Recall_gp', gp, sep = ''))) / (get(paste('Precision_gp', gp, sep = '')) + get(paste('Recall_gp', gp, sep = ''))))
	print('>>> Precision is:')
	print(get(paste('Precision_gp', gp, sep = '')))
	print('>>> Sensitivity is:')
	print(get(paste('Recall_gp', gp, sep = '')))
	print('>>> Specificity is:')
	print(get(paste('Spe_gp', gp, sep = '')))
	print('>>> AUC is:')
	print(get(paste('AUC_gp', gp, sep = '')))
	print('>>> F1 score is:')
	print(get(paste('F1_gp', gp, sep = '')))

	print(paste('>>> Model training of group', gp, 'finished.'))
	print('############################################')

	if (get(paste('F1_gp', gp, sep = '')) > F1_best){
		gp_best <- gp
		Pred_best <- get(paste('Pred_gp', gp, sep = ''))
		Prob_best <- get(paste('Prob_gp', gp, sep = ''))
		Pred_prob_best <- get(paste('Pred_prob_gp', gp, sep = ''))
		Prediction_best <- get(paste('Prediction_gp', gp, sep = ''))
		Perf.roc_best <- get(paste('Perf.roc_gp', gp, sep = ''))
		Perf.auc_best <- get(paste('Perf.auc_gp', gp, sep = ''))
		F1_best <- get(paste('F1_gp', gp, sep = ''))
		Model_best <- get(paste('Model_gp', gp, sep = ''))
		Precision_best <- get(paste('Precision_gp', gp, sep = ''))
		Recall_best <- get(paste('Recall_gp', gp, sep = ''))
		Spe_best <- get(paste('Spe_gp', gp, sep = ''))
		AUC_best <- get(paste('AUC_gp', gp, sep = ''))
	}

}

print('+++++++++++++++++++++++++++++++++++++++++++++++++++++')
if (exists('Model_best')){
	print('>>> Best group:')
	print(paste('Group', gp_best))
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


# register parallel front-end
cl <- makeCluster(cores); registerDoParallel(cl)

# stop cluster and register sequntial front end
stopCluster(cl); registerDoSEQ()

# save R data
save(list = objects(), file=paste(cell, PM, 'train_2fea.RData', sep = '_'))


### END
