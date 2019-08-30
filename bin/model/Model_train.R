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

# load data_train
file <- paste(cell, 'train', sep = '_')
data_train <- read.table(file, head = T)
data_train <- data.frame(data_train)
summary(data_train)
set.seed(seed)
cat('>>> Dimension of data matrix: \n')
dim(data_train)
n <- dim(data_train)[1]
FN <- dim(data_train)[2] - 5
col_name <- names(data_train)
fea_all <- col_name[5:(FN+4)]
data_train_mat <- data_train[,5:(FN+4)]
data_train$Type <- factor(as.character(data_train$Type))

# show which libraries were loaded  
cat('>>> Session Info: \n')
sessionInfo()

# register parallel front-end
cl <- makeCluster(cores); registerDoParallel(cl)

# Get feature importance
cat('>>> Start to sort importance ... ... \n')
ctrl <- trainControl(method = 'cv', number = 10 , savePredictions = "all", returnData = T, returnResamp = 'all', 
					 verboseIter = T,allowParallel = T)
Model_fea_all <- train(y = data_train$Type, x = data_train_mat, method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale"))
imp_all_fea <- varImp(Model_fea_all, scale = T)
sort_imp_all_fea <- data.frame(sortImp(imp_all_fea, FN))
cat('>>> Original importance (all feature): \n')
print(sort_imp_all_fea)

sort_fea_ori <- rownames(sort_imp_all_fea)
sort_fea_final <- sort_fea_ori[FN]
for (fn in rev(2:(FN-1)))
{
	cat('=====================================================\n')
	sort_fea_tmp <- sort_fea_ori[1:fn]
	cat(paste('>>> Top', fn, 'features (unsorted) are: \n', sep = ' '))
	print(sort_fea_tmp)
	assign(paste('Model_FN', fn, sep = ''), train(y = data_train$Type, x = data_train[sort_fea_tmp], method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale")))
	assign(paste('imp_FN', fn, sep = ''), varImp(get(paste('Model_FN', fn, sep = '')), scale = T))
	assign(paste('sort_imp_FN', fn, sep = ''), data.frame(sortImp(get(paste('imp_FN', fn, sep = '')), fn)))
	cat('>>> New ranking of features: \n')
	sort_fea_ori <<- c(rownames(get(paste('sort_imp_FN', fn, sep = ''))))
	print(sort_fea_ori)
	sort_fea_final <<- c(sort_fea_ori[fn], sort_fea_final)
	cat('=====================================================\n')
}
sort_fea_final <<- c(sort_fea_ori[1], sort_fea_final)
cat('>>> Task of sorting feature finished. \n')


for (fn in 1:FN)
{
	assign(paste('sort_fea', fn, sep = ''), sort_fea_final[1:fn])
	cat(paste('>>> Top', fn, 'sorted features are: \n', sep = ' '))
	print(get(paste('sort_fea', fn, sep = '')))
}

cat('>>> Write down feature ranking table ... ... \n')
cat(paste('    File name: ', 'Imp_', cell, '_', PM, "\n", sep = ''))
write.table(sort_fea_final, paste('Imp', cell, PM, sep = '_'), sep = '\n', col.names = F, row.names = F, quote = F)

# stop cluster and register sepuntial front end
stopCluster(cl); registerDoSEQ()

# save R data
rm(list = 'args')
save(list = objects(), file=paste(cell, PM, 'train.RData', sep = '_'))


### END
