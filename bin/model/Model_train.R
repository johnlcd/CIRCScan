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
print('>>> Dimension of data matrix: ')
dim(data_train)
n <- dim(data_train)[1]
FN <- dim(data_train)[2] - 5
col_name <- names(data_train)
fea_all <- col_name[5:(FN+4)]
data_train_mat <- data_train[,5:(FN+4)]
data_train$Type <- factor(as.character(data_train$Type))
sel_out <- 1:n
for (i in 1:4)
{
	assign(paste('indextest', i, sep = ''), sort(sample(sel_out, n*0.2)))
	assign(paste('indextrain', i, sep = ''), setdiff(1:n, get(paste('indextest', i, sep = ''))))
	sel_out <- setdiff(sel_out, get(paste('indextest', i, sep = '')))
	assign(paste('train', i, sep = ''), data.frame(data_train[get(paste('indextrain', i, sep = '')), ]))
	assign(paste('test', i, sep = ''), data.frame(data_train[get(paste('indextest', i, sep = '')), ]))
}
indextest5 <- sel_out
indextrain5 <- setdiff(1:n, indextest5)
train5 <- data.frame(data_train[indextrain5, ])
test5 <- data.frame(data_train[indextest5, ])

# show which libraries were loaded  
print('>>> Session Info: ')
sessionInfo()

# register parallel front-end
cl <- makeCluster(cores); registerDoParallel(cl)

# Get feature importance
print('>>> Start to sort importance ... ...')
ctrl <- trainControl(method = 'cv', number = 10 , allowParallel = T)
Model_fea_all <- train(y = data_train$Type, x = data_train_mat, method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale"))
imp_all_fea <- varImp(Model_fea_all, scale = T)
sort_imp_all_fea <- data.frame(sortImp(imp_all_fea, FN))
print('>>> Original importance (all feature):')
print(sort_imp_all_fea)

sort_fea_ori <- rownames(sort_imp_all_fea)
sort_fea_final <- sort_fea_ori[FN]
for (fn in rev(2:(FN-1)))
{
	print('=====================================================')
	sort_fea_tmp <- sort_fea_ori[1:fn]
	print(paste('>>> Top', fn, 'features (unsorted) are: ', sep = ' '))
	print(sort_fea_tmp)
	assign(paste('Model_FN', fn, sep = ''), train(y = data_train$Type, x = data_train[sort_fea_tmp], method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale")))
	assign(paste('imp_FN', fn, sep = ''), varImp(get(paste('Model_FN', fn, sep = '')), scale = T))
	assign(paste('sort_imp_FN', fn, sep = ''), data.frame(sortImp(get(paste('imp_FN', fn, sep = '')), fn)))
	print('>>> New ranking of features:')
	sort_fea_ori <<- c(rownames(get(paste('sort_imp_FN', fn, sep = ''))))
	print(sort_fea_ori)
	sort_fea_final <<- c(sort_fea_ori[fn], sort_fea_final)
	print('=====================================================')
}
sort_fea_final <<- c(sort_fea_ori[1], sort_fea_final)
print('>>> Task of sorting feature finished.')


for (fn in 1:FN)
{
	assign(paste('sort_fea', fn, sep = ''), sort_fea_final[1:fn])
	print(paste('>>> Top', fn, 'sorted features are: ', sep = ' '))
	print(get(paste('sort_fea', fn, sep = '')))
}

print('>>> Write down feature ranking table ... ...')
print(paste('    File name: ', 'Imp_', cell, '_', PM, sep = ''))
write.table(sort_fea_final, paste('Imp', cell, PM, sep = '_'), sep = '\n', col.names = F, row.names = F, quote = F)

# stop cluster and register sepuntial front end
stopCluster(cl); registerDoSEQ()

# save R data
save(list = objects(), file=paste(cell, PM, 'train.RData', sep = '_'))


### END
