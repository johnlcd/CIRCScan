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
summary(data_train)
feature_names <- names(data_train)
data_train$Type <- as.factor(data_train$Type)
set.seed(seed)
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
print('>>> Session Info: ')
sessionInfo()

# register parallel front-end
cl <- makeCluster(cores); registerDoParallel(cl)

# Get feature importance
for (j in 1:5)
{
	print('############################################')
	print(paste('>>> Group', j, ':'))

# load X and Y (this will be transferred to to train function)
	train <- get(paste('train', j, sep = ''))
	test <- get(paste('test', j, sep = ''))
	X = train[,5:(FN+4)]
	str(X)
	Y = as.factor(train$Type)
	
# train model
	ctrl <- trainControl(method = 'cv', number = 10, allowParallel = T)
	assign(paste('Model_gp', j, sep = ''), train(y = Y, x = X, method = PM, trControl = ctrl, prob.model = TRUE, preProc = c("center", "scale")))
	get(paste('Model_gp', j, sep = ''))
	
## estimate variable importance  
	assign(paste('importance_gp', j, sep = ''), varImp(get(paste('Model_gp', j, sep = '')), scale = FALSE))

## summarize importance  
	print('>>> Importance:')
	print(get(paste('importance_gp', j, sep = '')))
	assign(paste('sort_imp_gp', j, sep = ''), as.data.frame(sortImp(get(paste('importance_gp', j, sep = '')), FN)))
	assign(paste('top_imp_gp', j, sep = ''), as.data.frame(sortImp(get(paste('importance_gp', j, sep = '')), 10)))
	print('>>> Sorted importance:')
	print(get(paste('sort_imp_gp', j, sep = '')))
	print('>>> Top10 importance:')
	print(get(paste('top_imp_gp', j, sep = '')))
	write.table(get(paste('sort_imp_gp', j, sep = ''))[1], paste('Imp', cell, PM, j, sep = '_'), sep = '\t', col.names = F, row.names = T, quote = F)

## plot importance  
#	pdf(file = paste(cell, '_', PM, '_gp', j, '_imp.pdf', sep = ''))
	
#	plot(get(paste('importance_gp', j, sep = '')), main = paste('Importance of model ', PM, ' (group ', j, ')', sep = ''), xlab = 'Imporance', ylab = 'Feature')
#	plot(get(paste('importance_gp', j, sep = '')))

#	dev.off()
	
	print(paste('>>> group', j, 'finished.'))
	print('############################################')
}

# stop cluster and register sequntial front end
stopCluster(cl); registerDoSEQ()

# save R data
save(list = objects(), file=paste(cell, PM, 'train.RData', sep = '_'))


### END
