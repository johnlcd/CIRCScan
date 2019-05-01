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
FIP_list <- args[5]
seed <- 111

# load data_train
file <- paste(cell, 'exp_train', sep = '_')
data_train_all <- read.table(file, head = T)
data_train_all <- data.frame(data_train_all)
data_train_raw <- data_train_all
data_train_all$SRPBM <- log2(data_train_all$SRPBM)
#data_train <- data_train_all[(2^(data_train_all$SRPBM)>0.01),] # remove low expression circRNAs (SRPBM > -6)
data_train <- data_train_all
if (PM == 'glm') {
	data_train$SRPBM <- data_train$SRPBM/10 # value = log2(SRPBM)/10
}
summary(data_train)
set.seed(seed)
print('>>> Dimension of data matrix: ')
dim(data_train)
n <- dim(data_train)[1]
FN <- dim(data_train)[2] - 5
col_name <- names(data_train)
fea_all <- col_name[5:(FN+4)]
if (args[4] != 'all'){
	fea_sel <- c(strsplit(args[4], ',')[[1]])
} else {
	fea_sel <- fea_all
}
data_train_mat <- data_train[,5:(FN+4)]
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

print('>>> [1] Train model with circRNAs expression (SRPBM) of ENCODE Long Non-Poly(A) RNA-seq data ... ... ')
ctrl <- trainControl(method = 'cv', number = 10 , savePredictions = "all", returnData = T, returnResamp = 'all', 
					 verboseIter = T, allowParallel = T)
Model <- train(y = data_train$SRPBM, 
			   x = data_train_mat, 
			   method = PM, trControl = ctrl, metric = "RMSE", importance = T, preProc = c("center", "scale"))

## write the observed and predicted expression to file
obs_pred_exp <- data.frame(cbind(levels(data_train$Intron_pair), data_train$SRPBM, Model$finalModel$predicted))
colnames(obs_pred_exp) <- c('Intron_pair','true_exp', 'pred_exp')
write.table(obs_pred_exp, paste(cell, PM, 'train_pred_exp', sep = '_'), row.names = F, col.names = T, quote = F, sep = '\t')

## model performance (RMSE, R2)
mse <- mse(data_train$SRPBM, Model$finalModel$predicted)
rmse <- sqrt(mse)
SSE <- sum((Model$finalModel$predicted - data_train$SRPBM)^2)
SST <- sum((data_train$SRPBM - mean(data_train$SRPBM)) ^ 2)
R2 <- 1- SSE/SST
PCC <- cor.test(Model$finalModel$predicted, data_train$SRPBM, method = "pearson")
print('    Model performance ==> ')
print(paste('    Total RMSE: ', rmse, sep = ''))
print(paste('    Total R2: ', R2, sep = ''))
print('    Pearson\'s r (PCC): ')
print(PCC)

## feature importance
if (PM == 'rf') {
	imp <- importance(Model$finalModel)
	write.table(imp, paste(cell, PM, 'Imp_all', sep = '_'), row.names = T, col.names = T, quote = F, sep = '\t')
	print('    Importance ==> ')
	print(imp)
	sort_imp <- sort(importance(Model$finalModel)[,"%IncMSE"], decreasing = T)
	sort_imp <- data.frame(sort_imp)
	colnames(sort_imp) <- c('Importance')
	write.table(sort_imp, paste(cell, PM, 'sort_Imp', sep = '_'), row.names = T, col.names = F, quote = F, sep = '\t')
} else {
	imp <- varImp(Model)
	imp <- round(imp$importance, 2)
	imp <- data.frame(imp)
	order_imp <- order(imp[,'Overall'], decreasing = T)
	sort_fea_all <- rownames(imp)[order_imp]
	sort_imp <- data.frame(imp[,'Overall'][order_imp])
	rownames(sort_imp) <- sort_fea_all
	colnames(sort_imp) <- 'Importance'
	sort_imp <- data.frame(sort_imp)
	colnames(sort_imp) <- c('Importance')
}
print('    Sorted importance of MSE ==> ')
print(sort_imp)
if (PM == "rf") {
	pdf(file = paste(cell, PM, "Imp.pdf", sep = "_"))
	varImpPlot(Model$finalModel, type = 1, main = 'Feature Importance')
	dev.off()
}
print('    Summary of model: ')
print(Model)

# stop cluster and register sepuntial front end
stopCluster(cl); registerDoSEQ()

# load pred data
print('>>> [2] Predict circRNAs expression (SRPBM) with final model ... ... ')
data_pred <- read.table(paste(cell, 'exp_pred', sep = '_'), head = T)
data_pred.bed <- data_pred[1:4]
all_pred_mat <- data_pred[-(1:3)]
sel_pred_mat <- all_pred_mat[,fea_sel]
rownames(sel_pred_mat) <- data_pred$Intron_pair
print('>>> Number of intron pairs and all feature:')
print(paste(dim(all_pred_mat)[1], dim(all_pred_mat)[2]-2, sep = ';'))

# prediction by trained model
circ_pred <- predict(Model, sel_pred_mat)
circ_pred <- 2^circ_pred
circ_pred_df <- data.frame(circ_pred)
circ_pred_bed <- data.frame(cbind(data_pred.bed, circ_pred_df))
colnames(circ_pred_bed) <- c('Chr','Start','End','Intron_pair','SRPBM')
max_exp <- max(circ_pred_df)
min_exp <- min(circ_pred_df)
max_exp_bed <- circ_pred_bed[c(which(circ_pred_df == max_exp)),]
min_exp_bed <- circ_pred_bed[c(which(circ_pred_df == min_exp)),]
print('    The max expressed circRNA:')
print(max_exp_bed)
print('    The min expressed circRNA:')
print(min_exp_bed)

write.table(circ_pred_bed, paste(cell, PM, 'pred_exp.bed', sep = '_'), row.names = F, col.names = F, quote = F, sep = '\t')

# save R data
save(list = objects(), file=paste(cell, PM, 'train_pred_exp.RData', sep = '_'))

# validate and evaluate the results by known circRNA data
## read known circRNA FIP list
print('>>> [3] Validate and evaluate the expression level of circRNAs ... ... ')
CB_fip_list <- read.table(FIP_list, col.names='Intron_pair')
CB_fip_list <- CB_fip_list$Intron_pair
train_fip_list <- data_train_all$Intron_pair
report_fip_list <- union(train_fip_list, CB_fip_list)
col_name <- c('Chr', 'Start', 'End', 'Intron_pair', 'SRPBM')
all_exp_df <- data.frame(rbind(data_train_raw[,col_name], circ_pred_bed))
all_fip_list <- all_exp_df$Intron_pair
diff_fip_list <- setdiff(all_fip_list, report_fip_list) # un-reported circRNAs list
print(paste('    Total of', length(CB_fip_list), 'known circRNAs (circBase) FIPs among all', length(all_fip_list), 'FIPs.', sep = ' '))
CB_exp <- all_exp_df[all_exp_df$Intron_pair%in%CB_fip_list,]
diff_exp <- all_exp_df[all_exp_df$Intron_pair%in%diff_fip_list,]
write.table(CB_exp, paste(cell, PM, 'circbase_pred_exp', sep = '_'), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(diff_exp, paste(cell, PM, 'unreported_pred_exp', sep = '_'), row.names = F, col.names = F, quote = F, sep = '\t')


## T-test of circBase expression level and others value
print('>>> [4] T-test of circRNAs expression level between reported and un-reported groups ... ... ')
print(paste('    Mean expression level of all genome:', mean(all_exp_df$SRPBM), sep = ' '))
print(paste('    Mean value of reported and un-reported groups:', mean(CB_exp$SRPBM), ";", mean(diff_exp$SRPBM), sep = ' '))
CB_SRPBM <- CB_exp$SRPBM
diff_SRPBM <- diff_exp$SRPBM
tt_SRPBM <- t.test(CB_SRPBM, diff_SRPBM, paired = F)
print('    Summary of T-test: ')
print(tt_SRPBM)
print('    P-value: ')
print(tt_SRPBM$p.value)
if (tt_SRPBM$p.value < 0.05) {
	print('    Significant !!!')
} else {
	print('    Not significant.')
}

print('>>> Task DONE !!!')

# save R data
rm(list = 'args')
save(list = objects(), file=paste(cell, PM, 'train_pred_exp.RData', sep = '_'))


### END
