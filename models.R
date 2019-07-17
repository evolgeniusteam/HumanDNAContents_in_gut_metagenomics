## ---- create data.list ---- ##
library(SIAMCAT)
library(dplyr)
library(pROC)
library(doParallel)

load("data.Rdata")
load("sigfeatures.Rdata")
crc.meta <- read.delim("allCRC_metadata.txt", header = T, sep = "\t",as.is = T)
rownames(crc.meta) <- crc.meta$Sample_ID

repeat.rf.func <- function(feat.data, meta.data,num.folds = 10,num.resample = 10){
  ##it is not necessary to normalize feat.data when using randomforest model
  ##feat.data is relative abundance (sums = 1)
  rownames(meta.data) <- meta.data$Sample_ID
  siamcat <- siamcat(feat = feat.data, ## relative abundance, features as rows and samples as cols
                     meta = meta.data, ## metadata, rownames are sampleID
                     label = "Group", case = "CRC")
  siamcat.split <- create.data.split(siamcat, num.folds = num.folds,
                                     num.resample = num.resample)
  siamcat.train <- train.model(siamcat = siamcat.split,
                               method = "randomForest",
                               perform.fs = T, feature.type = "original",
                               param.fs = list(thres.fs = nrow(feat.data)-1,
                                               method.fs = "AUC"))
  siamcat.pred <- make.predictions(siamcat.train)
  ##importance scores in each model
  siamcat.weights.matrix <- weight_matrix(siamcat.pred)
  ##prediction scores
  siamcat.train.pred.label <- label(siamcat.pred)
  siamcat.train.pred.matrix <- pred_matrix(siamcat.pred)
  siamcat.train.pred.matrix <- siamcat.train.pred.matrix[names(siamcat.train.pred.label$label),]
  summ.stat = "mean"
  roc.mean <- roc(response = siamcat.train.pred.label$label,
                  predictor = apply(siamcat.train.pred.matrix, 1, summ.stat),
                  ci =T)
  return(list(models = siamcat.pred,
              features.imp = siamcat.weights.matrix,
              roc.mean = roc.mean))
}

cross.models <- function(project, metadata = metadata, rl.data = rl.data){
  meta.train <- metadata %>% filter(Project == project)
  feat.train <- rl.data[,meta.train %>% pull(Sample_ID)]
  meta.train <- data.frame(meta.train)
  rownames(meta.train) <- meta.train$Sample_ID
  ##train model
  models <- repeat.rf.func(feat.train, meta.train)
  return(models)
}

lodo.models <- function(project, metadata = metadata, rl.data = rl.data){
  meta.train <- metadata %>% filter(Project != project)
  feat.train <- rl.data[,meta.train %>% pull(Sample_ID)]
  meta.train <- data.frame(meta.train)
  rownames(meta.train) <- meta.train$Sample_ID
  ##train model
  models <- repeat.rf.func(feat.train, meta.train)
  return(models)
}

## ----- it would spend several days on modelling ----- ##
## ----- advise: Parallel Computing -----##
## ----- example -----##
cl <- makeCluster(7)
registerDoParallel(cl)
study <- unique(crc.meta$Project)
all.spe.cross <- foreach(project = study) %dopar% cross.models(project = project, metadata=crc.meta, rl.data = spe.data.df)
all.spe.lodo <- foreach(project = study) %dopar% lodo.models(project = project, metadata=crc.meta, rl.data = spe.data.df)
stopCluster(cl)

## ----- example ---------##
cl <- makeCluster(7)
registerDoParallel(cl)
study <- unique(crc.meta$Project)
spe.wt.data.df <- spe.data.df[spe.wt.list.spe.vector,]
wt.spe.cross <- foreach(project = study) %dopar% cross.models(project = project, metadata=crc.meta, rl.data = spe.wt.data.df )
wt.spe.lodo <- foreach(project = study) %dopar% lodo.models(project = project, metadata=crc.meta, rl.data = spe.wt.data.df )
stopCluster(cl)
