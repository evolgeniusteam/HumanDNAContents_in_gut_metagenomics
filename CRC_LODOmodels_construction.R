## --------------Create LODO models------------------##

args <- commandArgs(T)
library(ROCR)
library(dplyr)
library(randomForest)
library(caret)
library(purrr)
library(pROC)

load(args[1])
## LODO ----------
study <- unique(crc.meta$Project)

## ----create ten-time ten-fold models ----
## ----in total 100models -----------------
repeated.rf.confusion <- function(a, times = 10, k = 10){
  ##require group name is "Group"
  sample.split <- createMultiFolds(a$Group, times = times, k=k)
  a$Group <- as.factor(a$Group)
  rf.models <- list()
  pred.matrix <- list()
  confusion.result <- list()
  for (iter in names(sample.split)) {
    train.data <- a[sample.split[[iter]], ]
    test.data <- a[-sample.split[[iter]], ]
    rf.models[[iter]] <- randomForest(Group ~ .,train.data, proximity = T, importance = T, ntree = 500)
    predictions <- predict(rf.models[[iter]], test.data, type = "prob")
    pred.matrix[[iter]] <- predictions
    
    predictions.class <- predict(rf.models[[iter]], test.data, type= "response")
    confusion.result[[iter]] <- predictions.class
  }
  pred.matrix.list <- list()
  pred.df <- data.frame()
  pred.df.2 <- data.frame()
  
  fold.str <- paste("Fold", gsub(" ", "0", format(1:k)), sep = "")
  rep.str <- paste("Rep", gsub(" ", "0", format(1:times)), sep = "")
  for (rep in rep.str) {
    for (fold in fold.str) {
      iter <- paste(fold, rep, sep = ".")
      pred.df.2 <- data.frame(pred.matrix[[iter]], stringsAsFactors = F)
      pred.df <- rbind(pred.df, pred.df.2)
    }
    for (grp in names(pred.df)) {
      y <- data.frame(pred.df[,grp], stringsAsFactors = F);
      names(y) <- grp;
      y$Sample <- rownames(pred.df)
      pred.matrix.list[[grp]][[rep]] <- y
    }
    pred.df <- data.frame();
  }
  ## --- importance scores --- ##
  importance <- lapply(rf.models, function(data){ imp <- importance(data, type=1)})
  importance.df <- reduce(importance, cbind)
  importance.df <- data.frame(importance.df, stringsAsFactors = F)
  names(importance.df) <- names(rf.models)
  return(list(sample.split = sample.split,
              rf.models = rf.models,
              pred.matrix = pred.matrix,
              confusion.result = confusion.result,
              importance.df = importance.df,
              pred.matrix.list = pred.matrix.list
              #roc.df = roc.df,
              #roc.plot= roc.plot,
              #forest.auc = forest.auc,
  ))
}

## ----create LODO models -----------------
rf.funct2 <- function(pro, crc.meta = crc.meta, data.crc = data.crc){
  library(ROCR)
  library(dplyr)
  library(randomForest)
  library(caret)
  library(purrr)
  library(pROC)
  test.sample <- subset(crc.meta, Project %in% pro)$Sample_ID
  train.data <- data.crc[!rownames(data.crc) %in% test.sample, ]
  
  ##train model
  models <- repeated.rf.confusion(train.data,
                                  times = 10, k = 10)
  return(models)
}


## ----function for LODO prediction  -------------------##

lodo.pred.func <- function(rf.models = data$rf.models, metadata = crc.meta, test.data = test.data){
  pred.result <- list()
  
  for (iter in names(rf.models)) {
    predictions <- predict(rf.models[[iter]], test.data, type = "prob")
    predictions <- data.frame(predictions, stringsAsFactors = F)
    predictions$Sample_ID <- rownames(predictions)
    pred.result[[iter]] <- predictions
  }
  pred.result.df <- pred.result %>%
    reduce(rbind) %>% group_by(Sample_ID) %>%
    summarise(control_final = mean(CTR),
              case_final = mean(CRC))
  pred.result.df <- data.frame(pred.result.df, 
                               stringsAsFactors = F)
  
  rownames(pred.result.df) <- pred.result.df$Sample_ID
  
  names(pred.result.df)[2:3] <- c("CTR","CRC")
  pred.result.df$max <- apply(pred.result.df, 1, 
                              function(data){y <- names(which.max(data))})
  pred.result.df$max <- factor(pred.result.df$max,
                               levels = c("CTR","CRC"))
  
  pred.result.df$true <- test.data[rownames(pred.result.df), 
                                   "Group"]
  pred.result.df$true <- factor(pred.result.df$true,
                                levels = c("CTR","CRC"))
  
  forestpred <- ROCR::prediction(pred.result.df$CRC, 
                                 pred.result.df$true, 
                                 label.ordering = c("CTR","CRC"));
  forestperf <- performance(forestpred, "tpr", "fpr");
  forest.auc <- performance(forestpred, "auc")@y.values[[1]]
  
  return(list(pred.result.df =pred.result.df,
              forestpred = forestpred,
              forestperf = forestperf,
              forest.auc = forest.auc))
}

study <- unique(crc.meta$Project)

library(doParallel)
cl <- makeCluster(7)
registerDoParallel(cl)

spe.lodo.models <- list()
spe.lodo.pred <- list()

## get LODO models ----------##
data.crc <- spe.list[[1]]
a <- foreach(project = study) %dopar% rf.funct2(pro = project, crc.meta = crc.meta,  data.crc = data.crc)
spe.lodo.models[[1]] <- a

data.crc <- spe.list[[2]]
a <- foreach(project = study) %dopar% rf.funct2(pro = project, crc.meta = crc.meta,  data.crc = data.crc)
spe.lodo.models[[2]] <- a

data.crc <- spe.list[[3]]
a <- foreach(project = study) %dopar% rf.funct2(pro = project, crc.meta = crc.meta,  data.crc = data.crc)
spe.lodo.models[[3]] <- a



spe.lodo.nohdc.models <- list()
spe.lodo.nohdc.pred <- list()

spe.nohdc.list <- lapply(spe.list, function(data){
                           y <- data[, !names(data) %in% "HDC"];
                           return(y)
                           })
data.crc <- spe.nohdc.list[[1]]
a <- foreach(project = study) %dopar% rf.funct2(pro = project, crc.meta = crc.meta,  data.crc = data.crc)
spe.lodo.nohdc.models[[1]] <- a

data.crc <- spe.nohdc.list[[2]]
a <- foreach(project = study) %dopar% rf.funct2(pro = project, crc.meta = crc.meta,  data.crc = data.crc)
spe.lodo.nohdc.models[[2]] <- a

data.crc <- spe.nohdc.list[[3]]
a <- foreach(project = study) %dopar% rf.funct2(pro = project, crc.meta = crc.meta,  data.crc = data.crc)
spe.lodo.nohdc.models[[3]] <- a

##
save.image(file = paste0(args[2], "_LODO_models.RData"))
##

