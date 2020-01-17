
## -----CD external validation -----##
## -----taxonomic level ------------##
prjna384246.abund <- read.delim("prjna384246_meta.txt",
                                header = T, sep = "\t", as.is = T)
prjna384246.abund.result <- data.phylum.func(prjna384246.abund, i = 10)
prjna384246.abund.s <- prjna384246.abund.result$data.abund.s
rownames(prjna384246.abund.s) <- prjna384246.abund.s$Species
prjna384246.abund.s <- prjna384246.abund.s[,-1]/100

prjna384246.meta <- read.delim("PRJNA384246_runinfo",
                                header = T, sep = "\t", as.is = T)
prjna384246.meta <- prjna384246.meta[,c("Library_Name", 
                                        "Run", "Sample_Name",
                                        "Week_0", "Week_14", "age", "age_diagnosis", 
                                        "colitis", "current_immunomodulators",
                                        "current_steroids","duration_disease","female","host_disease",
                                        "prior_antitnf","prism_id","remission_week14",
                                        "remission_week30", "remission_week54")]
prjna384246.meta$Disease <- substr(prjna384246.meta$host_disease, 1, 2)

## ----HDC
prjna384246.reads <- read.delim("prjna384246_reads.txt",
                                header = F, sep = " ", as.is = T)
reads.func <- function(reads.data, i){
  ##original.abund.data is the raw abundant data from metaphlan2
  ##reads.data: one row is filename and the reads number*4
  ##all.metadata which contains information of a project
  ##i is the length of runID
  ##**dplyr isnt compatible with plyr
  names(reads.data) <- c("filenames","reads")
  reads.data$sample <- substr(reads.data$filenames, 1,i)
  print(nrow(reads.data))
  reads.data$reads <- as.numeric(reads.data$reads)
  reads.data$reads <- reads.data$reads/4
  reads.data$type <- substr(reads.data$filenames, i+14, i+18)
  reads.data.all <- reads.data %>% group_by(sample, type) %>% summarise(all = sum(reads))
  reads.data.all.dcast <- dcast(reads.data.all, sample ~ type)
  reads.data.all.dcast$all <- reads.data.all.dcast$clean+reads.data.all.dcast$conta
  reads.data.all.dcast$humanper <- reads.data.all.dcast$conta/reads.data.all.dcast$all
  reads.data.all.dcast$HDC <- reads.data.all.dcast$humanper *100
  return(list(reads.data= reads.data,
              reads.data.all.dcast = reads.data.all.dcast))
}

prjna384246.reads.result <- reads.func(prjna384246.reads, 
                                       i = 10)
prjna384246.HDC <- prjna384246.reads.result$reads.data.all.dcast
prjna384246.meta.merge <- merge(prjna384246.meta, 
                                prjna384246.HDC, 
                                by.x = "Run", 
                                by.y = "sample", 
                                all.x = T)

prjna384246.meta.cd <- subset(prjna384246.meta.merge, 
                              Disease %in% "CD")

## -----select CD patients data in week14
prjna384246.meta.cd.week14 <- subset(prjna384246.meta.cd, 
                                     !is.na(Week_14))
prjna384246.sample <- subset(prjna384246.meta.cd, 
                             !is.na(Week_14)) %>% 
  pull(Run)

## ------create group about treatment response
prjna384246.meta.cd.week14$Group <- "non_responde"
prjna384246.meta.cd.week14[which(prjna384246.meta.cd.week14$remission_week14 %in% "1"), "Group"] <- "responde"
rownames(prjna384246.meta.cd.week14) <- prjna384246.meta.cd.week14$Run

prjna384246.abund.s.t <- data.frame(t(prjna384246.abund.s),
                                    stringsAsFactors = F)
prjna384246.abund.s.t.match <- prjna384246.abund.s.t[, names(prjna384246.abund.s.t) %in% names(SRP057.abund.s.t)[1:223]]

##------fill lack data--------------##
not.match <- names(SRP057.abund.s.t)[1:223][!names(SRP057.abund.s.t)[1:223] %in% names(prjna384246.abund.s.t)]
for (i in not.match) {
  prjna384246.abund.s.t.match[,i] <- 0
} 

prjna384246.abund.s.final <- prjna384246.abund.s.t.match[rownames(prjna384246.abund.s.t.match) %in% prjna384246.sample,]
prjna384246.abund.s.final[,c("HDC", "Group")] <- prjna384246.meta.cd.week14[rownames(prjna384246.abund.s.final),c("HDC","Group")]

prjna384246.abund.s.list <- list(all.spe = prjna384246.abund.s.final,
                                 dif.spe = prjna384246.abund.s.final[, names(prjna384246.abund.s.final) %in% c(srp057.dif.spe, "Group", "HDC")],
                                 sig.spe = prjna384246.abund.s.final[, names(prjna384246.abund.s.final) %in% c(srp057.sig.spe, "Group", "HDC")], 
                                 all.spe.nohdc = prjna384246.abund.s.final[, !names(prjna384246.abund.s.final) %in% "HDC"],
                                 dif.spe.nohdc = prjna384246.abund.s.final[, names(prjna384246.abund.s.final) %in% c(srp057.dif.spe, "Group")],
                                 sig.spe.nohdc = prjna384246.abund.s.final[, names(prjna384246.abund.s.final) %in% c(srp057.sig.spe, "Group")])

prjna384246.pred.list <- list()
out.valid.pred <- function(rf.models, 
                           test.data = test.data){
  ## ---for external validation -----------##
  pred.result <- list()
  
  for (iter in names(rf.models)) {
    predictions <- predict(rf.models[[iter]], test.data, type = "prob")
    predictions <- data.frame(predictions, stringsAsFactors = F)
    predictions$Sample_ID <- rownames(predictions)
    pred.result[[iter]] <- predictions
  }
  pred.result.df <- pred.result %>%
    reduce(rbind) %>% group_by(Sample_ID)
  pred.result.df <- data.frame(pred.result.df, 
                               stringsAsFactors = F)
}
for (i in names(prjna384246.abund.s.list)) {
  prjna384246.pred.list[[i]] <- out.valid.pred(srp057.treat.rf.result[[i]]$rf.models, 
                                               test.data  = prjna384246.abund.s.list[[i]])
}


prjna384246.pred.df.list <- list()
prjna384246.pred.result.list <- list()

for (i in names(prjna384246.pred.list)) {
  prjna384246.pred.df.list[[i]] <- prjna384246.pred.list[[i]] %>%
    group_by(Sample_ID) %>%
    summarise(non = mean(non_responde),
              yes = mean(responde)) 
  prjna384246.pred.df.list[[i]] <- data.frame(prjna384246.pred.df.list[[i]], 
                               stringsAsFactors = F)
  
  rownames(prjna384246.pred.df.list[[i]]) <- prjna384246.pred.df.list[[i]]$Sample_ID
  
  names(prjna384246.pred.df.list[[i]])[2:3] <- c("non_responde","responde")
  prjna384246.pred.df.list[[i]]$max <- apply(prjna384246.pred.df.list[[i]], 1, 
                              function(data){y <- names(which.max(data))})
  prjna384246.pred.df.list[[i]]$max <- factor(prjna384246.pred.df.list[[i]]$max,
                               levels = c("non_responde","responde"))
  
  prjna384246.pred.df.list[[i]]$true <- prjna384246.abund.s.list[[i]][rownames(prjna384246.pred.df.list[[i]]), 
                                   "Group"]
  prjna384246.pred.df.list[[i]]$true <- factor(prjna384246.pred.df.list[[i]]$true,
                                levels = c("non_responde","responde"))
  
  forestpred <- ROCR::prediction( prjna384246.pred.df.list[[i]][,"responde"], 
                                  prjna384246.pred.df.list[[i]]$true, 
                                 label.ordering = c("non_responde","responde"));
  forestperf <- performance(forestpred, "tpr", "fpr");
  forest.auc <- performance(forestpred, "auc")@y.values[[1]]
  prjna384246.pred.result.list[[i]]$forestpred <- forestpred
  prjna384246.pred.result.list[[i]]$forestperf <- forestperf
  prjna384246.pred.result.list[[i]]$forest.auc <- forest.auc
  prjna384246.pred.result.list[[i]]$roc.df <- data.frame(fpr = forestperf@x.values[[1]], 
                                                         tpr = forestperf@y.values[[1]])
  prjna384246.pred.result.list[[i]]$roc.df$features <- i
}


prjna384246.roc.df <- data.frame(unlist(lapply(prjna384246.pred.result.list, 
                                               function(data){y <- data$forest.auc})),
                                 stringsAsFactors = F)
names(prjna384246.roc.df) <- "ROC"
prjna384246.roc.df$features <- rownames(prjna384246.roc.df)
prjna384246.roc.df$type <- c("all_spe", "dif_spe", "sig_spe",
                             "all_spe", "dif_spe", "sig_spe")


prjna384246.preform.df <- lapply(prjna384246.pred.result.list, 
                                 function(data){
                                   y <- data$roc.df
                                 }) %>% reduce(rbind)
  
prjna384246.preform.df$type <- substr(prjna384246.preform.df$features, 1, 7)
prjna384246.preform.df$HDC <- "including"
prjna384246.preform.df[which(prjna384246.preform.df$features %in% c("all.spe.nohdc","dif.spe.nohdc","sig.spe.nohdc")), "HDC"] <- "not including"

## -----external validation for predicting responses to therapies
prjna384246.preform.df %>%
  filter(HDC %in% "including") %>%
  mutate(type = factor(type,
                       levels = c("all.spe", "dif.spe", "sig.spe"))) %>%
  ggplot(aes(x = fpr, y = tpr)) + 
  geom_line(aes(color = type)) + 
  theme_classic() + 
  labs(x = "False Positve Rate", 
       y = "True Positive Rate", 
       title = "ROC plot") +
  theme(plot.title = element_text(hjust = .5, face = "bold", size = 11), 
        axis.text = element_text(size = 9)) +
  scale_color_manual(values = c("#e85a71","#56A902","#0072B2"))


## -----metabolic levels --------------------------##
## -----load pathways data of PRJNA384246 ---------##
prjna384246.path <- read.delim("prjna384246_pathabund_rela.txt",
                               header = T, sep = "\t", as.is = T)
names(prjna384246.path)[1] <- "pathway"
names(prjna384246.path)[2:172] <- substr(names(prjna384246.path)[2:172],
                                         1,10)
rownames(prjna384246.path) <- prjna384246.path$pathway

prjna384246.path.2 <- prjna384246.path[-grep(rownames(prjna384246.path), pattern = "|", fixed = T),]
prjna384246.path.2 <- prjna384246.path.2[-c(1:2),]

prjna384246.path.type <- data.frame(pathways = prjna384246.path.2$pathway, 
                                    stringsAsFactors = F)

prjna384246.path.type$Abb <- unlist(lapply(strsplit(prjna384246.path.type$pathways, 
                                            split = ":",
                                            fixed = T), 
                                   function(data){y <- data[[1]]}))

prjna384246.path.type$Abb.2 <- gsub(prjna384246.path.type$Abb, 
                                   pattern = "-", replacement = "_",
                                   fixed = T)

prjna384246.path.type$Abb.2 <- gsub(prjna384246.path.type$Abb.2, 
                                   pattern = "+", replacement = "_",
                                   fixed = T)
prjna384246.path.type$Abb.2 <- paste0("pathways_",
                                      prjna384246.path.type$Abb.2)

rownames(prjna384246.path.2) <- prjna384246.path.type$Abb.2
prjna384246.path.2 <- prjna384246.path.2[,-1]
prjna384246.path.t <- data.frame(t(prjna384246.path.2),
                                 stringsAsFactors = F)
prjna384246.path.t.match <- prjna384246.path.t[, names(prjna384246.path.t) %in% names(SRP057.abund.path.t)]

prjna384246.path.t.final <- prjna384246.path.t.match[rownames(prjna384246.path.t.match) %in% prjna384246.sample,]
prjna384246.path.t.final[,c("HDC", "Group")] <- prjna384246.meta.cd.week14[rownames(prjna384246.path.t.final), c("HDC","Group")]

prjna384246.abund.path.list <- list(all.path = prjna384246.path.t.final,
                                 dif.path = prjna384246.path.t.final[, names(prjna384246.path.t.final) %in% c(srp057.dif.path, "Group", "HDC")],
                                 sig.path = prjna384246.path.t.final[, names(prjna384246.path.t.final) %in% c(srp057.sig.path, "Group", "HDC")], 
                                 all.path.nohdc = prjna384246.path.t.final[, !names(prjna384246.path.t.final) %in% "HDC"],
                                 dif.path.nohdc = prjna384246.path.t.final[, names(prjna384246.path.t.final) %in% c(srp057.dif.path, "Group")],
                                 sig.path.nohdc = prjna384246.path.t.final[, names(prjna384246.path.t.final) %in% c(srp057.sig.path, "Group")])

prjna384246.pred.path.list <- list()
for (i in names(prjna384246.abund.path.list)) {
  prjna384246.pred.path.list[[i]] <- out.valid.pred(srp057.treat.path.rf.result[[i]]$rf.models, 
                                               test.data  = prjna384246.abund.path.list[[i]])
}


prjna384246.pred.path.df.list <- list()
prjna384246.pred.path.result.list <- list()

for (i in names(prjna384246.pred.path.list)) {
  prjna384246.pred.path.df.list[[i]] <- prjna384246.pred.path.list[[i]] %>%
    group_by(Sample_ID) %>%
    summarise(non = mean(non_responde),
              yes = mean(responde)) 
  prjna384246.pred.path.df.list[[i]] <- data.frame(prjna384246.pred.path.df.list[[i]], 
                                              stringsAsFactors = F)
  
  rownames(prjna384246.pred.path.df.list[[i]]) <- prjna384246.pred.path.df.list[[i]]$Sample_ID
  
  names(prjna384246.pred.path.df.list[[i]])[2:3] <- c("non_responde","responde")
  prjna384246.pred.path.df.list[[i]]$max <- apply(prjna384246.pred.path.df.list[[i]], 1, 
                                             function(data){y <- names(which.max(data))})
  prjna384246.pred.path.df.list[[i]]$max <- factor(prjna384246.pred.path.df.list[[i]]$max,
                                              levels = c("non_responde","responde"))
  
  prjna384246.pred.path.df.list[[i]]$true <- prjna384246.abund.path.list[[i]][rownames(prjna384246.pred.path.df.list[[i]]), 
                                                                      "Group"]
  prjna384246.pred.path.df.list[[i]]$true <- factor(prjna384246.pred.path.df.list[[i]]$true,
                                               levels = c("non_responde","responde"))
  
  forestpred <- ROCR::prediction( prjna384246.pred.path.df.list[[i]][,"responde"], 
                                  prjna384246.pred.path.df.list[[i]]$true, 
                                  label.ordering = c("non_responde","responde"));
  forestperf <- performance(forestpred, "tpr", "fpr");
  forest.auc <- performance(forestpred, "auc")@y.values[[1]]
  prjna384246.pred.path.result.list[[i]]$forestpred <- forestpred
  prjna384246.pred.path.result.list[[i]]$forestperf <- forestperf
  prjna384246.pred.path.result.list[[i]]$forest.auc <- forest.auc
  prjna384246.pred.path.result.list[[i]]$roc.df <- data.frame(fpr = forestperf@x.values[[1]], 
                                                         tpr = forestperf@y.values[[1]])
  prjna384246.pred.path.result.list[[i]]$roc.df$features <- i
}

unlist(lapply(prjna384246.pred.path.result.list, 
              function(data){y <- data$forest.auc}))
       
