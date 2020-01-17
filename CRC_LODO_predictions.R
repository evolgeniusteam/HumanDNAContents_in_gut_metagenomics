## ------prepare feature types -----##
all.spe.df.type <- subset(all.spe.df.merge, 
                          Species %in% names(spe.data))
all.spe.df.type$type <- "Others"
all.spe.df.type[which(all.spe.df.type$Species %in% crc.sig.spe), "type"] <- "HDC_related"
all.spe.df.type[which(all.spe.df.type$Species %in% crc.kw.spe), "type"] <- "Dif"
all.spe.df.type[which(all.spe.df.type$Species %in% "HDC"), "type"] <- "HDC"
all.spe.df.type[which(all.spe.df.type$Species %in% intersect(crc.sig.spe, crc.kw.spe)), "type"] <- "Both"
table(all.spe.df.type$type)

rownames(all.spe.df.type) <- all.spe.df.type$Species

## ------LODO models -------------##
## ------Based on species levels --------##

load("spe_LODO_models.RData")
names(spe.lodo.models) <- names(spe.list)
names(spe.lodo.models$all.spe.data) <- study
names(spe.lodo.models$kw.spe.data) <- study
names(spe.lodo.models$sig.spe.data) <- study


names(spe.lodo.nohdc.models) <- names(spe.nohdc.list)
spe.lodo.nohdc.models <- lapply(spe.lodo.nohdc.models, 
                                function(data){
                                  names(data) <- study;
                                  return(data)
                                })

## ------------- function for predictions -------------##
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

## ------AUROC for models with HDC + species -------
spe.lodo.pred <- list()

for (data in names(spe.list)) {
  spe.lodo.pred[[data]] <- list()
  data.crc <- spe.list[[data]]
  for (pro in study) {
    loso.sample <- subset(crc.meta, Project %in% pro)$Sample_ID
    test.data <- data.crc[rownames(data.crc) %in% loso.sample, ]
    spe.lodo.pred[[data]][[pro]] <- list()
    spe.lodo.pred[[data]][[pro]] <- lodo.pred.func(rf.models = spe.lodo.models[[data]][[pro]]$rf.models,
                                                   metadata = crc.meta,
                                                   test.data = test.data)
    
  }
}

spe.pred.matrix <- matrix(data = NA, 
                          nrow = 6, 
                          ncol = 7,
                          dimnames = list(c("allspe + HDC", "difspe + HDC", "sigspe + HDC",
                                            "allspe", "difspe", "sigspe"),
                                          study))

for (i in 1:3) {
  a <- lapply(spe.lodo.pred[[i]], function(data){
    y <- data$forest.auc
  }) %>% unlist
  spe.pred.matrix[i, ] <- a[colnames(spe.pred.matrix)]
}

##  ------AUROC for models with species -------

spe.lodo.nohdc.pred <- list()

for (data in names(spe.nohdc.list)) {
  spe.lodo.nohdc.pred[[data]] <- list()
  data.crc <- spe.nohdc.list[[data]]
  for (pro in study) {
    loso.sample <- subset(crc.meta, Project %in% pro)$Sample_ID
    test.data <- data.crc[rownames(data.crc) %in% loso.sample, ]
    spe.lodo.nohdc.pred[[data]][[pro]] <- list()
    spe.lodo.nohdc.pred[[data]][[pro]] <- lodo.pred.func(rf.models = spe.lodo.nohdc.models[[data]][[pro]]$rf.models,
                                                   metadata = crc.meta,
                                                   test.data = test.data)
    
  }
}

for (i in 1:3) {
  a <- lapply(spe.lodo.nohdc.pred[[i]], function(data){
    y <- data$forest.auc
  }) %>% unlist
  spe.pred.matrix[i+3, ] <- a[colnames(spe.pred.matrix)]
}

## -----LODO AUROC melt dataframe ---------------##
# model average
spe.lodo.df.melt <- data.frame(spe.pred.matrix, 
                           stringsAsFactors = F)
spe.lodo.df.melt$type <- rownames(spe.lodo.df.melt)
spe.lodo.df.melt <- melt(spe.lodo.df.melt, id.vars = "type")
names(spe.lodo.df.melt) <- c("SpeciesType","Project","AUC")

col.lodo.heatmap <- c("#C0C0C0","#ffbcbc","#ff7f7f","#FF0000")

spe.lodo.plot1 <- spe.lodo.df.melt %>% 
  mutate(SpeciesType=factor(SpeciesType,
                            levels = c("allspe + HDC", "allspe",
                                       "sigspe + HDC", "sigspe",
                                       "difspe + HDC", "difspe"))) %>%
  ggplot(aes(y=SpeciesType, x=Project, fill=AUC)) +
  geom_tile() + theme_bw() +
  # test in tiles
  geom_text(aes_string(label="format(AUC, digits=2)"), col='white', size=7)+
  # color scheme
  scale_fill_gradientn(colours = col.lodo.heatmap, limits=c(0.5, 1)) +
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Project') + ylab('Species Type') + 
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 3), guide=FALSE) +
  theme(axis.text = element_text(size = 12))

spe.lodo.plot2 <- spe.lodo.df.melt %>% 
  mutate(SpeciesType=factor(SpeciesType,
                            levels = c("allspe + HDC", "allspe",
                                       "sigspe + HDC", "sigspe",
                                       "difspe + HDC", "difspe"))) %>% 
  group_by(SpeciesType) %>% 
  summarise(AUROC=mean(AUC)) %>%
  ggplot(aes(y=SpeciesType, x=1, fill=AUROC)) + 
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUROC, digits=2)"), col='white', size=7)+
  scale_fill_gradientn(colours = col.lodo.heatmap, limits=c(0.5, 1), 
                       guide=FALSE) + 
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank()) + 
  xlab('Model Average') + ylab('')

spe.lodo.plot.all <- plot_grid(spe.lodo.plot1, spe.lodo.plot2, 
          rel_widths = c(5/6, 1/6), align = 'h')

## -----Importance dataframe --------------------##
spe.import.list <- lapply(spe.lodo.models$all.spe.data, 
                          function(data){
                            y <- data$importance.df;
                            y$features <- rownames(y);
                            y$importance <- rowMeans(y[,1:100], 
                                                     na.rm = T);
                            y.2 <- y[,c("features", "importance")];
                            return(y.2)
                          })
for (pro in names(spe.import.list)) {
  spe.import.list[[pro]]$project <- pro
}
spe.import.df <- spe.import.list %>% reduce(rbind)

spe.import.df.median <- spe.import.df %>%
  group_by(features) %>%
  summarise(median = median(importance),
            mean = mean(importance)) %>%
  arrange(desc(median))
spe.import.df.median <- data.frame(spe.import.df.median,
                                   stringsAsFactors = F)


## ----
spe.import.df.final <- merge(spe.import.df, 
                             all.spe.df.type, 
                             by.y = "Species", 
                             by.x = "features",
                             all.x = T)
                             
## ---- species importance plot -----------##

spe.import.df.final %>%
  filter(features %in% spe.import.df.median[1:30,]$features) %>%
  mutate(features = factor(features, 
                           levels = rev(spe.import.df.median[1:30,]$features))) %>%
  ggplot(aes(x = features, y = importance,
             color = type)) +
  geom_boxplot() +
  coord_flip() +
  scale_color_manual(values = c("#e85a71","#56A902","#6AAFE6", "#0072B2","#999999"))

save(spe.list, 
     spe.nohdc.list, 
     spe.lodo.models, 
     spe.lodo.nohdc.models, 
     spe.lodo.pred, 
     spe.lodo.nohdc.pred, 
     spe.lodo.df.melt, 
     spe.pred.matrix, 
     spe.import.df.final, 
     col.lodo.heatmap, 
     all.spe.df.type,
     spe.lodo.plot1, spe.lodo.plot2,
     file = "spe_LODO_models.RData")

rm(spe.lodo.models, spe.lodo.nohdc.models)
