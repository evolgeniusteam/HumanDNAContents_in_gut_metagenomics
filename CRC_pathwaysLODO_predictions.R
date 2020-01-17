
load("path_LODO_models.RData")
names(path.lodo.models) <- names(path.list)
names(path.lodo.nohdc.models) <- names(path.nohdc.list)

path.lodo.models <- lapply(path.lodo.models, 
                                function(data){
                                  names(data) <- study;
                                  return(data)
                                })

path.lodo.nohdc.models <- lapply(path.lodo.nohdc.models, 
                                function(data){
                                  names(data) <- study;
                                  return(data)
                                })

path.lodo.pred <- list()

for (data in names(path.list)) {
  path.lodo.pred[[data]] <- list()
  data.crc <- path.list[[data]]
  for (pro in study) {
    loso.sample <- subset(crc.meta, Project %in% pro)$Sample_ID
    test.data <- data.crc[rownames(data.crc) %in% loso.sample, ]
    path.lodo.pred[[data]][[pro]] <- list()
    path.lodo.pred[[data]][[pro]] <- lodo.pred.func(rf.models = path.lodo.models[[data]][[pro]]$rf.models,
                                                   metadata = crc.meta,
                                                   test.data = test.data)
  }
}

path.lodo.nohdc.pred <- list()

for (data in names(path.nohdc.list)) {
  path.lodo.nohdc.pred[[data]] <- list()
  data.crc <- path.nohdc.list[[data]]
  for (pro in study) {
    loso.sample <- subset(crc.meta, Project %in% pro)$Sample_ID
    test.data <- data.crc[rownames(data.crc) %in% loso.sample, ]
    path.lodo.nohdc.pred[[data]][[pro]] <- list()
    path.lodo.nohdc.pred[[data]][[pro]] <- lodo.pred.func(rf.models = path.lodo.nohdc.models[[data]][[pro]]$rf.models,
                                                         metadata = crc.meta,
                                                         test.data = test.data)
    
  }
}


path.pred.matrix <- matrix(data = NA, 
                          nrow = 6, 
                          ncol = 7,
                          dimnames = list(c("allpath + HDC", "difpath + HDC", "sigpath + HDC",
                                            "allpath", "difpath", "sigpath"),
                                          study))

for (i in 1:3) {
  a <- lapply(path.lodo.pred[[i]], function(data){
    y <- data$forest.auc
  }) %>% unlist
  path.pred.matrix[i, ] <- a[colnames(path.pred.matrix)]
}

for (i in 1:3) {
  a <- lapply(path.lodo.nohdc.pred[[i]], function(data){
    y <- data$forest.auc
  }) %>% unlist
  path.pred.matrix[i+3, ] <- a[colnames(spe.pred.matrix)]
}

# model average
path.lodo.df.melt <- data.frame(path.pred.matrix, 
                               stringsAsFactors = F)
path.lodo.df.melt$type <- rownames(path.lodo.df.melt)
path.lodo.df.melt <- melt(path.lodo.df.melt, id.vars = "type")
names(path.lodo.df.melt) <- c("SpeciesType","Project","AUC")

col.lodo.heatmap <- c("#C0C0C0","#ffbcbc","#ff7f7f","#FF0000")

path.lodo.plot1 <- path.lodo.df.melt %>% 
  mutate(SpeciesType=factor(SpeciesType,
                            levels = c("allpath + HDC", "allpath",
                                       "sigpath + HDC", "sigpath",
                                       "difpath + HDC", "difpath"))) %>%
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

path.lodo.plot2 <- path.lodo.df.melt %>% 
  mutate(SpeciesType=factor(SpeciesType,
                            levels = c("allpath + HDC", "allpath",
                                       "sigpath + HDC", "sigpath",
                                       "difpath + HDC", "difpath"))) %>% 
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

path.lodo.plot.all <- plot_grid(path.lodo.plot1, path.lodo.plot2, 
                               rel_widths = c(5/6, 1/6), align = 'h')


## -----Importance dataframe --------------------##
path.import.list <- lapply(path.lodo.models$all.path.data, 
                          function(data){
                            y <- data$importance.df;
                            y$features <- rownames(y);
                            y$importance <- rowMeans(y[,1:100], 
                                                     na.rm = T);
                            y.2 <- y[,c("features", "importance")];
                            return(y.2)
                          })
for (pro in names(path.import.list)) {
  path.import.list[[pro]]$project <- pro
}
path.import.df <- path.import.list %>% reduce(rbind)

## ----
all.path.df.type <- lapply(pathabund.filter.adenoma, rownames) %>%
  reduce(c) 

all.path.df.type <- unique(all.path.df.type)
all.path.df.type <- data.frame(all.path.df.type, 
                               stringsAsFactors = F)
names(all.path.df.type) <- "pathways"
all.path.df.type$Abb <- unlist(lapply(strsplit(all.path.df.type$pathways,
                                        split = ":", fixed = T), 
                               function(data){
                                y <- data[[1]]
                               }))
all.path.df.type$wholename <- unlist(lapply(strsplit(all.path.df.type$pathways,
                                               split = ":", fixed = T), 
                                      function(data){
                                        y <- data[2]
                                      }))
all.path.df.type <- all.path.df.type[-c(1:2), ]

all.path.df.type$Abb.2 <- gsub(all.path.df.type$Abb, 
                               pattern = "-", replacement = "_",
                               fixed = T)
all.path.df.type$Abb.2 <- gsub(all.path.df.type$Abb.2, 
                               pattern = "+", replacement = "_",
                               fixed = T)

all.path.df.type$Abb.2 <- paste0("pathways_", 
                                 all.path.df.type$Abb.2)

all.path.df.type <- subset(all.path.df.type, 
                           Abb.2 %in% names(path.data))
all.path.df.type[346, ] <- c("HDC","HDC","host DNA contents", 
                             "HDC")

all.path.df.type$type <- "Others"
all.path.df.type[which(all.path.df.type$Abb.2 %in% sig.path), "type"] <- "HDC related"
all.path.df.type[which(all.path.df.type$Abb.2 %in% kw.path), "type"] <- "Dif"
all.path.df.type[which(all.path.df.type$Abb.2 %in% intersect(sig.path, kw.path)), "type"] <- "Both"
all.path.df.type[which(all.path.df.type$Abb.2 %in% "HDC"), "type"] <- "HDC"
table(all.path.df.type$type)


path.import.df.final <- merge(path.import.df, 
                              all.path.df.type, 
                             by.x = "features", 
                             by.y = "Abb.2",
                             all.x = T)
## ---- species importance plot -----------##

path.import.df.median <- path.import.df.final %>%
  group_by(wholename) %>%
  summarise(median = median(importance),
            mean = mean(importance)) %>%
  arrange(desc(median))

path.import.df.median <- data.frame(path.import.df.median,
                                    stringsAsFactors = F)

path.import.df.final %>%
  filter(wholename %in% path.import.df.median[1:30,]$wholename) %>%
  mutate(wholename = factor(wholename, 
                           levels = rev(path.import.df.median[1:30,]$wholename))) %>%
  ggplot(aes(x = wholename, y = importance,
             color = type)) +
  geom_boxplot() +
  coord_flip() +
  scale_color_manual(values = c("#e85a71","#56A902","#6AAFE6", "#0072B2","#999999"))


save(path.list, 
     path.nohdc.list, 
     path.lodo.models, 
     path.lodo.nohdc.models, 
     path.lodo.pred, 
     path.lodo.nohdc.pred, 
     path.lodo.df.melt, 
     path.pred.matrix, 
     path.import.df.final, 
     col.lodo.heatmap, 
     all.path.df.type,
     path.lodo.plot1, path.lodo.plot2,
     file = "path_LODO_models.RData")
