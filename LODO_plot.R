## ----- example -----##
cl <- makeCluster(7)
registerDoParallel(cl)
study <- unique(crc.meta$Project)
all.spe.cross <- foreach(project = study) %dopar% cross.models(project = project, metadata=crc.meta, rl.data = spe.data.df)
all.spe.lodo <- foreach(project = study) %dopar% lodo.models(project = project, metadata=crc.meta, rl.data = spe.data.df)

spe.data.df <- spe.data.df[!rownames(spe.data.df) %in% "HDC", ]
all.spe.cross.nohdc <- foreach(project = study) %dopar% cross.models(project = project, metadata=crc.meta, rl.data = spe.data.df)
sig.spe.lodo.nohdc <- foreach(project = study) %dopar% lodo.models(project = project, metadata=crc.meta, rl.data = spe.data.df)

stopCluster(cl)
## ----- get lodo results of six datasets ----##


## -----LODO list ---- ##
spe.lodo.list <- c(sig.spe.lodo, sig.spe.lodo.nohdc, wt.spe.lodo, wt.spe.lodo.nohdc,all.spe.lodo, all.spe.lodo.nohdc)
## -----
study <- unique(crc.meta$Project)
names(spe.lodo.list) <- c(paste0(study, "_sigspe_HDC_rf"),
                          paste0(study, "_sigspe_rf"),
                          paste0(study, "_wtspe_HDC_rf"),
                          paste0(study, "_wtspe_rf"),
                          paste0(study, "_allspe_HDC_rf"),
                          paste0(study, "_allspe_rf"))
lodo.type <- c("sigspe_HDC","sigspe","wtspe_HDC","wtspe","allspe_HDC","allspe")

pred.matrix <- matrix(NA, nrow = nrow(crc.meta),ncol = 6,
                      dimnames = list(crc.meta$Sample_ID, 
                                      lodo.type))

for (project in study) {
  for (type in lodo.type) {
    siamcat <- spe.lodo.list[[paste0(project, "_", type, "_rf" )]]$models
    siamcat.imp <- rownames(spe.lodo.list[[paste0(project, "_", type, "_rf" )]]$features.imp)
    meta.test <- crc.meta %>% filter(Project == project)
    meta.test <- data.frame(meta.test)
    rownames(meta.test) <- meta.test$Sample_ID
    feat.test <- rl.data[siamcat.imp, rownames(meta.test)]
    
    siamcat.test <- siamcat(feat = feat.test)
    siamcat.test.pred <- make.predictions(siamcat = siamcat, 
                                          siamcat.holdout = siamcat.test, normalize.holdout = F)
    temp <-  rowMeans(pred_matrix(siamcat.test.pred))
    pred.matrix[names(temp), type] <- temp
  }
}

pred.matrix.df <- data.frame(pred.matrix, stringsAsFactors = F)
pred.matrix.df$Sample_ID <- rownames(pred.matrix.df)
pred.matrix.df.merge <- merge(pred.matrix.df, crc.meta, by = "Sample_ID", all = T)

spe.lodo <- data.frame()
for (project in study){
  for (type in lodo.type) {
    cases <- pred.matrix.df.merge %>%
      filter(Project == project) %>% 
      filter(Group=='CRC') %>% 
      pull(type) 
    controls <- pred.matrix.df.merge %>%
      filter(Project == project) %>% 
      filter(Group=='CTR') %>% 
      pull(type)                  
    temp <- roc(cases=cases, controls=controls)
    spe.lodo <- bind_rows(
      spe.lodo, tibble(study.test=project, AUC=c(temp$auc), type.2 = type ))
  }
}
names(spe.lodo) <- c("Project","AUC","SpeciesType")

col.lodo.heatmap2 <- c("#C0C0C0","#ffbcbc","#ff7f7f","#FF0000")
spe.lodo.plot1 <- spe.lodo %>% mutate(SpeciesType=factor(SpeciesType)) %>%
  ggplot(aes(y=SpeciesType, x=Project, fill=AUC)) +
  geom_tile() + theme_bw() +
  # test in tiles
  geom_text(aes_string(label="format(AUC, digits=3)"), col='white', size=7)+
  # color scheme
  scale_fill_gradientn(colours = col.lodo.heatmap2, limits=c(0.6, 1)) +
  # axis position/remove boxes/ticks/facet background/etc.
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x.top = element_text(angle=45, hjust=.1), 
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text = element_text(size = 11)) + 
  xlab('Project') + ylab('Species Type') + 
  scale_color_manual(values=c('#FFFFFF00', 'grey'), guide=FALSE) + 
  scale_size_manual(values=c(0, 3), guide=FALSE)
spe.lodo.plot2 <- spe.lodo %>% 
  group_by(SpeciesType) %>% 
  summarise(AUROC=mean(AUC)) %>%
  ggplot(aes(y=SpeciesType, x=1, fill=AUROC)) + 
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUROC, digits=3)"), col='white', size=7)+
  scale_fill_gradientn(colours = col.lodo.heatmap2, limits=c(0.6, 1), 
                       guide=FALSE) + 
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text = element_text(size = 11)) + 
  xlab('Model Average') + ylab('')

spe.lodo.plot3 <- spe.lodo %>% 
  group_by(Project) %>% 
  summarise(AUROC=mean(AUC)) %>%
  ggplot(aes(x=Project, y=1, fill=AUROC)) + 
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUROC, digits=3)"), col='white', size=7)+
  scale_fill_gradientn(colours = col.lodo.heatmap2, limits=c(0.6, 1), 
                       guide=FALSE) + 
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text = element_text(size = 11)) + 
  ylab('Model Average') + xlab('')

