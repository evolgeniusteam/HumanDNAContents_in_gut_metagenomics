options(stringsAsFactors = F)

crc.meta <- read.delim("allCRC_metadata.txt", header = T, sep = "\t", as.is = T)
study_spe_vector <- read.delim("study_spe_vector" , header = F, as.is = T)

spe.lodo.features <- lapply(all.spe.lodo, function(data){y <- data$features.imp; y.2 <- rowMeans(y)})
spe.lodo.features.df <- Reduce(cbind, spe.lodo.features)
spe.lodo.features.df <- data.frame(spe.lodo.features.df)
names(spe.lodo.features.df) <- paste0(study, "_allspe_HDC")

spe.lodo.features.df$Means <- rowMeans(spe.lodo.features.df[,1:7])
spe.lodo.features.df$Features <- rownames(spe.lodo.features.df)
spe.lodo.features.df <- spe.lodo.features.df %>% arrange(Means)

spe.lodo.features.df$type <- "Others"
spe.lodo.features.df[which(spe.lodo.features.df$Features %in% "HDC"), "type"] <- "HDC"
spe.lodo.features.df[which(spe.lodo.features.df$Features %in% spe.wt.list.spe.vector), "type"] <- "Dif"
spe.lodo.features.df[which(spe.lodo.features.df$Features %in% sig.spe.vector), "type"] <- "HDC correlated"
spe.lodo.features.df[which(spe.lodo.features.df$Features %in% intersect(spe.wt.list.spe.vector, sig.spe.vector)), "type"] <- "Both"

spe.lodo.features.df.melt <- melt(spe.lodo.features.df, id.vars = c("Features","Means","Type"))
spe.lodo.imp.plot <- spe.lodo.features.df.melt %>% 
  mutate(Features = factor(Features, levels = spe.lodo.features.df$Features)) %>%
  filter(Features %in% spe.lodo.features.df[1:30,]$Features) %>%
  ggplot(aes(x = Features, y = value, fill = type)) + geom_boxplot() +
  coord_flip() +
  labs(y = "Mean decrease in accuracy", title = "Importance Score") +
  scale_fill_discrete(breaks = unique(spe.lodo.features.df.melt$Type)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(),panel.border=element_blank(), 
        axis.line=element_line(size=.5,colour="black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", vjust = .5, hjust = .5)) + 
  scale_fill_manual(values = c("#e85a71","#56A902","#6AAFE6","#999999"))
