SRP057.abund <- read.csv("CD/abundance.csv")
names(SRP057.abund)[1] <- "ID"

SRP057.abund.result <- data.phylum.func(SRP057.abund, i =10)
SRP057.abund.s <- SRP057.abund.result$data.abund.s


SRP057.meta <- read.delim("CD/SRP057027_meta.txt", 
                          header = T, sep = "\t", 
                          as.is = T)
rownames(SRP057.meta) <- SRP057.meta$run_id

SRP057.abund.s <- SRP057.abund.s[,names(SRP057.abund.s) %in% c("Species",SRP057.meta$run_id)]
rownames(SRP057.abund.s) <- SRP057.abund.s$Species
SRP057.abund.s <- SRP057.abund.s[,-1]/100

## ----0.01% average, 0.1%
SRP057.abund.s.filter <- noise.removal(SRP057.abund.s)
SRP057.abund.s.t <- data.frame(t(SRP057.abund.s.filter),
                               stringsAsFactors = F)
SRP057.abund.s.t[,c("HDC", "Group")] <- SRP057.meta[rownames(SRP057.abund.s.t), c("our_HDC", "Sample.group")]

## -----HDC correlations----------##
srp057.cor.list <- list()
srp057.cor.list[["Baseline_control"]] <- rcorr(as.matrix(subset(SRP057.abund.s.t, Group %in% c("Baseline","Control"))[,1:224]),
                                               type = "spearman")
for (i in c("Week1","Week4","Week8")) {
  srp057.cor.list[[i]] <- rcorr(as.matrix(subset(SRP057.abund.s.t, Group %in% i)[,1:224]),
                                type = "spearman")
}

srp057.cor.r <- lapply(srp057.cor.list, function(data){
  y <- data.frame(data$r,
                  stringsAsFactors = F);
  y.P <- data.frame(data$P,
                  stringsAsFactors = F);
  y$features <- rownames(y);
  y.P$features <- rownames(y.P);
  y.final <- y[,c("HDC", "features")];
  y.P.final <- y.P[,c("HDC", "features")];
  names(y.final)[1] <- "rho";
  names(y.P.final)[1] <- "P_value";
  y.merge <- merge(y.final, 
                   y.P.final,
                   by = "features",
                   all = T)
  return(y.merge)
})

for (pro in names(srp057.cor.r)) {
  srp057.cor.r[[pro]]$type <- pro
}

## -----select species whose P-value of correlations with HDC were below 0.001 in Baselines and controls------##
select.spe <- srp057.cor.r$Baseline_control %>%
  filter(P_value < 0.001) %>%
  pull(features)

## -----46 species whose P-value of correlations with HDC were below 0.001
srp057.sig.spe <- select.spe

srp057.cor.r.df <- srp057.cor.r %>% reduce(rbind) %>%
  filter(features %in% select.spe)

## -----define groups ----------##
srp057.cor.r.df$Sig_level <- "<0.01"
srp057.cor.r.df[which(srp057.cor.r.df$P_value < 0.05 & srp057.cor.r.df$P_value >=0.01), "Sig_level"] <- "0.01~0.05"
srp057.cor.r.df[which(srp057.cor.r.df$P_value >= 0.05), "Sig_level"] <- ">=0.05"

## -----change names of species ----##
srp057.cor.r.df$features <- gsub(srp057.cor.r.df$features, 
                                 pattern = "_", replacement = " ",
                                 fixed = T)
srp057.cor.r.df$features <- gsub(srp057.cor.r.df$features, 
                                 pattern = "unclassified", replacement = "spp.",
                                 fixed = T)

srp057.cor.r.df$features_2 <- gsub(srp057.cor.r.df$features, pattern = " ", replacement = "_", fixed = T)
srp057.cor.r.df$features_2 <- gsub(srp057.cor.r.df$features_2, pattern = "spp.", replacement = "unclassified", fixed = T)

srp057.cor.r.sort <- srp057.cor.r.df %>%
  filter(type %in% "Baseline_control") %>%
  arrange(rho)

## -----HDC correlated species plot ------------##
srp057.cor.plot <- srp057.cor.r.df %>%
  mutate(features = factor(features, 
                           levels = srp057.cor.r.sort$features)) %>%
  mutate(Sig_level = factor(Sig_level,
                            levels = c("<0.01", "0.01~0.05", ">=0.05"))) %>%
  ggplot(aes(x = features, y = rho)) +
  geom_point(aes(shape = type, color = Sig_level), 
             size = 3) +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values = c("#e85a71","#56A902","#999999")) +
  scale_shape_manual(values = c(1,2,4, 5)) +
  geom_hline(aes(yintercept = 0), colour="#990000", linetype="dashed", alpha = .8) + 
  theme(axis.text = element_text(size = 11),
        plot.title = element_text(vjust = 0.5,hjust = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour = "lightgrey", 
                                      size = 0.25))

## -------identified differential species ------------##
x <- subset(SRP057.abund.s.t, Group %in% c("Baseline","Control"))

x.result <- data.frame(features = names(x),
                       p.value = NA)
features <- names(x)
rownames(x.result) <- features
features <- features[!features %in% "Group"]
for (i in features) {
  a <- wilcox.test(x[,i] ~ x$Group, x)
  x.result[i, "p.value"] <- a$p.value
}
srp057.dif.df <- x.result
srp057.dif.df <- subset(srp057.dif.df, !features %in% c("HDC", "Group"))

srp057.dif.df$features <- as.character(srp057.dif.df$features)
srp057.dif.df <- srp057.dif.df %>% 
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, 
                        method = "fdr"))
srp057.dif.spe <- srp057.dif.df %>%
  filter(fdr < 0.05) %>%
  pull(features)

length(intersect(srp057.dif.spe, srp057.sig.spe)) ## 31 matching species


## -----
SRP057.abund.s.t.base <- subset(SRP057.abund.s.t, Group %in% c("Baseline","Control"))
SRP057.abund.s.t.base[which(SRP057.abund.s.t.base$Group %in% "Baseline"), "Group"] <- "CD"
srp057.abund.base.list <- list(all.spe = SRP057.abund.s.t.base,
                               dif.spe = SRP057.abund.s.t.base[, names(SRP057.abund.s.t.base) %in% c(srp057.dif.spe, "Group", "HDC")],
                               sig.spe = SRP057.abund.s.t.base[, names(SRP057.abund.s.t.base) %in% c(srp057.sig.spe, "Group", "HDC")],
                               all.spe.nohdc = SRP057.abund.s.t.base[, !names(SRP057.abund.s.t.base) %in% "HDC"],
                               dif.spe.nohdc = SRP057.abund.s.t.base[, names(SRP057.abund.s.t.base) %in% c(srp057.dif.spe, "Group")],
                               sig.spe.nohdc = SRP057.abund.s.t.base[, names(SRP057.abund.s.t.base) %in% c(srp057.sig.spe, "Group")])

## ----cross-validation for SRP057027 species profile --------##
save(srp057.abund.base.list, 
     file = "SRP057027_base_spe_list.Rdata")

srp057.base.rf.result <- lapply(srp057.abund.base.list, repeated.rf.confusion)

load("srp057_base_spe_rfmodels.RData")
srp057.df.group <- data.frame(Group = srp057.abund.base.list$all.spe$Group, 
                              stringsAsFactors = F)
rownames(srp057.df.group) <- rownames(srp057.abund.base.list$all.spe)

roc.plot <- function(data, data.df.group){
  data.sig.result.pred <- lapply(data$pred.matrix.list,function(data){ y <- reduce(data, rbind);name <- names(y);names(y) <- c("Group","Sample"); y.2 <- y %>% group_by(Sample) %>% summarise(mean = mean(Group)); names(y.2) <- rev(name); return(y.2)})
  data.sig.result.pred.df <- data.sig.result.pred %>% reduce(merge, by = "Sample", all = T)
  data.sig.result.pred.df$max <- apply(data.sig.result.pred.df, 1, function(data){y <- names(which.max(data))})
  rownames(data.sig.result.pred.df) <- data.sig.result.pred.df$Sample
  data.sig.result.pred.df <- data.sig.result.pred.df[rownames(data.df.group), ]
  data.sig.result.pred.df$max <- factor(data.sig.result.pred.df$max)
  data.sig.result.pred.df$true <- data.df.group$Group
  data.sig.result.pred.df$true <- factor(data.sig.result.pred.df$true)
  
  forestpred <- prediction(data.sig.result.pred.df[,3], data.sig.result.pred.df$true);
  forestperf <- performance(forestpred, "tpr", "fpr");
  forest.auc <- performance(forestpred, "auc")@y.values[[1]]
  data.roc.df <- data.frame(fpr = forestperf@x.values[[1]], 
                            tpr = forestperf@y.values[[1]])
  
  data.roc.plot <- ggplot(data.roc.df, aes(x = fpr, y = tpr)) + geom_line() + theme_bw() + 
    labs(x = "False Positve Rate", y = "True Positive Rate", title = "ROC plot")  + 
    annotate( "text", x = .75, y = .75, label = paste0("AUC:",round(forest.auc, 2)), size = 4) +
    theme(plot.title = element_text(hjust = .5, face = "bold", size = 11), 
          axis.text = element_text(size = 9))
  return(list(pred.df =data.sig.result.pred.df,
              data.roc.df = data.roc.df,
              forest.auc = forest.auc,
              data.roc.plot = data.roc.plot))
}
srp057.base.rf.roc <- lapply(srp057.base.rf.result,
                             function(data){
                               y <- roc.plot(data, srp057.df.group)
                             })

## ----Importance scores for models based on all species+ HDC
srp057.base.import <- srp057.base.rf.result$all.spe$importance.df
srp057.base.import$features <- rownames(srp057.base.import)
srp057.base.import$Median <- apply(srp057.base.import[,1:100],1, median)

srp057.base.import.melt <- melt(srp057.base.import, 
                                id.vars = c("Median", "features"))

srp057.base.import.summary <- srp057.base.import %>%
  arrange(desc(Median))


srp057.spe.type <- data.frame(Species = names(SRP057.abund.s.t),
                              stringsAsFactors = F)
srp057.spe.type$type <- "Others"
srp057.spe.type[which(srp057.spe.type$Species %in% srp057.sig.spe), "type"] <- "HDC_related"
srp057.spe.type[which(srp057.spe.type$Species %in% srp057.dif.spe), "type"] <- "Dif spe"
srp057.spe.type[which(srp057.spe.type$Species %in% intersect(srp057.sig.spe, srp057.dif.spe)), "type"] <- "Both"
srp057.spe.type[which(srp057.spe.type$Species %in% "HDC"), "type"] <- "HDC"

srp057.spe.type <- subset(srp057.spe.type, !Species %in% "Group")
srp057.spe.type$Species_2 <- gsub(srp057.spe.type$Species, pattern = "_", replacement = " ", fixed = T)
srp057.spe.type$Species_2 <- gsub(srp057.spe.type$Species_2,  
                                 pattern = "unclassified", replacement = "spp.",
                                 fixed = T)


srp057.base.import.melt <- merge(srp057.base.import.melt, 
                                 srp057.spe.type, 
                                 by.x = "features", 
                                 by.y = "Species",
                                 all.x = T)

srp057.base.import.melt %>%
  filter(features %in% srp057.base.import.summary[1:30,]$features) %>%
  mutate(features = factor(features, 
                           levels = rev(srp057.base.import.summary[1:30,]$features))) %>%
  ggplot(aes(x = features, y = value, 
             color = type)) + 
  geom_boxplot() +
  coord_flip() +
  scale_color_manual(values = c("#e85a71","#56A902","#6AAFE6","#999999"))

