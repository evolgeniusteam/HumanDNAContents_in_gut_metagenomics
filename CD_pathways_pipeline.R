## -----pathways levels ------
## -----for distinguishing cases from controls -------##
## -----AND distinguishing remiss patients from non-responde patients
## -----------------##
##- -----models using pathways ----
srp057.abund.path <- pathabund.filter.adenoma$SRP057027[!rownames(pathabund.filter.adenoma$SRP057027) %in% c("UNMAPPED","UNINTEGRATED"),]
srp057.abund.path.filter <- noise.removal.path(srp057.abund.path, percent = 1e-06);

srp057.pathways.type <- data.frame(pathways = rownames(srp057.abund.path.filter), 
                                   stringsAsFactors = F)
srp057.pathways.type$Abb <- unlist(lapply(strsplit(srp057.pathways.type$pathways, 
                                     split = ":",
                                     fixed = T), 
                                   function(data){y <- data[[1]]}))

srp057.pathways.type$Abb.2 <- gsub(srp057.pathways.type$Abb, 
                                   pattern = "-", replacement = "_",
                                   fixed = T)

srp057.pathways.type$Abb.2 <- gsub(srp057.pathways.type$Abb.2, 
                                   pattern = "+", replacement = "_",
                                   fixed = T)
srp057.pathways.type$Abb.2 <- paste0("pathways_",
                                     srp057.pathways.type$Abb.2)

rownames(srp057.abund.path.filter) <- srp057.pathways.type$Abb.2


SRP057.abund.path.t <- data.frame(t(srp057.abund.path.filter),
                               stringsAsFactors = F)
SRP057.abund.path.t[,c("HDC", "Group")] <- SRP057.meta[rownames(SRP057.abund.path.t), c("our_HDC", "Sample.group")]

## -----------HDC related pathways -------------------##
srp057.cor.path.list <- list()
srp057.cor.path.list[["Baseline_control"]] <- rcorr(as.matrix(subset(SRP057.abund.path.t, 
                                                                     Group %in% c("Baseline","Control"))[,1:142]),
                                               type = "spearman")
for (i in c("Week1","Week4","Week8")) {
  srp057.cor.path.list[[i]] <- rcorr(as.matrix(subset(SRP057.abund.path.t, Group %in% i)[,1:142]),
                                type = "spearman")
}

srp057.cor.path.r <- lapply(srp057.cor.path.list, 
                            function(data){
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

for (pro in names(srp057.cor.path.r)) {
  srp057.cor.path.r[[pro]]$type <- pro
}

## -----select pathways whose P-value of correlationse with HDC <0.001 using baseline and controls------##
select.path <- srp057.cor.path.r$Baseline_control %>%
  filter(P_value < 0.001) %>%
  pull(features)

## -----P-value <0.001
srp057.sig.path <- select.path

srp057.cor.path.r.df <- srp057.cor.path.r %>% 
  reduce(rbind) %>%
  filter(features %in% select.path)

## -----define groups ----------##
srp057.cor.path.r.df$Sig_level <- "<0.01"
srp057.cor.path.r.df[which(srp057.cor.path.r.df$P_value < 0.05 & srp057.cor.path.r.df$P_value >=0.01), "Sig_level"] <- "0.01~0.05"
srp057.cor.path.r.df[which(srp057.cor.path.r.df$P_value >= 0.05), "Sig_level"] <- ">=0.05"

## -----47 consistently HDC-related pathways in four stage (Baseline+controls, week1, week4, week8)
srp057.sig.path_4 <- names(table(subset(srp057.cor.path.r.df, P_value <0.05)$features)[table(subset(srp057.cor.path.r.df, P_value <0.05)$features)==4])


srp057.cor.path.r.sort <- srp057.cor.path.r.df %>%
  filter(type %in% "Baseline_control") %>%
  arrange(rho)

## -----HDC correlated pathways plot ------------##
srp057.cor.path.plot <- srp057.cor.path.r.df %>%
  mutate(features = factor(features, 
                           levels = srp057.cor.path.r.sort$features)) %>%
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

## -------identified differential pathways ------------##
x <- subset(SRP057.abund.path.t, 
            Group %in% c("Baseline","Control"))

x.result <- data.frame(features = names(x),
                       p.value = NA)
features <- names(x)
rownames(x.result) <- features
features <- features[!features %in% "Group"]
for (i in features) {
  a <- wilcox.test(x[,i] ~ x$Group, x)
  x.result[i, "p.value"] <- a$p.value
}
srp057.dif.path.df <- x.result
srp057.dif.path.df <- subset(srp057.dif.path.df, !features %in% c("HDC", "Group"))

srp057.dif.path.df$features <- as.character(srp057.dif.path.df$features)
srp057.dif.path.df <- srp057.dif.path.df %>% 
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, 
                        method = "fdr"))
## ---7 differential pathways 
srp057.dif.path <- srp057.dif.path.df %>%
  filter(fdr < 0.05) %>%
  pull(features)

length(intersect(srp057.dif.path, srp057.sig.path)) ## 4 pathways 
length(intersect(srp057.dif.path, srp057.sig.path_4)) ## 2 pathways

## -----
SRP057.abund.path.base <- subset(SRP057.abund.path.t, Group %in% c("Baseline","Control"))
SRP057.abund.path.base[which(SRP057.abund.path.base$Group %in% "Baseline"), "Group"] <- "CD"
srp057.abund.path.base.list <- list(all.path = SRP057.abund.path.base,
                                    dif.path = SRP057.abund.path.base[, names(SRP057.abund.path.base) %in% c(srp057.dif.path, "Group", "HDC")],
                                    sig.path = SRP057.abund.path.base[, names(SRP057.abund.path.base) %in% c(srp057.sig.path, "Group", "HDC")],
                                    all.path.nohdc = SRP057.abund.path.base[, !names(SRP057.abund.path.base) %in% "HDC"],
                                    dif.path.nohdc = SRP057.abund.path.base[, names(SRP057.abund.path.base) %in% c(srp057.dif.path, "Group")],
                                    sig.path.nohdc = SRP057.abund.path.base[, names(SRP057.abund.path.base) %in% c(srp057.sig.path, "Group")])

## ----cross-validation for SRP057027 pathways profile --------##
save(srp057.abund.path.base.list, 
     file = "SRP057027_base_path_list.Rdata")

srp057.base.path.rf.result <- lapply(srp057.abund.path.base.list, repeated.rf.confusion)

load("srp057_base_path_rfmodels.RData")

srp057.base.path.rf.roc <- lapply(srp057.base.path.rf.result,
                             function(data){
                               y <- roc.plot(data, srp057.df.group)
                             })
unlist(lapply(srp057.base.path.rf.roc, function(data){y <- data$forest.auc}))

## ----Importance scores for models based on all species+ HDC
srp057.base.path.import <- srp057.base.path.rf.result$all.path$importance.df
srp057.base.path.import$features <- rownames(srp057.base.path.import)
srp057.base.path.import$Median <- apply(srp057.base.path.import[,1:100],1, median)

srp057.base.path.import.melt <- melt(srp057.base.path.import, 
                                id.vars = c("Median", "features"))

srp057.base.path.import.summary <- srp057.base.path.import %>%
  arrange(desc(Median))

## pathways type

srp057.pathways.type$type <- "Others"
srp057.pathways.type[which(srp057.pathways.type$Abb.2 %in% srp057.sig.path), "type"] <- "HDC_related"
srp057.pathways.type[which(srp057.pathways.type$Abb.2 %in% srp057.dif.path), "type"] <- "Dif spe"
srp057.pathways.type[which(srp057.pathways.type$Abb.2 %in% intersect(srp057.sig.path_4, srp057.dif.path)), "type"] <- "Both"
srp057.pathways.type[142, ] <- c("Human DNA contents", "HDC", "HDC", "HDC")
names(srp057.pathways.type)[1] <- "pathways"

srp057.base.path.import.melt <- merge(srp057.base.path.import.melt, 
                                 srp057.pathways.type, 
                                 by.x = "features", 
                                 by.y = "Abb.2",
                                 all.x = T)

srp057.base.path.import.melt %>%
  filter(features %in% srp057.base.path.import.summary[1:30,]$features) %>%
  mutate(features = factor(features, 
                           levels = rev(srp057.base.path.import.summary[1:30,]$features))) %>%
  ggplot(aes(x = features, y = value, 
             color = type)) + 
  geom_boxplot() +
  coord_flip() +
  scale_color_manual(values = c("#e85a71","#56A902","#6AAFE6", "#0072B2","#999999"))

## cross-validation classifiers for predicting treatment response
## treatment data 
srp057.abund.path.treat <- SRP057.abund.path.t[rownames(SRP057.abund.path.t) %in% rownames(srp057.treat.df.group.2), ]
srp057.abund.path.treat$Group <- srp057.treat.df.group.2[rownames(srp057.abund.path.treat),"Group"]

srp057.abund.path.treat.list <- list(all.path = srp057.abund.path.treat,
                                     dif.path = srp057.abund.path.treat[, names(srp057.abund.path.treat) %in% c(srp057.dif.path, "Group", "HDC")],
                                     sig.path = srp057.abund.path.treat[, names(srp057.abund.path.treat) %in% c(srp057.sig.path, "Group", "HDC")],
                                     all.path.nohdc = srp057.abund.path.treat[, !names(srp057.abund.path.treat) %in% "HDC"],
                                     dif.path.nohdc = srp057.abund.path.treat[, names(srp057.abund.path.treat) %in% c(srp057.dif.path, "Group")],
                                     sig.path.nohdc = srp057.abund.path.treat[, names(srp057.abund.path.treat) %in% c(srp057.sig.path, "Group")])

## ----cross-validation for patients species profile --------##
save(srp057.abund.path.treat.list, 
     file = "SRP057027_treat_path_list.Rdata")

srp057.treat.path.rf.result <- lapply(srp057.abund.path.treat.list, repeated.rf.confusion)

load("srp057_treat_path_rfmodels.RData")
srp057.treat.path.rf.roc <- lapply(srp057.treat.path.rf.result,
                                  function(data){
                                    y <- roc.plot(data, srp057.treat.df.group.2)
                                  })
unlist(lapply(srp057.treat.path.rf.roc, function(data){y <- data$forest.auc}))

