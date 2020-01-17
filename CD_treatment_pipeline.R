## ---------------load metadata------------------------##
SRP057.meta.treat <- subset(SRP057.meta,Sample.Type %in% "Crohn")
##----------------select samples whose FCP in baseline were over 250ug/g
SRP057.meta.treat.week <- subset(SRP057.meta.treat, 
                                   Sample.group %in% "Baseline" & FCP >250)
SRP057.meta.treat.sample <- SRP057.meta.treat.week$sample
rownames(SRP057.meta.treat.week) <- SRP057.meta.treat.week$sample
SRP057.meta.treat.week <- SRP057.meta.treat.week[rownames(SRP057.meta.treat.week) %in% names(table(SRP057.meta$sample)[table(SRP057.meta$sample)==4]), ]

SRP057.meta.treat.week$week_1 <- "0"
sample <- subset(SRP057.meta.treat, Sample.group %in% "Week1" & sample %in% SRP057.meta.treat.sample & FCP <= 250)$sample
sample <- as.character(sample)
SRP057.meta.treat.week[rownames(SRP057.meta.treat.week) %in% sample, "week_1"] <- 1

SRP057.meta.treat.week$week_4 <- "0"
sample <- subset(SRP057.meta.treat, Sample.group %in% "Week4" & sample %in% SRP057.meta.treat.sample & FCP <= 250)$sample
sample <- as.character(sample)
SRP057.meta.treat.week[rownames(SRP057.meta.treat.week) %in% sample,
                       "week_4"] <- 1

SRP057.meta.treat.week$week_8 <- "0"
sample <- subset(SRP057.meta.treat, Sample.group %in% "Week8" & sample %in% SRP057.meta.treat.sample & FCP <= 250)$sample
sample <- as.character(sample)

SRP057.meta.treat.week[rownames(SRP057.meta.treat.week) %in% sample,
                       "week_8"] <- 1

SRP057.meta.treat.week.melt <- SRP057.meta.treat.week[,c("sample","week_1","week_4","week_8")]
names(SRP057.meta.treat.week.melt)[2:4] <- c("Week1","Week4","Week8" )
SRP057.meta.treat.week.melt <- melt(SRP057.meta.treat.week.melt, 
                                    id.vars = "sample")

srp057.treat.df.group <- merge(subset(SRP057.meta, sample %in% rownames(SRP057.meta.treat.week) & !Sample.group %in% "Baseline")[,c("sample","run_id", "Sample.group")],
                               SRP057.meta.treat.week.melt, by.x= c("sample", "Sample.group"),
                               by.y = c("sample","variable"),
                               all.x = T)
srp057.treat.df.group$Group <- "non_responde"
srp057.treat.df.group[which(srp057.treat.df.group$value %in% "1"), "Group"] <- "responde"

srp057.treat.df.group.2 <- data.frame(Group = srp057.treat.df.group$Group,
                                      stringsAsFactors = F)
rownames(srp057.treat.df.group.2) <- srp057.treat.df.group$run_id


srp057.abund.treat.s <- SRP057.abund.s.t[rownames(SRP057.abund.s.t) %in% rownames(srp057.treat.df.group.2), ]
srp057.abund.treat.s$Group <- srp057.treat.df.group.2[rownames(srp057.abund.treat.s),"Group"]

srp057.abund.treat.list <- list(all.spe = srp057.abund.treat.s,
                               dif.spe = srp057.abund.treat.s[, names(srp057.abund.treat.s) %in% c(srp057.dif.spe, "Group", "HDC")],
                               sig.spe = srp057.abund.treat.s[, names(srp057.abund.treat.s) %in% c(srp057.sig.spe, "Group", "HDC")],
                               all.spe.nohdc = srp057.abund.treat.s[, !names(srp057.abund.treat.s) %in% "HDC"],
                               dif.spe.nohdc = srp057.abund.treat.s[, names(srp057.abund.treat.s) %in% c(srp057.dif.spe, "Group")],
                               sig.spe.nohdc = srp057.abund.treat.s[, names(srp057.abund.treat.s) %in% c(srp057.sig.spe, "Group")])

## ----cross-validation for patients using species profile --------##
save(srp057.abund.treat.list, 
     file = "SRP057027_treat_spe_list.Rdata")
srp057.treat.rf.result <- lapply(srp057.abund.treat.list, repeated.rf.confusion)

load("srp057_treat_spe_rfmodels.RData")
srp057.treat.rf.roc <- lapply(srp057.treat.rf.result,
                             function(data){
                               y <- roc.plot(data, srp057.treat.df.group.2)
                             })

unlist(lapply(srp057.treat.rf.roc, function(data){y <- data$forest.auc}))

## ----Importance scores for treating models
## ----based on HDC related species+ HDC
srp057.treat.import <- srp057.treat.rf.result$sig.spe$importance.df
srp057.treat.import$features <- rownames(srp057.treat.import)
srp057.treat.import$Median <- apply(srp057.treat.import[,1:100],1, median)

srp057.treat.import.melt <- melt(srp057.treat.import, 
                                id.vars = c("Median", "features"))

srp057.treat.import.summary <- srp057.treat.import %>%
  arrange(desc(Median))


srp057.treat.import.melt <- merge(srp057.treat.import.melt, 
                                 srp057.spe.type, 
                                 by.x = "features", 
                                 by.y = "Species",
                                 all.x = T)

srp057.treat.import.melt %>%
  filter(features %in% srp057.treat.import.summary[1:30,]$features) %>%
  mutate(features = factor(features, 
                           levels = rev(srp057.treat.import.summary[1:30,]$features))) %>%
  ggplot(aes(x = features, y = value, 
             color = type)) + 
  geom_boxplot() +
  coord_flip() +
  scale_color_manual(values = c("#e85a71","#6AAFE6","#0072B2"))

