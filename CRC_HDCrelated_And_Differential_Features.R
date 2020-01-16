library(dplyr)
library(Hmisc)
library(reshape2)
load("CRC_basicinfo.Rdata")

## load: crc.abund.s.list, metadata, species information

crc.abund.spe.filter <- lapply(crc.abund.s.list, 
                               function(data){ 
                                 y <- data[, colnames(data) %in% rownames(crc.meta)]; 
                                 return(y)})

noise.removal <- function(data, percent = 0.01, percent2 = 0.001, top=NULL){
  ## rows are species, cols are samples
  Matrix <- data
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  rowmax <- apply(Matrix_1, 1, max)
  bigones2 <- rowmax >= percent2 
  print(percent2)
  Matrix_2 <- Matrix_1[bigones2,]
  return(Matrix_2)
}

## ------remove species with low abundance
crc.abund.spe.filter <- lapply(crc.abund.spe.filter, 
                               noise.removal)
                               
## ------remove species which were not bacteria
crc.abund.spe.final <- list()
for (pro in names(crc.abund.spe.filter)) {
  crc.abund.spe.final[[pro]] <- crc.abund.spe.filter[[pro]][rownames(crc.abund.spe.filter[[pro]]) %in% all.spe.df.merge$Species, ]
}

crc.abund.merge.list <- list()
for (pro in names(crc.abund.spe.final)) {
  a <- data.frame(t(crc.abund.spe.final[[pro]]),
                  stringsAsFactors = F)
                  
  crc.abund.merge.list[[pro]] <- cbind(crc.meta[rownames(a), "HDC"],
                                       a)
  names(crc.abund.merge.list[[pro]])[1] <- "HDC"
  crc.abund.merge.list[[pro]] <- crc.abund.merge.list[[pro]][!rownames(crc.abund.merge.list[[pro]]) %in% "features", ]
  for (i in 1:ncol(crc.abund.merge.list[[pro]])) {
    crc.abund.merge.list[[pro]][,i] <- as.numeric(crc.abund.merge.list[[pro]][,i])
  }
}

## -------HDC correlation -----#

## ---- correlation function ----##
cor.func <- function(data, type = "CD", project = "PRJNA389280"){
  ##data, first column is humanper 
  cut_sig <- function(p) {
    out <- cut(p, breaks = c(0, 0.01,0.05,1), include.lowest = T, 
               labels = c("**", "*", ""))
    return(out)
  }
  ##-----require data, first column is humanper----
  data.cor <- rcorr(as.matrix(data), type = "spearman")
  data.coef <- data.cor$r[1,-1]
  data.P <- data.cor$P[1,-1]
  data.P.nona <- na.omit(data.P)
  data.coef.nona <- na.omit(data.coef)
  data.sig <- data.frame(apply(data.frame(data.P.nona), 2, cut_sig))
  rownames(data.sig) <- names(data.P.nona)
  ##samples in week 0
  ##coef results bewtween humanDNA% and taxa abundance
  data.cor.sum <- cbind(data.coef.nona, data.P.nona, data.sig)
  names(data.cor.sum) <- c("coef","P_value","SigLevel")
  data.cor.sum[,"Species"] <- rownames(data.cor.sum)
  data.cor.sum[,"Type"] <- type
  data.cor.sum[,"Project"] <- project
  sig.taxa <- subset(data.cor.sum,! SigLevel %in% "")
  return(list(data.cor.sum = data.cor.sum,
              sig.taxa = sig.taxa))
}


crc.abund.s.filt.cor <- list()
for (pro in names(crc.abund.merge.list)) {
  crc.abund.s.filt.cor[[pro]] <- cor.func(crc.abund.merge.list[[pro]], type = "CRC", project = pro)
}

crc.cor.df <- lapply(crc.abund.s.filt.cor, function(data){
  y <- data$data.cor.sum;
  return(y);}) %>%reduce(rbind)

cor.data.sig <- subset(crc.cor.df, P_value<0.05)

## -----delete species whose correlation didn't keep a consistent trend
cor.data.posi <- subset(cor.data.sig, coef>0)
cor.data.nega <- subset(cor.data.sig, coef<0)
delete.spe <- intersect(cor.data.nega$Species, cor.data.posi$Species)
cor.data.sig <- subset(cor.data.sig, !Species %in% delete.spe)

cor.data.sig.sum <- cor.data.sig %>% group_by(Species, Type) %>% summarise(count= n())
cor.data.sig.sum.2 <- dcast(cor.data.sig.sum, Species ~Type)

## -----26 HDC-related species -------##
crc.sig.spe <- subset(cor.data.sig.sum.2, CRC>1 )$Species

# ------Differential species -----------##
crc.kw <- list()

for (project in names(crc.abund.spe.final)) {
  y <- crc.abund.spe.final[[project]]
  x <- data.frame(t(y), stringsAsFactors = F)
  print(rownames(x))
  x.meta <- subset(crc.meta, Sample_ID %in% rownames(x))
  rownames(x.meta) <- x.meta$Sample_ID
  x <- x[rownames(x.meta),]
  x$Group <- x.meta$Group
  x$Group <- factor(x$Group, levels = c("CTR","CRC"))
  x.result <- data.frame(features = names(x),
                         p.value = NA,
                         greater= NA,
                         less = NA)
  features <- names(x)
  rownames(x.result) <- features
  features <- features[!features %in% "Group"]
  for (i in features) {
    a <- wilcox.test(x[,i] ~ x$Group, x)
    x.result[i, "p.value"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "greater")
    x.result[i, "greater"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "less")
    x.result[i, "less"] <- a$p.value
  }
  x.result$Project <- project
  crc.kw[[project]] <- x.result
}
## -----adjusted P-value --------------------##
crc.kw.fdr <- crc.kw %>% reduce(rbind) %>% filter(!features %in% "Group")
crc.kw.fdr <- crc.kw.fdr %>% group_by(Project) %>% arrange(p.value) %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr"))

crc.kw.fdr$trend <- "Control"
crc.kw.fdr[which(crc.kw.fdr$less<0.05), "trend"] <- "Case"

crc.kw.sum <- crc.kw.fdr %>% filter(fdr <0.05) %>% group_by(features, trend) %>% summarise(count = n())

## -----delete species whose alteration didn't keep a consistent trend
crc.kw.CTR <- subset(crc.kw.sum, trend %in% "Control")
crc.kw.CRC <- subset(crc.kw.sum, trend %in% "Case")
delete.spe <- intersect(crc.kw.CTR$features, crc.kw.CRC$features)
##rm(crc.kw.sum.2)
crc.kw.sum.3 <- subset(crc.kw.sum, count >1)
crc.kw.sum.3 <- data.frame(crc.kw.sum.3, stringsAsFactors = F)
crc.kw.sum.3$features <- as.character(crc.kw.sum.3$features)
crc.kw.sum.3$Disease <- "CRC"
crc.kw.spe <- crc.kw.sum.3$features ## 16 differential species -----------

## -------overall species data -----------##
abund.s.list <- lapply(crc.abund.spe.final, 
                       function(data){
  data$features <- rownames(data);
  return(data)
})
abund.s.list.df <- abund.s.list %>% reduce(merge, 
                                           by = "features", 
                                           all =T)
abund.s.list.df[is.na(abund.s.list.df)] <- 0
rownames(abund.s.list.df) <- abund.s.list.df$features
abund.s.list.df <- abund.s.list.df[, -1]

## ----312 species , 711 samples -----------##
spe.data <- data.frame(t(abund.s.list.df), 
                       stringsAsFactors = F)
spe.data[, c("HDC", "Group")] <- crc.meta[rownames(spe.data), c("HDC", "Group")]
spe.data$Group <- factor(spe.data$Group, 
                         levels = c("CTR","CRC"))
                         

## ---------create species list -----------##
spe.list <- list(all.spe.data = spe.data, 
                 kw.spe.data = spe.data[, names(spe.data) %in% c(crc.kw.spe, "HDC", "Group")],
                 sig.spe.data = spe.data[,names(spe.data) %in%  c(crc.sig.spe, "HDC", "Group")])



##---------Metabolic level------------------------
##---------remove pathways with low abundances
noise.removal.path <- function(data, percent = 0.1, top=NULL){
  Matrix <- data
  zero.num <- apply(Matrix, 1, function(data){y <- length(data[data ==0])})
  col.num <- ncol(data)
  bigones <- zero.num<=0.15*col.num
  Matrix_1 <- Matrix[bigones,]
  
  rowmax <- apply(Matrix_1, 1, max)
  bigones2 <- rowmax >= percent 
  Matrix_2 <- Matrix_1[bigones2,]
  
  return(Matrix_2)
}

pathabund.filter.adenoma.2 <- lapply(pathabund.filter.adenoma[1:7],
                                     function(data){
                                       y <- data[, names(data) %in% rownames(crc.meta)];
                                       y <- y[!rownames(y) %in% c("UNMAPPED","UNINTEGRATED"),]
                                       y.2 <- noise.removal.path(y, percent = 1e-06);
                                       return(y.2)
                                     }) 
for (names in names(pathabund.filter.adenoma.2)) {
  rownames(pathabund.filter.adenoma.2[[names]]) <- unlist(lapply(strsplit(rownames(pathabund.filter.adenoma.2[[names]]), 
                                                                          split = ":", fixed = T), 
                                                                 function(y){x <- y[[1]]}))
  rownames(pathabund.filter.adenoma.2[[names]]) <- gsub(rownames(pathabund.filter.adenoma.2[[names]]), 
                              pattern = "-", replacement = "_", fixed = T)
  rownames(pathabund.filter.adenoma.2[[names]]) <- gsub(rownames(pathabund.filter.adenoma.2[[names]]), 
                              pattern = "+", replacement = "_", fixed = T)
  
  pathabund.filter.adenoma.2[[names]]$pathway <- rownames(pathabund.filter.adenoma.2[[names]])
}

## -------overall pathways data -----------##
pathabund.filter.adenoma2.df <- pathabund.filter.adenoma.2 %>% 
  reduce(merge, by = "pathway", all = T)
pathabund.filter.adenoma2.df[is.na(pathabund.filter.adenoma2.df)] <- 0
rownames(pathabund.filter.adenoma2.df) <- paste0("pathways_",
                                                 pathabund.filter.adenoma2.df$pathway)
pathabund.filter.adenoma2.df <- pathabund.filter.adenoma2.df[,-1]
path.data <- data.frame(t(pathabund.filter.adenoma2.df),
                        stringsAsFactors = F)
path.data[, c("HDC", "Group")] <- crc.meta[rownames(path.data), c("HDC", "Group")]
path.data$Group <- factor(path.data$Group, 
                         levels = c("CTR","CRC"))

for (i in 1:346) {
  path.data[,i] <- as.numeric(path.data[,i])
}

##--------HDC related pathways----------------##
crc.path.abund.merge.list <- list()
for (pro in names(pathabund.filter.adenoma.2)) {
  a <- data.frame(t(pathabund.filter.adenoma.2[[pro]]),
                  stringsAsFactors = F)
  crc.path.abund.merge.list[[pro]] <- cbind(crc.meta[rownames(a), "HDC"],
                                       a)
  names(crc.path.abund.merge.list[[pro]])[1] <- "HDC"
  crc.path.abund.merge.list[[pro]] <- crc.path.abund.merge.list[[pro]][!rownames(crc.path.abund.merge.list[[pro]]) %in% "pathway", ]
  for (i in 1:ncol(crc.path.abund.merge.list[[pro]])) {
    crc.path.abund.merge.list[[pro]][,i] <- as.numeric(crc.path.abund.merge.list[[pro]][,i])
  }
}
names(crc.path.abund.merge.list)[5] <- "PRJNA447983_cohort1"

## -------HDC correlation -----#
crc.abund.path.filt.cor <- list()
for (pro in names(crc.path.abund.merge.list)) {
  crc.abund.path.filt.cor[[pro]] <- cor.func(crc.path.abund.merge.list[[pro]], 
                                             type = "CRC", 
                                             project = pro)
}

crc.path.cor.df <- lapply(crc.abund.path.filt.cor, function(data){
  y <- data$data.cor.sum;
  return(y);}) %>%reduce(rbind)

crc.cor.path.sig <- subset(crc.path.cor.df, P_value<0.05)

## -----delete pathways whose correlation didn't keep a consistent trend
cor.path.posi <- subset(crc.cor.path.sig, coef>0)
cor.path.nega <- subset(crc.cor.path.sig, coef<0)
delete.spe <- intersect(cor.path.posi$Species, cor.path.nega$Species)
crc.cor.path.sig <- subset(crc.cor.path.sig, !Species %in% delete.spe)

crc.path.sig.sum <- crc.cor.path.sig %>% 
  group_by(Species, Type) %>% summarise(count= n())
crc.path.sig.sum.2 <- dcast(crc.path.sig.sum, Species ~Type)


## -------
crc.sig.path_3 <- subset(cor.path.sig.sum.2, CRC>2)$Species ## 11个
crc.sig.path <- subset(cor.path.sig.sum.2, CRC>1)$Species ## 40个

## --------53 differential pathways ---------##
crc.path.kw<- list()

for (project in names(pathabund.filter.adenoma.2)) {
  y <- pathabund.filter.adenoma.2[[project]]
  x <- data.frame(t(y), stringsAsFactors = F)
  print(rownames(x))
  x.meta <- subset(crc.meta, Sample_ID %in% rownames(x))
  rownames(x.meta) <- x.meta$Sample_ID
  x <- x[rownames(x.meta),]
  for (i in 1:ncol(x)) {
    x[,i] <- as.numeric(x[,i])
  }
  x$Group <- x.meta$Group
  x$Group <- factor(x$Group, levels = c("CTR","CRC"))
  x.result <- data.frame(features = names(x),
                         p.value = NA,
                         greater= NA,
                         less = NA)
  features <- names(x)
  rownames(x.result) <- features
  features <- features[!features %in% "Group"]
  for (i in features) {
    a <- wilcox.test(x[,i] ~ x$Group, x)
    x.result[i, "p.value"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "greater")
    x.result[i, "greater"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "less")
    x.result[i, "less"] <- a$p.value
  }
  x.result$Project <- project
  crc.path.kw[[project]] <- x.result
}

##-------------adjusted P-value-----------------##
path.kw.fdr <- crc.path.kw %>% reduce(rbind) %>% filter(!features %in% "Group")
path.kw.fdr <- path.kw.fdr %>% group_by(Project) %>% arrange(p.value) %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr"))

path.kw.fdr$trend <- "Control"
path.kw.fdr[which(path.kw.fdr$less<0.05), "trend"] <- "Case"

path.kw.sum <- path.kw.fdr %>% 
  filter(fdr <0.05) %>% 
  group_by(features, trend) %>% 
  summarise(count = n())

## -----delete pathways whose alteration didn't keep a consistent trend

path.kw.CTR <- subset(path.kw.sum, trend %in% "Control")
path.kw.CRC <- subset(path.kw.sum, trend %in% "Case")
delete.spe <- intersect(path.kw.CTR$features, 
                        path.kw.CRC$features)

path.kw.sum.2 <- subset(path.kw.sum,
                        ! features %in% delete.spe)

path.kw.sum.3 <- subset(path.kw.sum.2, count >1)
path.kw.sum.3 <- data.frame(path.kw.sum.3, 
                            stringsAsFactors = F)
path.kw.sum.3$features <- as.character(path.kw.sum.3$features)
path.kw.sum.3$Disease <- "CRC"

crc.kw.path <- path.kw.sum.3$features

kw.path <- paste0("pathways_", crc.kw.path)
sig.path <- paste0("pathways_", crc.sig.path)


## ---------create pathways list -----------##
path.list <- list(all.path.data = path.data, 
                 kw.path.data = path.data[, names(path.data) %in% c(kw.path, "HDC", "Group")],
                 sig.path.data = path.data[,names(path.data) %in%  c(sig.path, "HDC", "Group")])


save.image("article1_basic.RData")

## --------for models ---------------------##
save(spe.list, crc.meta, 
     file = "spe_list.Rdata")
save(path.list, crc.meta, 
     file = "path_list.Rdata")
