## ---- load filtered data ----- ##
load("filterdata.Rdata")
crc.meta <- read.delim("allCRC_metadata.txt", header = T, sep = "\t",as.is = T)

## ---- correlation function ----##
cor.func <- function(data, type = "CRC", project = "PRJEB6070"){
  ## ---- HDC as first column
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
  ## 
  data.cor.sum <- cbind(data.coef.nona, data.P.nona, data.sig)
  names(data.cor.sum) <- c("coef","P_value","SigLevel")
  data.cor.sum[,"Species"] <- rownames(data.cor.sum)
  data.cor.sum[,"Type"] <- type
  data.cor.sum[,"Project"] <- project
  sig.taxa <- subset(data.cor.sum,! SigLevel %in% "")
  return(list(data.cor.sum = data.cor.sum,
              sig.taxa = sig.taxa))
}

## -----species level -----##
crc.abund.s.filter.t <- crc.abund.s.filter
for (pro in names(crc.abund.s.filter)) {
  a <- data.frame(t(crc.abund.s.filter[[pro]]), stringsAsFactors = F);
  a.merge <- cbind(subset(crc.meta, Project %in% pro)[rownames(a),"HDC"], a)
  names(a.merge)[1] <- "HDC"
  crc.abund.s.filter.t[[pro]] <- a.merge
}

crc.abund.s.filt.cor <- list()
for(pro in names(crc.abund.s.filter.t)){
  crc.abund.s.filt.cor[[pro]] <- cor.func(crc.abund.s.filter.t[[pro]], type = "CRC", project = pro)
}
crc.abund.s.filt.sig <- lapply(crc.abund.s.filt.cor, function(data){y <- data$sig.taxa})
all.CRC.cor.sig <- reduce(crc.abund.s.filt.sig, rbind)
all.CRC.cor.posi <- all.CRC.cor.sig %>% filter(coef > 0)
all.CRC.cor.nega <- all.CRC.cor.sig %>%  filter(coef < 0)
delete.spe <- intersect(all.CRC.cor.posi$Species, all.CRC.cor.nega$Species)
all.CRC.cor.sig <- all.CRC.cor.sig %>% filter(!Species %in% delete.spe)

##----species appeared twice or more ------
all.CRC.cor.sig.species <- names(table(all.CRC.cor.sig$Species)[table(all.CRC.cor.sig$Species)>1])

all.CRC.cor.summa <- all.CRC.cor.sig %>% group_by(Species) %>% 
  summarise(median = median(coef)) %>% 
  arrange(median)
## ---- load features which are diagnostic signatures in meta-analysis
study_spe_vector <- read.delim("study_spe_vector" , header = F, as.is = T)
all.CRC.cor.sig$study <- "not include"
all.CRC.cor.sig[which(all.CRC.cor.sig$Species %in% study_spe_vector$V1), "study"] <- "include"
CRC.sigspe.plot <- all.CRC.cor.sig %>%
  filter(Species %in% all.CRC.cor.sig.species) %>%
  mutate(Species = factor(Species, levels = all.CRC.cor.summa$Species)) %>%
  mutate(study = factor(study, levels = c("not include","include"))) %>%
  ggplot(aes(x = Species, y = coef))  +
  geom_point(size = 4, alpha = .7, aes(shape = Project, color = study)) +
  scale_shape_manual(values = c(0,1,3,4,5,20,8)) +
  scale_color_manual(values = c("black","darkred")) +
  coord_flip()  + theme_classic()+
  geom_hline(aes(yintercept = 0), colour="#990000", linetype="dashed", alpha = .8) + 
  theme(axis.text = element_text(size = 11),
        plot.title = element_text(vjust = 0.5,hjust = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour = "lightgrey"))
        
##----pathways level ----------##
pathabund.filter.t <- pathabund.filter2
for (pro in names(pathabund.filter2)) {
  x <- pathabund.filter2[[pro]]
  x$path <- rownames(x)
  a <- data.frame(t(x), stringsAsFactors = F)
  names(a) <- a["path",]
  a <- a[!rownames(a) %in% "path",]
  a.merge <- cbind(subset(crc.meta, Project %in% pro)[rownames(a),"HDC"], a)
  names(a.merge)[1] <- "HDC"
  pathabund.filter.t[[pro]] <- a.merge
}

pathabund.corresult <- list()
for(pro in names(pathabund.filter.t)){
  pathabund.corresult[[pro]] <- cor.func(pathabund.filter.t[[pro]], type = "CRC", project = pro)
}
pathabund.sig.list <- lapply(pathabund.corresult, function(data){y <- data$sig.taxa})
pathabund.cor.sig <- reduce(pathabund.sig.list, rbind)
names(pathabund.cor.sig)[4] <- "pathways"
path.cor.posi <- pathabund.cor.sig %>% filter(coef > 0)
path.cor.nega <- pathabund.cor.sig %>%  filter(coef < 0)
delete.path <- intersect(path.cor.posi$pathways, path.cor.nega$pathways)
pathabund.cor.sig <- pathabund.cor.sig %>% filter(!pathways %in% delete.path)
pathabund.cor.sig$Abb <- unlist(lapply(strsplit(pathabund.cor.sig$pathways, split = ":", fixed = T), function(data){y <- data[1]}))

## ----- HDC correlated pathways appeared three times or more ------
pathabund.cor.sig.2 <- subset(pathabund.cor.sig, Abb %in% names(table(pathabund.cor.sig$Abb)[table(pathabund.cor.sig$Abb)>2]))

## -----
pathabund.cor.summa <- pathabund.cor.sig %>% group_by(Abb) %>% 
  summarise(median = median(coef)) %>% 
  arrange(median)
study_pathways_vector <- read.delim("study_pathways_vector" , header = F, as.is = T)

pathabund.cor.sig$study <- "not include"
pathabund.cor.sig[which(pathabund.cor.sig$Abb %in% study_pathways_vector$V1), "study"] <- "include"
CRC.sigpath.plot <- pathabund.cor.sig %>%
  filter(Abb %in% pathabund.cor.sig.2$Abb) %>%
  mutate(Abb = factor(Abb, levels = pathabund.cor.summa$Abb)) %>%
  mutate(study = factor(study, levels = c("not include","include"))) %>%
  ggplot(aes(x = Abb, y = coef))  +
  geom_point(size = 4, alpha = .7, aes(shape = Project, color = study)) +
  scale_shape_manual(values = c(0,1,3,4,5,20,8)) +
  scale_color_manual(values = c("black","darkred")) +
  coord_flip()  + theme_classic()+
  geom_hline(aes(yintercept = 0), colour="#990000", linetype="dashed", alpha = .8) + 
  theme(axis.text = element_text(size = 11),
        plot.title = element_text(vjust = 0.5,hjust = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour = "lightgrey"))

## ---- HDC correlated features -----##

sig.spe <- lapply(crc.abund.s.filt.sig, rownames)
sig.spe.vector <- Reduce(c, sig.spe)
sig.spe.vector <- unique(sig.spe.vector)

sig.path <- lapply(pathabund.sig.list, rownames)
sig.path.vector <- Reduce(c, sig.path)
sig.path.vector <- unique(sig.path.vector)
