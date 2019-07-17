##  ------ wilcox test ------- ##
crc.meta <- read.delim("allCRC_metadata.txt", header = T, sep = "\t",as.is = T)
rownames(crc.meta) <- crc.meta$Sample_ID
load("filterdata.Rdata")

spe.wt.list <- list()
for (pro in names(crc.abund.s.filter)) {
  x.t <- data.frame(t(crc.abund.s.filter[[pro]]), stringsAsFactors = F);
  x.t$Group <- subset(crc.meta, Project %in% pro)[rownames(a),"Group"]
  x.t.result <- data.frame(Species = rownames(crc.abund.s.filter[[pro]]),
                           p.value = NA)
  for (i in 1:(ncol(x.t)-1)) {
    a <- wilcox.test(x.t[,i] ~ x.t$Group, x.t)
    x.t.result[i, "p.value"] <- a$p.value
  }
  spe.wt.list[[pro]] <- x.t.result
}
spe.wt.list.spe <- lapply(spe.wt.list, function(x){y <- subset(x, p.value<0.05)$Species; y <- as.character(y)})
spe.wt.list.spe.vector <- Reduce(c, spe.wt.list.spe)
spe.wt.list.spe.vector <- unique(spe.wt.list.spe.vector)

spe.wt.data <- lapply(crc.abund.s.filter, function(x){
  y <- x[spe.wt.list.spe.vector, ];
  y[is.na(y)] <- 0;
  rownames(y) <- spe.wt.list.spe.vector;
  return(y)
})

spe.hdc.data <- lapply(crc.abund.s.filter, function(x){
  y <- x[sig.spe.vector, ];
  y[is.na(y)] <- 0;
  rownames(y) <- sig.spe.vector;
  return(y)
})

pathabund.wt.list <- list()
for (project in names(pathabund.filter2)) {
  x <- pathabund.filter2[[project]]
  x.meta <- subset(crc.meta, Project %in% project)
  x.t <-data.frame( t(x), stringsAsFactors = F)
  x.meta <- x.meta[rownames(x.t),]
  x.t$Group <- x.meta$Group
  x.t.result <- data.frame(pathways = rownames(x),
                           p.value = NA)
  for (i in 1:(ncol(x.t)-1)) {
    a <- wilcox.test(x.t[,i] ~ x.t$Group, x.t)
    x.t.result[i, "p.value"] <- a$p.value
  }
  pathabund.wt.list[[project]] <- x.t.result
}
pathabund.wt.list.path <- lapply(pathabund.wt.list, function(x){y <- subset(x, p.value<0.05)$pathway; y <- as.character(y)})
pathabund.wt.list.path.vector <- Reduce(c, pathabund.wt.list.path)
pathabund.wt.list.path.vector  <- unique(pathabund.wt.list.path.vector )

path.wt.data <- lapply(pathabund.filter2, function(x){
  y <- x[pathabund.wt.list.path.vector, ];
  y[is.na(y)] <- 0;
  rownames(y) <- pathabund.wt.list.path.vector;
  return(y)
})

path.hdc.data <- lapply(pathabund.filter2, function(x){
  y <- x[sig.path.vector, ];
  y[is.na(y)] <- 0;
  rownames(y) <- sig.path.vector;
  return(y)
})

## -----
all.features <- lapply(crc.abund.s.filter, rownames)
all.features.vector <- unique(Reduce(c, all.features))
spe.data <- lapply(crc.abund.s.filter, function(x){
  y <- x[all.features.vector,]
  y[is.na(y)] <- 0;
  rownames(y) <- all.features.vector;
  return(y)
})

all.features <- lapply(pathabund.filter2, rownames)
all.features.vector <- unique(Reduce(c, all.features))
pathways.data <- lapply(pathabund.filter2, function(x){
  y <- x[all.features.vector,]
  y[is.na(y)] <- 0;
  rownames(y) <- all.features.vector;
  return(y)
})


data.list <- list(spe.data = spe.data,
                  pathways.data  = pathways.data ,
                  spe.wt.data = spe.wt.data,
                  spe.hdc.data = spe.hdc.data,
                  path.wt.data = path.wt.data,
                  path.hdc.data = path.hdc.data)
save(data.list, file = "data.Rdata")
