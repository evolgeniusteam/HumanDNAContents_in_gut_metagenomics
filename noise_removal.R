## prepare metadata
crc.meta <- read.delim("allCRC_metadata.txt", header = T, sep = "\t",as.is = T)
crc.meta$HDC <- crc.meta$humanper *100
rownames(crc.meta) <- crc.meta$Sample_ID

## prepare taxonomic data and metabolic data
load("species.Rdata")
load("path.Rdata")

abund.s.CRC <- lapply(all.abund.list.parse, function(data){y <- subset(data, Kingdom %in% "Bacteria"  & is.na(Strain) & !is.na(Species))[,-c(1:7,9)]; rownames(y) <- y$Species;  y <- data[,-1]/100; return(y)})

pathabund.2 <- pathabund
pathabund.2[["PRJNA447983_cohort1"]] <- pathabund.2[["PRJNA447983_cohort1"]][,subset(crc.meta, Project %in% "PRJNA447983_cohort1")$Sample_ID]
pathabund.2[["PRJEB6070"]] <- pathabund.2[["PRJEB6070"]][,subset(crc.meta, Project %in% "PRJEB6070")$Sample_ID]
pathabund.2[["PRJEB12449"]] <- pathabund.2[["PRJEB12449"]][,subset(crc.meta, Project %in% "PRJEB12449")$Sample_ID]
pathabund.2[["PRJEB7774"]] <- pathabund.2[["PRJEB7774"]][,subset(crc.meta, Project %in% "PRJEB7774")$Sample_ID]

## ---- noise.removal -----##
## ---- remove low abundance data ---- ##
noise.removal <- function(data, percent = 0.01, percent2 = 0.001, top=NULL){
  ##the input should be a data with features as rows and samples as columns.
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

crc.abund.s.filter <- lapply(abund.s.CRC, noise.removal)

## ---- remove stratified info ----##
pathabund.filter <- lapply(pathabund.2, function(x){
  y <- x[-grep(rownames(x), pattern = "|", fixed = T),]
})
lapply(pathabund.filter, dim)
## ---- remove low abundance data ---- ##
pathabund.filter2 <- lapply(pathabund.filter, function(x){
  y <- noise.removal(x, percent = 0, percent2 = 1e-06)
  y <- y[!rownames(y) %in% c("UNMAPPED","UNINTEGRATED"),]
})
lapply(pathabund.filter2, dim)

save(crc.abund.s.filter, pathabund.filter2, file = "filterdata.Rdata")
