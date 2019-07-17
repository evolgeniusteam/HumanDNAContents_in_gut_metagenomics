## ----- example -----##
cl <- makeCluster(7)
registerDoParallel(cl)
study <- unique(crc.meta$Project)
all.spe.cross <- foreach(project = study) %dopar% cross.models(project = project, metadata=crc.meta, rl.data = spe.data.df)
all.spe.lodo <- foreach(project = study) %dopar% lodo.models(project = project, metadata=crc.meta, rl.data = spe.data.df)
stopCluster(cl)

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
