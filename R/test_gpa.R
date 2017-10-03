

alldat <- readRDS("/Users/bbuchsbaum/analysis/hyper/ventral_surface/ROI_11121_alldat.RDS")

vmat <- do.call(cbind, lapply(alldat$vdat, function(x) x$mat))
svmat <- scale(vmat)

nmat <- do.call(cbind, lapply(alldat$ndat, function(x) x$mat))
snmat <- scale(nmat)

both_mat <- rbind(svmat,snmat)
groups <- sapply(vdat, function(x) ncol(x$mat))
gpa_res <- GPA(as.data.frame(both_mat), group=groups)
