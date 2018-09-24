library(devtools)
library(hdf5r)
library(tibble)
library(memoise)

bad <- c(1015, 1017,1029, 1023)
base_path <- "/Users/bbuchsbaum/analysis/hyper/all_surface_h5/"

sids <- scan(paste0(base_path, "/sids"), "")
sids <- sids[!sids %in% bad]
fnames <- paste0(base_path, paste0(sids, "_surface_data.h5"))

dfsim <- read.table(paste0(base_path, "sim_summary.txt"), header=TRUE)

dfsim$rpc1 <- .25 + (dfsim$pc1 - min(dfsim$pc1)) / ((max(dfsim$pc1) - min(dfsim$pc1)) * 2)
simtab <- as.list(dfsim$rpc1)
names(simtab) <- dfsim$label

distlab <- function(x1,x2) {
  if (substr(x1,1,3) != substr(x2,1,3)) {
    1
  } else if (x1 == x2) {
    0
  } else {
    ad <- adist(x1,x2)[,]
    if (ad == 1) {
      label <- substr(x1, 1, nchar(x2)-2)
      ret <- 1 - simtab[[label]]
    } else {
      1
    }
  }
}


simlab <- function(x1,x2) {
  if (substr(x1,1,3) != substr(x2,1,3)) {
    0
  } else if (x1 == x2) {
    1
  } else {
    ad <- adist(x1,x2)[,]
    if (ad == 1) {
      label <- substr(x1, 1, nchar(x2)-2)
      ret <- simtab[[label]]
    } else {
      0
    }
  }
}

#simlab <- memoise(simlab)

simfun <- function(x1, x2) {

  sapply(1:length(x1), function(i) {
    simlab(x1[i], x2[i])
  })
}

distfun <- function(x1, x2) {

  sapply(1:length(x1), function(i) {
    distlab(x1[i], x2[i])
  })
}


save_roi <- function(rnum) {

  Xl <- do.call(rbind, lapply(1:length(sids), function(i) {
    fn <- fnames[i]
    s <- sids[i]
    fh5 <- H5File$new(fn, mode="r")
    labs <- fh5[[paste0(rnum, "/nback/labels")]][]
    mat <- fh5[[paste0(rnum, "/nback/data")]][,]
    smat <- scale(t(scale(t(mat))))
    ret <- tibble(sid=s, labs=list(labs), X=list(smat))
    fh5$close()
    ret
  }))

  saveRDS(Xl, paste0("~/Dropbox/code/neuroca/testdata/hyper_", rnum, "_nback.rds"))

}


Xs <- Xl[["X"]]
labs <- Xl[["labs"]]
