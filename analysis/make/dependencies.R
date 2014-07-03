#!/usr/bin/env Rscript


pkgs.cran <- c("ape", "phytools", "geiger", "phylolm", "foreach",
               "doMC", "MASS", "nlme", "ggplot2", "reshape2", "RColorBrewer",
               "plyr", "scales", "lattice", "knitr", "markdown")

pkgs.github <- c("richfitz/sowsear"="0.1-1")


pkgs <- installed.packages()

## eventually it will be better to monitor what version we are using
## will add this later
to.install <- setdiff(pkgs.cran, rownames(pkgs))

## installation message
msg.install <- (if (length(to.install) == 0) character(0) else
                paste("installing:", paste(to.install, collapse=", ")))

if (length(msg.install) > 0){
    message(paste(msg.install, collapse="\n"))
    install.packages(to.install)
}


## github
pkgs.github.short <- sub("^.+/", "", names(pkgs.github))

tr <- structure(names(pkgs.github), names=pkgs.github.short)

to.install <- setdiff(pkgs.github.short, rownames(pkgs))
installed <- intersect(pkgs.github.short, rownames(pkgs))

to.upgrade <-
  installed[numeric_version(pkgs[installed, "Version"]) <
            numeric_version(pkgs.github[tr[installed]])]

msg.install <- (if (length(to.install) == 0) character(0) else
                paste("installing:", paste(to.install, collapse=", ")))
msg.upgrade <- (if (length(to.upgrade) == 0) character(0) else
                paste("upgrading:", paste(to.upgrade, collapse=", ")))
msg <- c(msg.install, msg.upgrade)
if (length(msg) > 0) {
  message(paste(msg, collapse="\n"))
  install.packages("devtools")
  library(devtools)
  for (pkg in tr[c(to.install, to.upgrade)]) {
    install_github(pkg)
  }
}
