
library(rwrfhydro)
library(data.table)
library(ggplot2)
library(plyr)
library(sensitivity)
library(randtoolbox)

#########################################################
# SETUP
#########################################################

source("calib_utils.R")

# Metrics
metrics <- c("obj", "nse", "nselog", "cor", "rmse", "bias", "kge", "fdc")


#########################################################
# MAIN CODE
#########################################################

# First run so need to initialize
ReadNamelist("namelist.calib")
cyclecount <- 0

# Setup value lists from paramBnds
message("Setup lists")

if (SA_method == "DELSA") {
  
  # Setup value lists from paramBnds
  par.ranges <- list()
  for(i in 1:nrow(paramBnds)) {
    par.ranges[[paramBnds[i, "param"]]] <- unname(unlist(paramBnds[i, c("min", "max")]))
  }
  X0 <- sensitivity::parameterSets(par.ranges, samples = SA_sample_size, method = c(SA_par_gen_method))
  varprior = sapply(par.ranges, diff)^2/12
  
  X = do.call(rbind, lapply(1:ncol(X0), function(i) {
    X2i = X0
    X2i[, i] = X2i[, i] * SA_perturb
    X2i
  }))
  
  x_all = as.data.frame(rbind(X0, X))
  names(x_all) <- names(par.ranges)
  
}
# add the index as the first columns, since it existed in Aubrey s way of coding 
x_all <- cbind.data.frame(id = c(1:nrow(x_all)), x_all)

# Initialize parameter archive DF
message("Initialize parameter archive")
x_archive <- as.data.frame(matrix(, nrow=1, ncol=ncol(x_all)+length(metrics)))
names(x_archive) <- c(names(x_all), metrics)

# Output parameter set
message("Output parameter set")
write.table(x_all, file="params_new.txt", row.names=FALSE, sep=" ")

# Save and exit
save.image("proj_data_SENS.Rdata")
quit("no")

