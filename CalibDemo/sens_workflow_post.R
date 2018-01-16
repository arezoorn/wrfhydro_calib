library(rwrfhydro)
library(data.table)
library(ggplot2)
library(plyr)
library(boot)
#########################################################
# SETUP
#########################################################

source("calib_utils.R")

# Multi-core
parallelFlag <- TRUE
ncores <- 32
if (parallelFlag && ncores>1) {
  library(doParallel)
  cl <- makeForkCluster(ncores)
  registerDoParallel(cl)
}


#########################################################
# MAIN CODE
#########################################################

if (file.exists("proj_data_SENS.Rdata")) { 
  load("proj_data_SENS.Rdata")
} else {
  message("No proj_data_SENS.Rdata file found")
}

# Load obs so we have them for next iteration
load(obsFile)
obsDT <- obsDT[!is.na(obs),]

# convrt the hourly obs to daily obs
obsDT$Date <- rwrfhydro::CalcDateTrunc(obsDT$POSIXct)
setkey(obsDT, Date)
obsDT.d <- obsDT[, list(obs = mean(obs, na.rm = TRUE)), by = "Date"]

# Find the index of the gage
rtLink <- ReadRouteLink(rtlinkFile)
rtLink <- data.table(rtLink)
linkId <- which(trimws(rtLink$gages) %in% siteId)

# Initialize chrtout
if (!exists("chrt.d.all")) chrt.d.all <- data.table()
if (!exists("chrt.h.all")) chrt.h.all <- data.table()

for (cyclecount in 1:nrow(x_all)) {
  # Read model out and calculate performance metric
  outPath <- paste0(runDir, "/SENS_RESULTS/OUTPUT", cyclecount)
  print(outPath)
  
  # Read files
  message("Reading model out files.")
  system.time({
    filesList <- list.files(path = outPath,
                            pattern = glob2rx("*.CHRTOUT_DOMAIN*"),
                            full.names = TRUE)
    filesListDate <- as.POSIXct(unlist(plyr::llply(strsplit(basename(filesList),"[.]"), '[',1)), format = "%Y%m%d%H%M", tz = "UTC")
    whFiles <- which(filesListDate >= startDate)
    filesList <- filesList[whFiles]
    if (length(filesList) == 0) stop("No matching files in specified directory.")
    chrt <- as.data.table(plyr::ldply(filesList, ReadChFile, linkId, .parallel = parallelFlag))
  })

  # add the chrt data to the chrt.h.all and calculate the stats for the hourly time step
  chrt[, site_no := siteId]
  setkey(chrt, "site_no", "POSIXct")
  setkey(obsDT, "site_no", "POSIXct")
  chrt.h <- merge(chrt, obsDT, by=c("site_no", "POSIXct"), all.x=FALSE, all.y=FALSE)
  chrt.h$id <- cyclecount
# chrt.h$tag <- x_all$tag[cyclecount] We do not have any tag anymore
  chrt.h.all <- rbindlist(list(chrt.h.all, chrt.h))

  # Assess performance for the hourly flow
  F_new <- objFn(chrt.h$q_cms, chrt.h$obs)
  statNse <- rwrfhydro::Nse(chrt.h$q_cms, chrt.h$obs)
  statNseLog <- rwrfhydro::NseLog(chrt.h$q_cms, chrt.h$obs)
  statCor <- cor(chrt.h$q_cms, chrt.h$obs)
  statRmse <- rwrfhydro::Rmse(chrt.h$q_cms, chrt.h$obs)
  statBias <- sum(chrt.h$q_cms - chrt.h$obs, na.rm=TRUE)/sum(chrt.h$obs, na.rm=TRUE) * 100
  statKge <- Kge(chrt.h$q_cms, chrt.h$obs)
  chrt.h <- CalcFdc(chrt.h, "q_cms")
  chrt.h <- CalcFdc(chrt.h, "obs")
#  statFdc <- integrate(splinefun(as.data.frame(chrt.h)[,"q_cms.fdc"], as.data.frame(chrt.h)[,"q_cms"], method='natural'), 0.05, 0.95, subdivisions=1000)$value
  statFdc <- NA

  # Archive results
  x_archive_h[cyclecount,] <- c(x_all[cyclecount,], F_new, statNse, statNseLog, statCor, statRmse, statBias, statKge, statFdc)
  print(x_archive_h[cyclecount,])

########################## ####### DAILY CALCULATIONS ###################################################### 
  # Convert to daily
  chrt.d <- Convert2Daily(chrt)
  obsDT.d[, site_no := siteId]
  chrt.d[, site_no := siteId]
  # Merge
  setkey(chrt.d, "site_no", "Date")
  setkey(obsDT.d, "site_no", "Date")
  chrt.d <- merge(chrt.d, obsDT.d,  by=c("site_no", "Date"), all.x=FALSE, all.y=FALSE)
  chrt.d$id <- cyclecount
# chrt.d$tag <- x_all$tag[cyclecount] We do not have tag here.
  chrt.d.all <- rbindlist(list(chrt.d.all, chrt.d))
  
  # Assess performance
  F_new <- objFn(chrt.d$q_cms, chrt.d$obs)
  statNse <- rwrfhydro::Nse(chrt.d$q_cms, chrt.d$obs)
  statNseLog <- rwrfhydro::NseLog(chrt.d$q_cms, chrt.d$obs)
  statCor <- cor(chrt.d$q_cms, chrt.d$obs)
  statRmse <- rwrfhydro::Rmse(chrt.d$q_cms, chrt.d$obs)
  statBias <- sum(chrt.d$q_cms - chrt.d$obs, na.rm=TRUE)/sum(chrt.d$obs, na.rm=TRUE) * 100
  statKge <- Kge(chrt.d$q_cms, chrt.d$obs)
  chrt.d <- CalcFdc(chrt.d, "q_cms")
  chrt.d <- CalcFdc(chrt.d, "obs")
  statFdc <- integrate(splinefun(as.data.frame(chrt.d)[,"q_cms.fdc"], as.data.frame(chrt.d)[,"q_cms"], method='natural'), 0.05, 0.95, subdivisions=1000)$value
  
  # Archive results
  x_archive[cyclecount,] <- c(x_all[cyclecount,], F_new, statNse, statNseLog, statCor, statRmse, statBias, statKge, statFdc)
  print(x_archive[cyclecount,])
  
}

# Interim save
save.image("proj_data_SENS.Rdata")

################################ DELSA Calculations for each Metric at both hourly and daily time step

if (SA_method == "DELSA") {
  delsaFirst <- list()
  for (timeStep in c("hourly", "daily")) {
  for (metric in metrics)  {
    if (timeStep == "daily") x <- list(y = x_archive[, metric], X0 = X0, X = rbind(X0, X), varprior = varprior)
    if (timeStep == "hourly") x <- list(y = x_archive_h[, metric], X0 = X0, X = rbind(X0, X), varprior = varprior)
    
    id <- deparse(substitute(x))
    
    Kpar = ncol(x$X0)
    Nsamp = nrow(x$X0)
    vartot = rep(0, Nsamp)
    delsafirst = deriv = varfir = matrix(NA, ncol = Kpar, nrow = Nsamp)
    out <- as.numeric(x$y)
    for (rsamp in 1:Nsamp) {
      for (jpar in 1:Kpar) {
        idx.pert = Nsamp * jpar + rsamp
        deriv[rsamp, jpar] = (out[idx.pert] - out[rsamp])/(x$X[idx.pert,
                                                               jpar] - x$X[rsamp, jpar])
        varfir[rsamp, jpar] = (deriv[rsamp, jpar]^2) * (x$varprior[jpar])
        vartot[rsamp] = vartot[rsamp] + varfir[rsamp, jpar]
        if (jpar == Kpar) {
          for (jjpar in 1:Kpar) delsafirst[rsamp, jjpar] = varfir[rsamp,
                                                                  jjpar]/vartot[rsamp]
        }
      }
    }
    colnames(delsafirst) = colnames(x$X)
    delsaFirst[[timeStep]][[metric]]$delsafirst = delsafirst
    assign(id, x, parent.frame())
  }
  }
}


# the default plots from the Sensitivity Packages
writePlotDir <- paste0(runDir, "/plots")
dir.create(writePlotDir)
obj = x

#plot1 # these plots are only provided for the daily timestep and the objective function as the metric
temp = as.data.frame(delsaFirst$daily$obj$delsafirst)
names(temp) <-  names(x_all)[2:ncol(x_all)]
temp$id <- 1:nrow(temp)
temp = reshape2::melt(temp, id.var = "id")
gg <- ggplot2::ggplot(data = temp, ggplot2::aes(x = value,
                                                 colour = variable)) + 
  ggplot2::stat_ecdf() +
  ggplot2::scale_x_continuous("DELSA results for first order sensitivity") +
  ggplot2::scale_y_continuous("Cum. frequency") +
  ggplot2::labs(title = "CDF of first order sensitivity across parameter space")
ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_CDF_DELSA.png"),
       plot=gg, units="in", width=16, height=8, dpi=300)

#plot2
temp$y <- obj$y[temp$id]

temp2 = as.data.frame(obj$X0)
names(temp2) <-  names(x_all)[2:ncol(x_all)]
temp2$id <- 1:nrow(temp2)
temp2 = reshape2::melt(temp2, id.var = "id")
temp2$x <- temp2$value
temp2$value <- NULL
temp = merge(temp, temp2)

gg <- ggplot2::ggplot(data = temp) + ggplot2::geom_point(ggplot2::aes(x = value,
                                                                       y = y)) + 
  ggplot2::scale_x_continuous(name = "DELSA first order sensitivity") +
  ggplot2::scale_y_continuous(name = "Model output") +
  ggplot2::facet_wrap(~variable, scales = "free") +
  ggplot2::labs(title = "First order sensitivity as related to model response")
ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_DELSA_model_response.png"),
       plot=gg, units="in", width=16, height=8, dpi=300)

#plot3
gg <- ggplot2::ggplot(data = temp) + ggplot2::geom_point(ggplot2::aes(y = value,
                                                                      x = x, colour = y)) + 
  ggplot2::scale_y_continuous(name = "DELSA first order sensitivity") +
  ggplot2::scale_x_continuous(name = "Parameter value") +
  ggplot2::scale_color_continuous(name = "Model response") +
  ggplot2::facet_wrap(~variable, scales = "free") +
  ggplot2::labs(title = "First order sensitivity as as related to parameter value")
ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_DELSA_parameter_value.png"),
       plot=gg, units="in", width=16, height=8, dpi=300)

# Let s do a bootstrap resampling, I want to do this for all the metrics and both temporal resolutions
Quantile <- function(data, indices, SA_quantileFrac = 0.9) {
  d <- data[indices] # allow boot to select sample
  quantileNo <- quantile(d, SA_quantileFrac) #calcualte the quantile
  return(quantileNo)
}

bootRes <- data.table()
for (timeStep in c("hourly", "daily")) {
   for (metric in setdiff(metrics, "fdc")) {
       for (param in 1:(ncol(x_all)-1)) {
            results <- boot(data=delsaFirst[[timeStep]][[metric]]$delsafirst[, param],
		  statistic=Quantile, 
                  R=SA_bootstrap_replicates)
            bootRes <- rbindlist(list(bootRes, data.table(delsaFirst = results$t[,1],
						     timeStep = timeStep, metric = metric, 
						     parameter = names(x_all)[param+1])))
           }
     }
}

# add the plots 
gg <- ggplot(bootRes, aes(parameter, delsaFirst)) + geom_boxplot()+
	facet_grid(metric~timeStep)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_DELSA_uncertainty_estimate.png"),
       plot=gg, units="in", width=16, height=8, dpi=300)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Summary plots

# Plot setup
ggPalette <- gg_color_hue(14)

plotGroups <- list(soil=c('bexp', 'dksat', 'smcmax', 'refkdt', 'slope', 'RETDEPRTFAC', 'LKSATFAC'),
                   other=c('Zmax', 'Expon', 'CWPVT', 'VCMX25', 'MP', 'HVT', 'MFSNO'))

# Hydrographs
gg <- ggplot(data=chrt.d.all) +
  geom_line(aes(x=Date, y=q_cms, color=id, group=id), lwd=0.6) +
  geom_point(aes(x=Date, y=obs, group=id), color = "black") + 
  scale_y_log10() +
  ggtitle(paste0("Model Sensitivity: ", chrt.d.all$site_no[1]))
ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_hydrograph.png"),
       plot=gg, units="in", width=16, height=8, dpi=300)

# gg <- ggplot(data=chrt.d.all[tag %in% plotGroups[["soil"]],], aes(x=POSIXct, y=q_cms, color=tag, group=id)) +
#   geom_line(lwd=0.6) +
#   scale_y_log10() +
#   facet_wrap(~ tag) +
#   annotate(geom='line', x=chrt.d.all[tag=="control",]$POSIXct, y=chrt.d.all[tag=="control",]$q_cms, color='black', lwd=0.3) +
#   ggtitle(paste0("Model Sensitivity: ", chrt.d.all$site_no[1]))
# ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_soilp_sens.png"),
#        plot=gg, units="in", width=16, height=8, dpi=300)
# 
# 
# gg <- ggplot(data=chrt.d.all[tag %in% plotGroups[["other"]],], aes(x=POSIXct, y=q_cms, color=tag, group=id)) +
#   geom_line(lwd=0.6) +
#   scale_y_log10() +
#   facet_wrap(~ tag) +
#   annotate(geom='line', x=chrt.d.all[tag=="control",]$POSIXct, y=chrt.d.all[tag=="control",]$q_cms, color='black', lwd=0.3) +
#   ggtitle(paste0("Model Sensitivity: ", chrt.d.all$site_no[1]))
# ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_otherp_sens.png"),
#        plot=gg, units="in", width=16, height=8, dpi=300)

# Metric plots
# tmp<-x_archive[,c("tag", metrics)]
# tmp<-melt(tmp, id='tag')
# tmpmax <- ddply(tmp, .(tag, variable), function(x) {max(x$value)})
# tmpmin <- ddply(tmp, .(tag, variable), function(x) {min(x$value)})
# tmpall <- merge(tmpmin, tmpmax, by=c("tag", "variable"))
# names(tmpall)[3:4]<-c("min", "max")
# tmpall <- plyr::join(tmpall, data.frame(tag=unique(tmpall$tag), id=seq(1:length(unique(tmpall$tag)))), by="tag")
# 
# gg <- ggplot(data=tmpall, aes(x=tag, fill=tag)) + 
#   geom_rect(aes(x=tag, xmin=id-0.45, xmax=id+0.45, ymin=min, ymax=max, color=tag), alpha=0.8) + 
#   facet_wrap(~variable, scales="free_y") +
#   ggtitle(paste0("Metric Sensitivity: ", chrt.d.all$site_no[1])) +
#   theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
# ggsave(filename=paste0(writePlotDir, "/", chrt.d.all$site_no[1], "_metric_sens.png"),
#        plot=gg, units="in", width=16, height=8, dpi=300)


# Save and exit
if (parallelFlag) stopCluster(cl)
save.image("proj_data_SENS.Rdata")
quit("no")



