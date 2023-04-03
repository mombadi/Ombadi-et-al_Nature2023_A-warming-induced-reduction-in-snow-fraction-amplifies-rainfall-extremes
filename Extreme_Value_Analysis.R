#================================================
# Reproducibility script: GEV analysis
# Ombadi et al., "A warming-induced reduction in
#     snow fraction amplifies rainfall extremes"
# Mark Risser
# mdrisser@lbl.gov
# Lawrence Berkeley National Laboratory
# March, 2023
#================================================

# Libraries
library(climextRemes)
library(raveio)
library(R.matlab)
library(stringr)

# Set working directory to folder with .mat files
# setwd(<local directory>)

# Files to loop over
pr_directories <- list.files()
# > pr_directories
# [1] "AWI_hist_ams_lp.mat"        "AWI_ssp585_ams_lp.mat"      "BCC_hist_ams_lp.mat"       
# [4] "BCC_ssp585_ams_lp.mat"      "CMCC_hist_ams_lp.mat"       "CMCC_ssp585_ams_lp.mat"    
# [7] "ECEarth3_hist_ams_lp.mat"   "ECEarth3_ssp585_ams_lp.mat" "GFDL_hist_ams_lp.mat"      
# [10] "GFDL_ssp585_ams_lp.mat"     "MPI_hist_ams_lp.mat"        "MPI_ssp585_ams_lp.mat"     
# [13] "MRI_hist_ams_lp.mat"        "MRI_ssp585_ams_lp.mat"      "TAI_hist_ams_lp.mat"       
# [16] "TAI_ssp585_ams_lp.mat" 

# Extract information
files_df <- data.frame(
  Model = sapply(str_split(pr_directories, "_"), function(x){x[[1]]}),
  Expt = sapply(str_split(pr_directories, "_"), function(x){x[[2]]})
)

# Wrapper to calculate return values from hist for fixed duration
hist_rvs <- function(ams_d){
  Nlon <- dim(ams_d)[2]
  Nlat <- dim(ams_d)[1]
  gev_array <- array(NA, dim=c(Nlat,Nlon,3))
  rv_array <- array(NA, dim=c(Nlat,Nlon,4))
  for(i in 1:Nlon){
    if(i %% 10 == 0) cat(i, " ")
    for(j in 1:Nlat){
      for(d in 1){
        if(sum(is.na(ams_d[j,i,])) == 0){
          possibleError <- tryCatch(
            test.fit <- fit_gev( y = ams_d[j,i,], 
                                 getParams = TRUE,
                                 returnPeriod = c(2,5,10,20),
                                 optimArgs = list( control = list(maxit=2000) ) ),
            error=function(e) e
          )
          if(inherits(possibleError, "error")){
            next
          } 
          
          if( test.fit$info$failure == FALSE ){ # For successful optimization
            rv_array[j,i,] <- test.fit$returnValue
            gev_array[j,i,] <- test.fit$mle
          } else{
            # Try initialization
            possibleError2 <- tryCatch(
              test.fit2 <- fit_gev( y = ams_d[j,i,], 
                                    getParams = TRUE,
                                    returnPeriod = c(2,5,10,20),
                                    initial = list(location = mean(ams_d[j,i,]),
                                                   scale = sd(ams_d[j,i,]),
                                                   shape = 0.1),
                                    optimArgs = list( control = list(maxit=2000) ) ),
              error=function(e) e
            )
            if(inherits(possibleError2, "error")){
              next
            } 
            if( test.fit2$info$failure == FALSE ){ # For successful optimization
              rv_array[j,i,] <- test.fit2$returnValue
              gev_array[j,i,] <- test.fit2$mle
            }
          }
        }
      }
    }
  }
  return(list(gev_array,rv_array))
  
}
# Wrapper to calculate return values from hist for fixed duration
f_ssp_rvs_rps <- function(ams_d, hist_rv_array){
  Nlon <- dim(ams_d)[2]
  Nlat <- dim(ams_d)[1]
  gev_array <- array(NA, dim=c(Nlat,Nlon,3))
  rv_array <- array(NA, dim=c(Nlat,Nlon,4))
  rp_array <- array(NA, dim=c(Nlat,Nlon,4))
  for(i in 1:Nlon){
    if(i %% 10 == 0) cat(i, " ")
    for(j in 1:Nlat){
      for(d in 1){
        if(sum(is.na(ams_d[j,i,])) == 0){
          possibleError <- tryCatch(
            test.fit <- fit_gev( y = ams_d[j,i,], 
                                 getParams = TRUE,
                                 returnPeriod = c(2,5,10,20),
                                 returnValue = hist_rv_array[j,i,],
                                 optimArgs = list( control = list(maxit=2000) ) ),
            error=function(e) e
          )
          if(inherits(possibleError, "error")){
            next
          } 
          
          if( test.fit$info$failure == FALSE ){ # For successful optimization
            gev_array[j,i,] <- test.fit$mle
            rv_array[j,i,] <- test.fit$returnValue
            rp_array[j,i,] <- exp(test.fit$logReturnPeriod)
          } else{
            # Try initialization
            possibleError2 <- tryCatch(
              test.fit2 <- fit_gev( y = ams_d[j,i,], 
                                    getParams = TRUE,
                                    returnPeriod = c(2,5,10,20),
                                    returnValue = hist_rv_array[j,i,],
                                    initial = list(location = mean(ams_d[j,i,]),
                                                   scale = sd(ams_d[j,i,]),
                                                   shape = 0.1),
                                    optimArgs = list( control = list(maxit=2000) ) ),
              error=function(e) e
            )
            if(inherits(possibleError2, "error")){
              next
            } 
            if( test.fit2$info$failure == FALSE ){ # For successful optimization
              gev_array[j,i,] <- test.fit2$mle
              rv_array[j,i,] <- test.fit2$returnValue
              rp_array[j,i,] <- exp(test.fit2$logReturnPeriod)
            }
          }
        }
      }
    }
  }
  return(list(rv_array,rp_array,gev_array))
  
}

#================================================
# Warming: fro Extended Data Table 2
#================================================
# AWI-CM-1-1-MR & 2071 - 2100 & 1950 - 1979 & 4.78  \\
# BCC-CSM2-MR  & 2071 - 2100 & 1950 - 1979 & 4.36  \\
# CMCC-CM2-SR5  & 2071 - 2100 & 1950 - 1979 & 5.63 \\
# EC-Earth3 & 2071 - 2100 & 1950 - 1979 & 5.88   \\
# GFDL-ESM4 & 2071 - 2100 & 1950 - 1979 & 3.86   \\
# MPI-ESM1-2-HR & 2071 - 2100 & 1950 - 1979 & 3.86 \\
# MRI-AGCM-3-2-H & 2071 - 2100 & 1950 - 1979 & 4.94  \\
# TaiESM1 & 2071 - 2100 & 1950 - 1979 & 6.64 \\
warming_df <- data.frame(
  Model = (unique(files_df$Model)),
  Warming = c(4.78, 4.36, 5.63, 5.88, 3.86, 3.86, 4.94, 6.64)
)

#================================================
# Loop over models / durations
#================================================
r_vec <- c(2,5,10,20) # return times
d_vec <- c(3,12,24)   # durations
for(m in (unique(files_df$Model))){
  if(m %in% c("GFDL","MRI")){ # Daily only
    print(m)
    hist_ams_data <- read_mat(pr_directories[files_df$Model == m & files_df$Expt == "hist"])[[1]]
    ssp_ams_data <- read_mat(pr_directories[files_df$Model == m & files_df$Expt == "ssp585"])[[1]]
    
    # baseline fits
    baseline_rv <- hist_rvs(ams_d = hist_ams_data)
    # CC-only fits
    cc_rvs_rps <- f_ssp_rvs_rps(ams_d = hist_ams_data*1.07, hist_rv_array = baseline_rv)
    # ssp fits
    ssp_rvs_rps <- f_ssp_rvs_rps(ams_d = ssp_ams_data, hist_rv_array = baseline_rv[[2]])
    
    # Write to .mat files
    for(r in 1:4){
      mat_obj <- abind::abind(baseline_rv[[2]][,,r], ssp_rvs_rps[[1]][,,r], ssp_rvs_rps[[2]][,,r], cc_rvs_rps[[2]][,,r], along = 3)
      writeMat(con = paste0(m, "_24hour_", r_vec[r], "year.mat"), x = mat_obj)
    }
    
  } else{ # 3hr, 12hr, daily
    for(d in 1:3){
      print(m)
      print(d)
      hist_ams_data <- read_mat(pr_directories[files_df$Model == m & files_df$Expt == "hist"])[[1]]
      ssp_ams_data <- read_mat(pr_directories[files_df$Model == m & files_df$Expt == "ssp585"])[[1]]
      
      # baseline fits
      baseline_rv <- hist_rvs(hist_ams_data[,,d,])
      # CC-only fits
      cc_rvs_rps <- f_ssp_rvs_rps(hist_ams_data[,,d,]*1.07, baseline_rv)
      # ssp fits
      ssp_rvs_rps <- f_ssp_rvs_rps(ssp_ams_data[,,d,], baseline_rv[[2]])
      
      # Write to .mat files
      for(d in 1:3){
        for(r in 1:4){
          mat_obj <- abind::abind(baseline_rv[[2]][,,r], ssp_rvs_rps[[1]][,,r], ssp_rvs_rps[[2]][,,r], cc_rvs_rps[[2]][,,r], along = 3)
          writeMat(con = paste0(m, "_", d_vec[d], "hour_", r_vec[r], "year.mat"), x = mat_obj)
        }
      }
    }
  }
}

