########################################################################
#:: "rolcor_estim_1win" function from the R package NonParRolCor       #
########################################################################
#:: The function "rolcor_estim_1win" estimates the rolling window      #
#:: correlation for only one window and their statistical significance #
#:: estimated through a non-parametric computing-intensive method.     #
#:: Methods inspired/derived from R. Telford (2013):                   #
#:: <https://quantpalaeo.wordpress.com/2013/01/04/> and from Polanco-  #  
#:: Martinez & Lopez-Martinez (2021), "A non-parametric method to test # 
#:: the statistical significance in rolling window correlations, and   #
#:: applications to ecological time series", Ecological Informatics 64,# 
#:: 101379, https://doi.org/10.1016/j.ecoinf.2021.101379               #
########################################################################
#:: Some pieces of code come from the R "RolWinMulCor" package by      #
#:: Josue M. Polanco-Martinez (2020, Ecological Informatics 60, 101163)# 
#:: https://doi.org/10.10162/j.ecoinf.2020.101163, specifically, from  # 
#:: the function "rolwincor_1win."				       # 
########################################################################

########################################################################
#:: Led, designed and programmed by Josué M. Polanco-Martinez,         #
#:: Email: josue.m.polanco@gmail.com                                   #
#:: The parallelized adaptation was carried out by José L.             #
#:: López-Martínez, Email: jose.lopez@correo.uady.mx                   #
########################################################################
#   Copyright (C) 2022 by Josué M. Polanco-Martínez and 	       #
#   José L. López-Martínez            	                               #
#   This file/code is part of the R package NonParRolCor               #
########################################################################
#
#   NonParRolCor is free software: 
#   you can redistribute it and/or modify it under the terms of the GNU 
#   General Public License as published by the Free Software 
#   Foundation, either version 3 of the License, or (at your option) 
#   any later version.
#
#   NonParRolCor is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with NonParRolCor. If not, see <http://www.gnu.org/licenses/>.
#
#####################################################################

 rolcor_estim_1win <- function(inputdata, CorMethod="pearson", widthwin=3, 
                       Align="center", rmltrd=TRUE, Scale=TRUE, MCSim=1000, 
                       Np=2, prob=0.95)  { 

 #---------------------------------------------------------------------#
 # "rolcor_estim_1win" function estimates the rolling (running) window 
 # correlation (coefficients and their respective statical significance) 
 # between TWO time series (bi-variate case) sampled on identical time 
 # points for ONLY ONE window-length.
 #---------------------------------------------------------------------#
 # Description of variables (INPUTS): 
 # 1.  inputdata:  matrix of 3 columns: time, variable 1 (e.g. "X"), 
 # 		   variable 2 (e.g. "Y"). Please verify that the inputdata
 #                 is a matrix object. 
 # 2.  CorMethod:  method used to estimate the correlations, by default 
 #                 is "pearson" but other options ("spearman" and "kendall")
 #                 are possible (please see at: R>?cor.test). 
 # 3. widthwin:    window's size to compute the rolling window correlations, 
 #   		   this MUST be >= 3. 
 # 4. Align:       to align the rolling object, NonParRolCor uses 
 #                 the options: "left", "center" and "right." Please 
 #                 note that if "widthwin" is even it's not possible to 
 #   	    	   to use the option "center" (please look at: R>?running).
 # 5. rmltrd:      remove (by default is "TRUE"; please use "FALSE"
 #                 otherwise) linear trend in the data analysed.   
 # 6. Scale:       scale (by default is "TRUE"; please use "FALSE" 
 # 		   otherwise): "normalized" or "standardized" the data 
 # 		   analysed. 
 # 7. MCSim:       Number of Monte-Carlo simulations to permute the second 
 #                 time series, by default MCSim=1000.  
 # 8. Np:         Number of cores, by default is 2 (please verify the 
 #		   number of cores of your computer. It's not advisable 
 #                 to use the maximum number of cores). 
 # 9. prob:       Numeric vector of probabilities with values in [0,1],
 # 		   by default prob=0.95 (p=0.05), please look at 
 # 		   R?quantile for extra/more information. 	
 #------------------------------------------------------------------------#

 #------------------------------------------------------------------------#
  # Check 1: number of columns  - MUST be three columns: 
  #          time, variable X and variable Y 
  if (dim(inputdata)[2] != 3) 
   stop("\n The input data MUST be an array/matrix with 3 columns (first  
    column the time and the second and third the variables under analysis.  
    Thank you for using our NonParRolCor package. \n ")

  # Check 2: the times should be regular/evenly spaced - no gaps! 
  #Deltat <- diff(inputdata[,1])  # Deltat is the temporal resolution! 
  #if (length(unique(Deltat)) != 1)
   cat("\n W A R N I N G: The input data must be regular (evenly spaced 
   time). Otherwise, please, consider to address this drawback before 
   using NonParRolCor. We recommend our BINCOR package and method 
   (also in CRAN; https://cran.r-project.org/package=BINCOR), but other 
   packages and methods can be used. Thank you so much for using our 
   NonParRolCor package. \n")

############################################################################
  # Check 3: widthwin MUST be odd if Align="center" otherwise 
  # gtools::running will not work! 
  if (widthwin %% 2 == 0 & Align == "center" | widthwin < 3) { 
   stop(paste("\n 'widthwin' is EVEN and 'Align' has been defined as 
    `center' or 'widthwin' is < 3.  Thank you for using NonParRolCor
     package. \n"))
  }
############################################################################

  # Check 4: removing linear trend - if rmltrd=TRUE
  if(isTRUE(rmltrd)) {
   ts1.tmp   <- cbind(inputdata[,1], c(pracma::detrend(inputdata[,2])))
   ts2.tmp   <- cbind(inputdata[,2], c(pracma::detrend(inputdata[,3])))
   inputdata <- cbind(ts1.tmp[,1], ts1.tmp[,2], ts2.tmp[,2])  
  } 
 
  # Check 5: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
  if(isTRUE(Scale)) { 
   inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
  }
  # ----------------------------------------------------------------------- 

  # Transforming input data to time series objects! 
  NL   <- dim(inputdata)[1]
  freq <- length((inputdata[,1]))/length(unique(inputdata[,1]))
  ts1  <- ts(inputdata[,2], start=c(inputdata[1,1],1), 
             end=c(inputdata[NL,1],freq), deltat=1/freq) 
  ts2  <- ts(inputdata[,3], start=c(inputdata[1,1],1), 
             end=c(inputdata[NL,1],freq), deltat=1/freq)

  # ------------------------------------------------------------------
  # Computing the rolling window (running) correlation and CRITVAL  

  # Array to save the correlation coefficients and CRITVAL
  the_matrixCOR  <- array(NA, c(1, NL-2)) # cor. coefficients 
  CRITVAL        <- array(NA, c(1, 2))    # critical of maximum correlation
  
  # Auxiliary   function 
  cor_pval.fun <- function(ts1, ts2){
   res_estim   <- cor.test(ts1, ts2, method=CorMethod)
                         c(correlation=res_estim$estimate) 
  }

  ############# To declare the number of cores to be used  ##############  
  ############# to parallelize the computations     	  ##############
  cl <- makeCluster(Np)  # number of cores 
  registerDoParallel(cl) # to register the parallel backend as 
  #registerDoSNOW(cl)    # to register the SNOW parallel backend, i.e.  
      		         # the cluster object "cl" 
 
  cor_pval_run <- gtools::running(ts1, ts2, fun=cor_pval.fun,
                    width=widthwin, align=Align)

  ######################################################################
  # Statistically significance test via MC sims: from Telford (2013), 
  # Polanco-Martinez (2020), Polanco-Martinez and Lopez-Martinez (2021). 
  ######################################################################
  res.tmp  <- foreach(k=1:MCSim, .combine='c') %dopar%{
     ts2.tmp <- sample(ts2)
     z.tmp   <- max(abs(gtools::running(ts1, ts2.tmp, fun=cor_pval.fun,
                    width=widthwin, align=Align)))
  }
  
  CRITVAL[1,] <- cbind(widthwin, quantile(res.tmp, prob=prob))

  stopCluster(cl) # To stop the cluster (parallelized process) 
  
 ######################################################################

  # ESto no tiene nada que ver con la estadistica, mover de lugar! 
  time.runcor  <- time(ts1)
  Nrun         <- length(time.runcor) 
 
######################################################################
  # if widthwin is even  
  if(widthwin %% 2 == 0 &  Align == "left") { 
   left_win <- widthwin/2 
   righ_win <- widthwin/2 
  }
  if(widthwin %% 2 == 0 &  Align == "right") { 
   left_win <- widthwin/2 + 1  
   righ_win <- widthwin/2 - 1  
  }
  # if widthwin is odd 
  if(widthwin %% 2 != 0 &  Align == "left") { 
   left_win <- floor(widthwin/2)
   righ_win <- ceiling(widthwin/2)
  }
  if(widthwin %% 2 != 0 &  Align == "center") { 
   left_win <- (widthwin - 1)/2 + 1  
   righ_win <- (widthwin - 1)/2 
  }
  if(widthwin %% 2 != 0 &  Align == "right") { 
   left_win <- ceiling(widthwin/2) + 1 
   righ_win <- floor(widthwin/2) - 1 
  }
######################################################################

 # Numerical output 
 namesLS <- c("Correlation_coefficients", "CRITVAL",
	      "CorMethod", "widthwin", "left_win", "righ_win") 
 LIST <- list(cbind(time.runcor[left_win:(Nrun-righ_win)], cor_pval_run), 
              CRITVAL, CorMethod, widthwin, left_win, righ_win)
 names(LIST) <- namesLS

 return(LIST)
 
} # End function 


