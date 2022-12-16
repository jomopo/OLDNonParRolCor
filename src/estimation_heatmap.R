########################################################################
#:: "rolcor_estim_heatmap" function from the R package NonParRolCor    #
########################################################################
#:: The function "rolcor_estim_heatmap" estimates the rolling window   #
#:: correlation for several window-lengths and their statistical       #
#:: significance through a non-parametric computing-intensive method.  #
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
#:: the "rolwincor_heatmap" function".   			       #
########################################################################

########################################################################
#:: Led, designed and programmed by Josué M. Polanco-Martinez,         #
#:: Email: josue.m.polanco@gmail.com                                   #
#:: The parallelized adaptation was carried out by José L.             #
#:: López-Martínez, Email: jose.lopez@correo.uady.mx                   #
########################################################################
#   Copyright (C) 2021/2022 by Josué M. Polanco-Martínez and 	       #
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
#   NonParRolCor is distributed in the hope that it will be 
#   useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
#   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with NonParRolCor If not, see <http://www.gnu.org/licenses/>.
#
########################################################################

 rolcor_estim_heatmap <- function(inputdata, CorMethod="pearson", 
                            typewidthwin="FULL", widthwin_1=3, 
                            widthwin_N=dim(inputdata)[1], Align="center", 
                            rmltrd=TRUE, Scale=TRUE, MCSim=1000, 
	                    prob=0.95, Np=2)  

 { # Open the main function

 #---------------------------------------------------------------------#
 #                Description of arguments (INPUTS):                   #
 #---------------------------------------------------------------------#
 # 1. inputdata:  Matrix of 3 columns: time, first variable  (e.g. "X")
 #                and second variable (e.g. "Y").
 # 2. CorMethod:  Method used to estimate the correlations, by default 
 #                is "pearson" but other two options ("spearman" and
 #                "kendall") can be used (please look at: R>?cor.test).
 # 3. typewidthwin: "FULL" is to estimate the windows from 2, 4, ..., 
 #                to dim(inputdata)[1]) if Align is equal to "left" or 
 #                "right", or from 3, 5,..., to dim(inputdata)[1]) if  
 #                Align is "center". The other option is "PARTIAL",  
 # 		  please you should take into account that widthwin_1 &
 #                widthwin_N MUST be ODD if the option Align is "center".
 # 4. widthwin_1: First value for the size of the windows (by default is  
 # 	 	  3) when the option typewidthwin="PARTIAL" is selected. 
 # 5. widthwin_N: Last value for the size of the windows (by default is 
 #                dim(inputdata)[1], i.e. number of datums in inputdata) 
 # 		  when the option typewidthwin="PARTIAL" is selected.
 # 6. Align:      To align the rolling object, NonParRolCor uses three 	
 #                options: "left, "center", and "right" (please look at: 
 # 		  R>?running). However, there are some restrictions, 
 #                which have been described lines above. We recommend
 # 	          to use the "center" option to ensure that variations in 
 #                the correlations are aligned with the variations in the 
 # 	          relationships of the variables under study, rather than 
 #                being shifted to left/right (Polanco-Martínez 2019/20), 
 # 		  but this imply that the window-lengths MUST be ODD. 
 # 7. rmltrd:     Remove (by default is TRUE; please use "FALSE" 
 #                otherwise) the linear trend in the data analysed. 
 # 8. Scale:      Scale (by default is "TRUE"; please use "FALSE" 
 #   	          otherwise) or "normalize" or "standardize" the data 
 #		  analysed. 
 # 9. MCSim:      Number of Monte-Carlo simulations to permute the second 
 #                time series, by default MCSim=1000.         
 # 10. prob:      Numeric vector of probabilities with values in [0,1], 
 # 	          by default prob=0.95 (p=0.05), please look at 
 #                R?quantile for extra/more information. 	       
 # 11. Np:        Number of cores, by default is 2 (please verify the 
 #   	          number of cores of your computer. It's not advisable 
 #                to use the maximum number of cores). 
 #---------------------------------------------------------------------#

 #---------------------------------------------------------------------- 
 # 			    Inputs "Checks"  
 #---------------------------------------------------------------------- 
 # Check 1: number of columns, inputdata MUST contain three columns: 
 #          time, first variable (e.g. "X") & second variable (e.g. "Y")  
 if (dim(inputdata)[2] != 3) 
  stop("\n The input data MUST be an array or matrix with 3 columns 
   (first column the time and the second and third columns the 
   variables under analysis. Thank you for using our NonParRolCor 
   package. \n")

 # "Check" 2: the time steps MUST be regular/evenly - no gaps! 
 #Deltat <- diff(inputdata[,1])  # Deltat is the temporal resolution! 
 #if (length(unique(Deltat)) != 1)
  cat("\n W A R N I N G: The input data must be regular (evenly spaced 
   time). Otherwise, please, consider to address this drawback before 
   using NonParRolCor. We recommend you our BINCOR package and method 
   (also in CRAN; https://cran.r-project.org/package=BINCOR), but other 
   packages and methods can be used. Thank you so much for using our 
   NonParRolCor package. \n")

 #######################################################################
 # Check 3.1: if Align="center" only estimate windows of the form 2p + 1
 # That is, widthwin_1 and widthwin_N MUST be odd numbers 
 if (typewidthwin == "PARTIAL" & Align == "center") { 
  if (widthwin_1 %% 2 == 0 || widthwin_N %% 2 == 0) { 
   stop("\n widthwin_1 or widthwin_N is/are EVEN number/s and these  
   (both) MUST be ODD. Thank you for using our NonParRolCor package. \n")
  }

 # Check 3.2: initial and final values for the window lengths! 
 if (widthwin_1 == widthwin_N) {
   stop("\n The initial and final value for the window-lengths are 
    the same. Please, modify these values. Thank you for using 
    our NonParRolCor package. \n") 
  }
 }

 #######################################################################
 # Check 4: removing linear trend - if rmltrd=TRUE
 if(isTRUE(rmltrd)) {
  ts1.tmp   <- cbind(inputdata[,1], c(pracma::detrend(inputdata[,2])))
  ts2.tmp   <- cbind(inputdata[,1], c(pracma::detrend(inputdata[,3])))
  inputdata <- cbind(ts1.tmp[,1], ts1.tmp[,2], ts2.tmp[,2])  
 } 

 # Check 5: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
  if(isTRUE(Scale)) { 
  inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
 }
 #----------------------------------------------------------------------

 #---------------------------------------------------------------------- 
 # 		Transforming input data to time series object 
 #---------------------------------------------------------------------- 
 NL   <- dim(inputdata)[1]
 freq <- length((inputdata[,1]))/length(unique(inputdata[,1]))
 ts1  <- ts(inputdata[,2], start=c(inputdata[1,1],1), 
            end=c(inputdata[NL,1],freq), deltat=1/freq)
 ts2  <- ts(inputdata[,3], start=c(inputdata[1,1],1), 
	    end=c(inputdata[NL,1],freq), deltat=1/freq)

 time.runcor <- time(ts1)
 Nrun        <- length(time.runcor)   

 #----------------------------------------------------------------------
 # 	Procedure to estimate the windows and the number of windows 
 # 	to compute the rolling window correlations
 # ------------------------------------------------------------------
 # "typewidthwin" indicates how the rolling windows are computed: for 
 # the "FULL" option the window correlations are computed for all the 
 # window-lengths from 2 or 3 to NL (number of datums in inputdata). 
 # For the "PARTIAL" option the window correlations are computed from  
 # widthwin_1 to widthwin_N. "nwin" is the maximum number of windows.
 # NL is the number of elements of the time series under study.
 #----------------------------------------------------------------------

 #
 if (Align == "left" || Align == "right")  { 
  if (typewidthwin == "FULL") {
   if (NL %%2 == 0) { 
    nwin <- NL/2 - 1 
   } else {
    if (NL %%2 != 0)   { 
     nwin <- floor(NL/2) 
    } 
   } 
   Rwidthwin <- seq(4, NL, 2) # length(Rwidthwin) = nwin 
  }
  if (typewidthwin == "PARTIAL") { 
   nwin      <- length(seq(widthwin_1, widthwin_N, 2)) 
   Rwidthwin <- seq(widthwin_1, widthwin_N, 2) # 
  }
 }  

 #
 if(Align == "center") { 
  if (typewidthwin == "FULL") { 
   if (NL %%2 == 0) { 
    nwin <- NL/2 - 1 
   } else {
     if (NL %%2 != 0)   { 
     nwin <- floor(NL/2) 
     } 
   } 
    Rwidthwin <- seq(3, NL, 2) # length(Rwidthwin) = nwin 
   }
 
  if (typewidthwin == "PARTIAL") { 
   nwin      <- length(seq(widthwin_1, widthwin_N, 2)) 
   Rwidthwin <- seq(widthwin_1, widthwin_N, 2) # length(Rwidthwin) = nwin 
  } 
 }
 #----------------------------------------------------------------------

 #----------------------------------------------------------------------
 # 	Array to save the correlation coefficients and CRITVAL 
 #----------------------------------------------------------------------
 the_matrixCOR  <- array(NA, c(nwin, NL-2)) # cor. coefficients 
 CRITVAL        <- array(NA, c(nwin, 2))    # critical of maximum correlation
 
 #############	        Auxiliary   function              ##############
 # Function to estimate the correlation coefficients 
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

 #############	       START:    The big-for 	          ##############
 for (w in 1:nwin) { 
  cor_pval_run <- gtools::running(ts1, ts2, fun=cor_pval.fun,
                   width=Rwidthwin[w], align=Align)

  ######################################################################
  # Statistically significance test via MC sims: from Telford (2013), 
  # Polanco-Martinez (2020), Polanco-Martinez and Lopez-Martinez (2021). 
  ######################################################################
  res.tmp  <- foreach(k=1:MCSim, .combine='c') %dopar%{
     ts2.tmp <- sample(ts2)
     z.tmp   <- max(abs(gtools::running(ts1, ts2.tmp, fun=cor,
                    width=Rwidthwin[w], align=Align)))
  }
  
  CRITVAL[w,] <- cbind(Rwidthwin[w], quantile(res.tmp, prob=prob))
 ######################################################################

 ######################################################################
 #############	These pieces of code are used to align   ##############
 #############  the correlation coefficients, and the    ##############
 #############  CRITVALs to their corresponding "times"  ############## 
  # if widthwin is even and left or right  
  if(Align == "left" || Align == "right") { 
   left_win <- Rwidthwin[w]/2 
   righ_win <- Rwidthwin[w]/2 
  }
  # if widthwin is odd and center 
  if(Align == "center") { 
   left_win <- (Rwidthwin[w] - 1)/2   
   righ_win <- (Rwidthwin[w] - 1)/2 + 1 
  }
 
  #####################################################################
  the_matrixCOR[w,left_win:(Nrun-righ_win)] <- cor_pval_run 

 }
 #############	END:	The big-for 	#############

 stopCluster(cl) # To stop the cluster (parallelized process) 

 #----------------------------------------------------------------------
 #############        Setting up the outputs              ############## 
 namesLS     <- c("matcor", "CRITVAL", "NoWindows", "Windows", 
                  "left_win", "righ_win", "MCSim", "prob")  
 LIST        <- list(the_matrixCOR, CRITVAL, nwin, Rwidthwin, 
	             left_win, righ_win, MCSim, prob)
 names(LIST) <- namesLS
 return(LIST)  

} # Close the main function 
 
