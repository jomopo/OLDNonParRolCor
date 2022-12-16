########################################################################
#:: "plot_rolcol_estim_1win" function from the R package NonParRolCor  #
########################################################################
#:: The function "plot_rolcor_estim_1win" plot the output of the       #
#:: function "rolcor_estim_1win" from the R package NonParRolCor.      #
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

 plot_rolcor_estim_1win <- function(inputdata, corcoefs, CRITVAL, 
              widthwin, left_win, righ_win, varX="X", varY="Y", 
              coltsX="black", coltsY="blue", rmltrd=TRUE, Scale=TRUE, 
              HeigWin1=2.05, HeigWin2=2.75, colCOEF="black", CEXLAB=1.15, 
	      CEXAXIS=1.05, LWDtsX=1, LWDtsY=1, LWDcoef=1, 
              colCRITVAL="black", pchCRIVAL=16) 
 { 

 #---------------------------------------------------------------------#
 # "plot_rolcol_estim_1win" uses the outputs of the function 
 # rolcor_estim_1win to plot in ONLY ONE window-length the correlation 
 # coefficients and their statistic significance. 
 #---------------------------------------------------------------------#
 #########################################################################
 # Parameters (from 1 to 11) that MUST be provided: 
 # 1. inputdata: same data that were used in the rolcor_estim_1win function. 
 # 2. corcoefs:  rolwincor_1win's output named "Correlation_coefficients"
 # 3. CRITVAL:   Critical values obtained with the rolcor_estim_1win function
 #               to define the statistical significance. 
 # 4. widthwin:  rolwincor_1win's output named "widthwin", which indicates
 # 	     	 the window's size to compute the rolling window correlations. 
 # 5. left_win:  rolwincor_1win's output named "left_win", which is used 
 #     	  	 in time domain when the cor. coef are plotted. 
 # 6. righ_win:  rolwincor_1win's output named "righ_win", which is used   
 #     	  	 in time domain when the cor. coef are plotted. 
 # 7. varX:    name of the first ("independent") variable, e.g., "X". 
 # 8. varY:    name of the second ("dependent") variable , e.g., "Y". 
 # 9. coltsX:  the color to be used to plot the first variable, by the 
 #             default the color is "black". 	
 # 10. coltsY: the color to be used to plot the second variable, by 
#              default the color is "blue".  
 #########################################################################
 # 11. rmltrd:  remove (by default is TRUE; please use "FALSE" 
 #              otherwise) linear trends in the data analysed. 
 # 12. Scale:   scale (by default is "TRUE"; please use "FALSE" 
 #   	        otherwise) or "normalize" or "standardize" the data 
 #	        analysed. 
 # 13. HeigWin1: window proportion window to plot the variables (>R?layout).
 # 14. HeigWin2: window proportion to plot the correlation coefficients. 
 # 15. colCOEF: the color to be used to plot the correlation coefficients, 
 #		by default the color is "black". 
 # 16. CEXLAB:  the size (cex) for the labels, the default value is 1.15. 
 # 17. CEXAXIS: the size (cex.axis) for the axis, the default value is 1.05. 
 # 18. LWDtsX:  the line width/s (lwd) for the first variable, it has a 
 #              default value of 1. 
 # 19. LWDtsY:  the line width (lwd) for second variable, by default is 1. 
 # 20. LWDcoef: the line width (lwd) to be used to plot the correlation 
 #    	        coefficients, it has a default value of 1. 
 # 21. colCRITVAL: color used to plot the correlation coefficients that 
 # 		   are statistically significant.  
 # 22: pchCRIVAL:  pch used to plot the correlation coefficients that are 
 # 		   statistically significant.  
 #---------------------------------------------------------------------#

 #---------------------------------------------------------------------#
  # Check 1: removing linear trend - if rmltrd=TRUE
  if(isTRUE(rmltrd)) {
   dtrd_tmp  <- pracma::detrend(inputdata[,-1])
   inputdata <- cbind(inputdata[,1], dtrd_tmp) 
  } 

  # Check 2: scaling data: [X_i - mean(X_i)]/sd(X-i), mean=0 & sd=1
  if(isTRUE(Scale)) { 
   inputdata <- cbind(inputdata[,1], apply(inputdata[,-1], 2, scale))
  }

 # Auxiliary code to be used in the plots! 
  # Transforming input data to time series objects! 
  #Deltat <- diff(inputdata[,1]) 
  NL     <- dim(inputdata)[1]
  NP     <- dim(inputdata)[2]
  freq   <- length((inputdata[,1]))/length(unique(inputdata[,1]))
  ts1    <- ts(inputdata[,2], start=c(inputdata[1,1],1), 
            end=c(inputdata[NL,1],freq), deltat=1/freq) # ts1 = X
  ts2    <- ts(inputdata[,3:NP], start=c(inputdata[1,1],1), 
	    end=c(inputdata[NL,1],freq), deltat=1/freq) #ts2 = X 
  Z      <- cbind(ts1, ts2) # To be used in the running cor./fun. (rollapply)
  #
  time.runcor <- time(ts1)
  Nrun        <- length(time.runcor) 
  #
  rc.ts1_ts2         <- corcoefs[,2]

  # Setting up the graphical parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar=c(3.90, 4.05, 2.2, 3.3) + 0.1) 
  layout(matrix(c(1,2), 2, 1, byrow=FALSE), heights=c(HeigWin1, HeigWin2))
 #---------------------------------------------------------------------#

 #---------------------------------------------------------------------#
 ####################################################################
 # Setting up the graphical parameters
 ####################################################################
 #### pending to use open par and close, it's very important for CRAN
 ####################################################################
 # Plot data 
 plot(ts1, t="l", col=coltsX, las=1, xlab="", ylab="", xaxt="n", 
  yaxt="n", xaxs="i", yaxs="i", xlim=c(time(ts1)[1], time(ts1)[NL]), 
  main=paste(varX, " vs. ", varY, sep=""), lwd=LWDtsX)
 points(ts2, t="l", col=coltsY, las=1, xlab="", ylab="", lwd=LWDtsY)
 axis(1, at=pretty(time.runcor), labels=pretty(time.runcor), 
  cex.axis=CEXAXIS)
 axis(2, at=pretty(ts1), labels=pretty(ts1), col=coltsX, las=1, 
  cex.axis=CEXAXIS)
 axis(4, at=pretty(ts2), labels=pretty(ts2), col=coltsY, las=1, 
  cex.axis=CEXAXIS, cex.lab=CEXLAB)
 mtext(1, text="Time", line=2.5, cex=CEXLAB)
 mtext(2, text=varX, col=coltsX, line=2, cex=CEXLAB, las=3)
 mtext(4, text=varY, col=coltsY, line=2, cex=CEXLAB, las=3)
 # Plot coef. cor. values 
 id.cor <- which(abs(rc.ts1_ts2) > CRITVAL[,2], arr.ind=TRUE)
 plot(time.runcor[left_win:(Nrun-righ_win)], rc.ts1_ts2, t="l", yaxt="n", 
  xaxs="i", yaxs="i", xlim=c(range(time.runcor)[1], range(time.runcor)[2]),  
  ylim=c(-1,1), cex.lab=CEXLAB, cex.axis=CEXAXIS, lwd=LWDcoef,  
  xlab="Time", ylab="Dynamic correlation coefficient", las=1, 
   col="black", main=paste(varX, " vs. ", 
   varY, " (", "window-length= ", widthwin, ")", sep="")) 
 points(time.runcor[left_win:(Nrun-righ_win)][id.cor], rc.ts1_ts2[id.cor], 
  t="p", col=colCRITVAL, pch=pchCRIVAL) 
 axis(2, seq(-1, 1, 0.20), las=1, cex.axis=CEXAXIS)
 abline(h=0, col=1, lty=1)

 
 } # End of the function plot_rolcol_estim_1win 
