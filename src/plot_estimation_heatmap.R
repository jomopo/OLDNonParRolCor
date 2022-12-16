########################################################################
#: "plot_rolcor_estim_heatmap" function from the R package NonParRolCor#
########################################################################
#:: The function "plot_rolcor_estim_heatmap" plot the time series under#
#:: study and a heat map of the rolling window correlation coefficients#
#:: and its statistical significance that is estimated through a       #
#:: non-parametric computing-intensive method. Methods inspired/derived#
#:: from R. Telford (2013):                                            #        
#:: <https://quantpalaeo.wordpress.com/2013/01/04//> and from Polanco- #
#:: Martinez & Lopez-Martinez (2021), "A non-parametric method to test # 
#:: the statistical significance in rolling window correlations, and   #
#:: applications to ecological time series", Ecological Informatics 64,# 
#:: 101379, https://doi.org/10.1016/j.ecoinf.2021.101379               #
########################################################################
#:: Some pieces of code come from the R "RolWinMulCor" package by      #
#:: Josue M. Polanco-Martinez (2020, Ecological Informatics 60, 101163)# 
#:: https://doi.org/10.10162/j.ecoinf.2020.101163, specifically, from  # 
#:: the "plot_heatmap.R" function.		                       # 
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

 plot_rolcor_estim_heatmap <- function(inputdata, corcoefs, CRITVAL, 
                               Rwidthwin="", typewidthwin="", widthwin_1=3, 
                               widthwin_N=dim(inputdata)[1], varX="X", 
                               varY="Y", coltsX="black", coltsY="blue", 
                               LWDtsX=1, LWDtsY=1, CEXLAB=1.15, CEXAXIS=1.05) 

 { # Open the main function

 #---------------------------------------------------------------------#
 #                Description of arguments (INPUTS):                   #
 #---------------------------------------------------------------------#
 # 1. inputdata: The same data matrix (time, first and second variable)#
 #               that was used in the "estimation_NonParRolCor" function. 
 # 2. corcoefs:  estimation_NonParRolCor's output named 
 #               "Correlation_coefficients".
 # 3. CRITVAL:   The critical values computed through the function 
 #               "estimation_NonParRolCor".
 # 4. Rwidthwin: estimation_NonParRolCor's output named "widthwin", 
 #               which indicates the window's size to compute the 
 #               rolling window correlations. 
 # 5. typewidthwin: to plot all the possible window-lengths ("FULL") or 
 # 		 only a band of window-lengths ("PARTIAL"). 
 # 6. widthwin_1: First value for the size of the windows (by default is  
 # 	 	 3) when the option typewidthwin="PARTIAL" is selected. 
 # 7. widthwin_N: Last value for the size of the windows (by default is 
 #               dim(inputdata)[1], i.e. number of datums in inputdata) 
 # 		 when the option typewidthwin="PARTIAL" is selected.
 #######################################################################
 # Parameters that are not strictly necessary to define:
 # 8.  varX:    Name of the first variable, e.g. "X". 
 # 9.  varY:    name of the second variable, e.g. "Y". 
 # 10. coltsX:  The color to be used to plot the first (independent) 
 #              variable, by the default the color is "black". 	
 # 11. coltsY:  The color to be used to plot the second (dependent)
 #  	        variable, by default the color is "blue".  
 # 12. LWDtsX:  The line width (lwd) for the first variable, it has  
 #		a default value of 1.  
 # 13. LWDtsY:  The line width (lwd) for the second variable, by default is 1. 
 # 14. CEXLAB:  The size (cex) for the labels, the default value is 1.15. 
 # 15. CEXAXIS: The size (cex.axis) for the axis, the default value is 1.05. 
 #---------------------------------------------------------------------#

 #---------------------------------------------------------------------#
 #               Auxiliary code to be used in the plots 
 #---------------------------------------------------------------------#
 if (typewidthwin == "FULL") { 
   rango     <- range(Rwidthwin)
   NoYaxis   <- floor(length(Rwidthwin)/5)  
 } 
 if (typewidthwin == "PARTIAL") { 
   nwin      <- length(seq(widthwin_1, widthwin_N, 2)) 
   Rwidthwin <- seq(widthwin_1, widthwin_N, 2) # length(Rwidthwin) = nwin 
  if (nwin <= 10) { 
   rango     <- range(Rwidthwin)
   NoYaxis   <- 2 
  }
  if (nwin > 10) { 
   rango     <- range(Rwidthwin)
   NoYaxis   <- floor(length(Rwidthwin)/5) 
  }  
 }   
 nwin <- length(Rwidthwin)

 #---------------------------------------------------------------------#
 #         Transforming the input data to time series objects
 #---------------------------------------------------------------------#
 NL     <- dim(inputdata)[1]
 freq   <- length((inputdata[,1]))/length(unique(inputdata[,1]))
 # Please note that we scale (standardized/normalized) 
 ts1    <- ts(scale(inputdata[,2]), start=c(inputdata[1,1],1), 
              end=c(inputdata[NL,1],freq), deltat=1/freq) 
 ts2    <- ts(scale(inputdata[,3]), start=c(inputdata[1,1],1), 
              end=c(inputdata[NL,1],freq), deltat=1/freq)
 #
 time.runcor <- time(ts1) 
 Nrun        <- length(time.runcor) 
 #
 the_matrixCOR <- corcoefs

 #---------------------------------------------------------------------#
 #                Setting up the graphical parameters
 #---------------------------------------------------------------------#
 oldpar <- par(no.readonly = TRUE)
 on.exit(par(oldpar))
 par(mar=c(2.35, 4.95, 2.2, 3.5) + 0.1) 
 layout(matrix(c(1,2,3), 3, 1, byrow=FALSE), heights=c(2, 3.5, 0.65))
 # Palette!
 Ncol     <- length(Rwidthwin)
 # Number of colors MUST be equal to length(Rwidthwin) or Ncol 
 Palette  <- colorspace::diverge_hcl(4*Ncol, c=c(100,0), l=c(50,90), power=1.3)
 # Colorbar! 
 rangev   <- seq(min(the_matrixCOR, na.rm=TRUE), max(the_matrixCOR, 
                 na.rm=TRUE), length.out=11)
 rangebar <- matrix(rangev, nrow=1, ncol=11, byrow=TRUE)
 # To make transparent a color (e.g. green)
 makeTransparent <- function(someColor, alpha=0.15) scales::alpha(someColor, alpha)

 #---------------------------------------------------------------------#
 # Plotting the time series (data input) and the heat map (correlations)
 #---------------------------------------------------------------------#
 # Plot data (plot 1)
 to_ylim <- range(range(ts1), range(ts2))
 plot(ts1, t="l", col=coltsX, las=1, xlab="", ylab="", xaxt="n", 
  yaxt="n", xaxs="i", yaxs="i", xlim=c(time(ts1)[1], time(ts1)[NL]),  
  ylim=to_ylim, main=paste(varX, " vs. ", varY, sep=""), lwd=LWDtsX)
 points(ts2, t="l", col=coltsY, las=1, xlab="", ylab="", lwd=LWDtsY)
 axis(1, at=seq(time.runcor[1], time.runcor[Nrun], 
  length.out=length(time.runcor)), labels=inputdata[,1])
 axis(2, at=pretty(ts1), labels=pretty(ts1), col=coltsX, las=1, 
  cex.axis=CEXAXIS)
 axis(4, at=pretty(ts2), labels=pretty(ts2), col=coltsY, las=1, 
  cex.axis=CEXAXIS, cex.lab=CEXLAB)
 mtext(1, text="Time", line=2.75, cex=CEXLAB)
 mtext(2, text=varX, col=coltsX, line=2.5, cex=CEXLAB, las=3)
 mtext(4, text=varY, col=coltsY, line=2.5, cex=CEXLAB, las=3)
 # -------------------------------------------------------------
 # Multi time-scale window correlation (plot 2)
 # To take into account the statistical significance in the image-plot! 
 new_the_matrixCOR <- the_matrixCOR
 if(1){
  for (k in 1:nwin) { 
   id.xy <- which(abs(the_matrixCOR[k,]) <= CRITVAL[k,2], arr.ind=TRUE) 
   new_the_matrixCOR[k,id.xy] <- rep(NA, length(id.xy)) 
  } 
 }
 image(t(new_the_matrixCOR), xlab="", ylab="", las=1, 
  col=Palette, xaxt="n", yaxt="n")
 contour(t(the_matrixCOR), add=TRUE, drawlabels=TRUE)
 axis(1, at=(seq(0,1,length.out=length(time.runcor))), 
  labels=inputdata[,1])
 # Set up Y axis 
 # To verify if NoYaxis is odd or even number! 
 if(NoYaxis %% 2 == 0) { 
 # Case 1: 
 NoYaxis_CASE1 <- NoYaxis  
 to_at_CASE1   <- seq(rango[1], rango[2], by=NoYaxis_CASE1) 
 length_CASE1  <- length(to_at_CASE1)
 } else {
 # Case 2: 
 NoYaxis_CASE2 <- NoYaxis - 1 
 to_at_CASE2   <- seq(rango[1], rango[2], by=NoYaxis_CASE2) 
 length_CASE2  <- length(to_at_CASE2)
 }
 # 
 if (NoYaxis %% 2 == 0) { 
 axis(2, at=seq(0, 1, length.out=length_CASE1), 
  labels=to_at_CASE1, cex.axis=0.95*CEXAXIS, las=1)
 } else { 
   axis(2, at=seq(0, 1, length.out=length_CASE2), 
    labels=to_at_CASE2, cex.axis=0.95*CEXAXIS, las=1)
  }
 mtext(1, text="Time", line=2.75, cex=CEXLAB)
 mtext(2, text="Window-length (time-scales)", line=3.35, cex=CEXLAB)
 # -------------------------------------------------------------
 # Colorbar (plot 3) 
 image(z=t(rangebar), axes=FALSE, col=Palette, frame.plot=TRUE,
  yaxt="n", xaxt="n") 
 axis(1, at=round(seq(0,1,length.out=11),2), labels=round(rangebar, 
  digits=2), cex.axis=CEXAXIS, las=1)

} # End of the main function 
