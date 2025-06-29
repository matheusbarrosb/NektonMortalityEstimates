
iLCCC <- function(mids = NULL,
                  catch = NULL,
                  K = NULL,
                  Linf = NULL,
                  t0 = NULL,
                  binsize = NULL,
                  ex.points = 1,
                  plot = FALSE,
                  plot.yup.lim = NULL,
                  plot.ylow.lim = NULL,
                  plot.xlow.lim = NULL,
                  plot.xup.lim = NULL,
                  N_runs = 1000) {
  
  if (is.matrix(catch)) {
    stop(noquote("Catch must be arranged as a vector"))
  }
  
  if (length(mids) != length(catch)) {
    stop(noquote("Catch and midlengths must be vectors of same length"))
  }
  
  mids <- as.vector(mids)
  catch <- as.vector(catch)
  #upper and lower limits of each size class
  class.min <- mids - (binsize/2)
  class.max <- mids + (binsize/2)
  
  
  if (t0 <- FALSE) {
    t0 <- exp(-0.3922 - 0.2752*log(Linf) - 1.038*log(K))
  } else {
    t0 <- t0
  }
  
  rel.age <- ifelse(mids <= Linf, t0 - (1/K)*log(1 - mids/Linf), NA) # age at any given size class
  dt <- ifelse(mids <= Linf, (t0 - (1/K)*log(1 - class.max/Linf)) - (t0 - (1/K)*log(1 - class.min/Linf)), NA) # amount of time from one size class to another
  # if a midlength is larger than the asympotic size, it is simply excluded from subsequent analysis
  # as calculating dt with midlenths > Linf leads to NaN values
  
  lnNdt <- ifelse(is.na(dt), NA, log((catch+1)/dt)) # only calculate logged catch if the divisor is != NA
  
  df <- data.frame(lnNdt, rel.age)
  df <- na.exclude(df)
  df[df == 0] <- NA #exclude zero values
  
  yvar <- as.numeric(df$lnNdt)
  xvar <- df$rel.age
  
  selection <- c(which(yvar == max(yvar)) + 1,
                 which(xvar == max(xvar)) - ex.points) # select the data rows after the mode (highest point) to fit the line
  # if you suspect of undersampling, modify to c(which(yvar == max(yvar)) + 1,
  # which(xvar == max(xvar))-1) or -2 or more, this should exclude the largest size classes in the sample from the regression
  
  df.cc <- as.data.frame(cbind(xvar,yvar))
  df.selec.cc <- df.cc[selection[1]:selection[2],]# creates a data frame with only the selected rows
  
  df.selec.cc[is.finite(rowSums(df.selec.cc)),] # select only finite values
  df.selec.cc[is.numeric(rowSums(df.selec.cc)),] # exclude NaN measurements
  
  # some Infs and NaNs can pop up during the process when logarithms end up as imaginary or complex numbers
  
  df.selec.cc <- subset(df.selec.cc, df.selec.cc$yvar > 1)
  
  lm.cc <- lm(yvar ~ xvar,
              data = df.selec.cc)
  
  reg.output <- summary(lm.cc)
  intercept.reg <- reg.output$coefficients[1]
  slope.reg <- -reg.output$coefficients[2]
  se.slope.reg <- reg.output$coefficients[4]
  
  # bootstrap to obtain C.I.s for the total mortality
  
  boot.Z <- NULL # empty object to store the Z replicates
  
  for (i in 1:N_runs) {
    
    resamp.df = df.selec.cc[sample(1:nrow(df.selec.cc),
                                   nrow(df.selec.cc),
                                   replace = TRUE),]
    
    resamp.df <- resamp.df[apply(resamp.df, 1,
                                 function(row) all(row != 0)),]
    resamp.df[is.na(resamp.df) | resamp.df=="Inf"] <- NA
    resamp.df[is.na(resamp.df) | resamp.df=="NaN"] <- NA
    
    lm.boot <- lm(yvar ~ xvar,
                  data = resamp.df,
                  na.action = na.omit)
    
    boot.Z <- c(boot.Z, lm.boot$coefficients[2])
    
  }
  
  Z <- median(-boot.Z)
  Z.range <- range(-boot.Z)
  Z.CI <- quantile(-boot.Z, probs = c(0.05, 0.95), na.rm = TRUE)
  boot.Z <- as.vector(boot.Z)
  
  output <- list(intercept.reg,
                 Z, se.slope.reg,
                 Z.CI, Z.range,
                 as.vector(-boot.Z),
                 df.selec.cc)
  names(output) <- c("Intercept (a)", "Z",
                     "Standard Error for Z",
                     "Z bootstrapped 90% C.I.",
                     "Range",
                     "Z posterior distribution",
                     "Regression input")
  
  preds <- predict(lm.cc, pred.df = data.frame(x = df.selec.cc$xvar),
                   interval = "confidence")
  
  par(mfrow = c(1,2))
  
  if (plot == TRUE) {
    
    plot(yvar ~ xvar, xlab = "Relative age (years)",
         ylab = "ln(Catch)/dt",
         main = paste("Length-converted catch curve"),
         cex.main = 0.8, ylim = c(plot.ylow.lim, plot.yup.lim),
         xlim = c(plot.xlow.lim, plot.xup.lim))
    mtext(paste("Median Z =", round(Z,3),
                ", 95% C.I. =", round(Z.CI[1],3), "to", round(Z.CI[2],3)),
          side = 3, cex = 0.7)
    points(df.selec.cc$yvar ~ df.selec.cc$xvar, col = "black",
           pch = 20)
    abline(lm.cc)
    lines(df.selec.cc$xvar, preds[,3], lty = "dashed")
    lines(df.selec.cc$xvar, preds[,2], lty = "dashed")
    
    
    hist(-boot.Z, main = "Posterior distribution of Z",
         xlab = "Z", breaks = 40, col = 0, cex.main = 0.8)
    
  }
  
  return(output)
  
}
# 
# #test data
# md <- seq(10,80,4) #midlengths or size classes
# ct <- c(1,3,6,8,14,20,28,33,37,46,36,24,16,12,8,6,3,1) #catch per size class
# CC1 <- iLCCC(mids = md, catch = ct, K = 0.4, Linf = 90, t0 = FALSE, binsize = 4, plot = TRUE)

