# model5P <- function(x, y, sd=NULL, W.coef=0.25, Plot=FALSE, Title="", Xlab="", Ylab="", pCol = "indianred2", lCol = "royalblue3", Sub = "Weighted 5P logistic regr.", ...){
source('~/FredScripts/lm.rob.R')
.PL5 <- function(bottom, top, xmid, scal, s,  X){
  yfit <-bottom+(top-bottom)/(1+10^((xmid-X)*scal))^s
  return(yfit)
}
.scePL5 <- function(param, X, yobs, Weights, wcoef) {
  bottom <- param[1]
  top <- param[2]
  xmid <- param[3]
  scal <- param[4]
  s <- param[5]
  ytheo <- .PL5(bottom, top, xmid, scal, s, X)
  residus <- yobs - ytheo
  Weights <- (1/(residus^2))^(wcoef)
  return(sum(Weights*(yobs - ytheo)^2))
}
.scePL4 <- function(param, X, yobs, Weights, wcoef) {
  bottom <- param[1]
  top <- param[2]
  xmid <- param[3]
  scal <- param[4]
  s <- 1
  ytheo <- .PL5(bottom, top, xmid, scal, s, X)
  residus <- yobs - ytheo
  Weights <- (1/(residus^2))^(wcoef)
  return(sum(Weights*(yobs - ytheo)^2))
}
.initPars <- function(x, y){
  bottom = min(y)
  top = max(y)
  xmid = (max(x)+min(x))/2
  z <- (y - bottom)/(top - bottom)
  z[z==0] <- 0.01; z[z==1] <- 0.99
  scal = coef(lm(x ~ log(z/(1-z))))[2]
  scal <- as.numeric(scal)
  s = 1
  return(c(bottom, top, xmid, scal, s))
}
.getPars <- function(model){
  bottom <- model$estimate[1]
  top <- model$estimate[2]
  xmid<-model$estimate[3]
  scal <- model$estimate[4]
  s <- model$estimate[5]
  return(cbind.data.frame(bottom=bottom, top=top, xmid=xmid, scal=scal, s=s))
}
.getPerf <- function(y, yfit){
  test <- summary(lm(yfit ~ y))
  fstat <- test$fstatistic
  p <- pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
  goodness <- test$adj.r.squared
  err <- sqrt(1/(length(yfit)-2)*sum((yfit-y)^2))
  return(cbind.data.frame(goodness=goodness, err=err, p=p))
}

model5P <- function(x, y, sd=NULL, wcoef=0.25, method=c("p4", "p5"), #nfold=5,
                    Plot=TRUE, addLine = TRUE, addPoints = TRUE, addXinf = FALSE,
                    pCol = "aquamarine1", lCol = "red3", Sub = "Weighted 5P logistic regr.", ...){
  
  if(any(is.na(x) | is.na(y))){
    NAs <- union(which(is.na(x)), which(is.na(y)))
    x <- x[-NAs]
    y <- y[-NAs]
  }
  
  method <- match.arg(method)
  switch(method,
         p4 = {.sce <- .scePL4},
         p5 = {.sce <- .scePL5}
  )
  
  init <- .initPars(x, y)
  weights <- rep(1, length(y))
  best <- nlm(f = .sce, p = init, X = x, yobs = y, Weights=weights, wcoef=wcoef)
  
  # Best estimates
  bottom <- best$estimate[1]
  top <- best$estimate[2]
  xmid<-best$estimate[3]
  scal <- best$estimate[4]
  s <- best$estimate[5]
  
  # Estimation des valeurs
  newX <- seq(min(x)*0.75, max(x)*1.2, length=500)
  newY <- .PL5(bottom, top, xmid, scal, s, newX)
  yfit <- .PL5(bottom, top, xmid, scal, s, x)
  perf <- .getPerf(y, yfit)
  
  # coordinates of inflexion point
  Xflex = xmid + (1/scal)*log10(s)
  Yflex = bottom + (top - bottom)*(s/(s+1))^s
  
  # coordinates of rep = 0.5
  Y50 = 0.5*(top + bottom)
  X50 = xmid - 1/scal*log10(((top - bottom)/(Y50 - bottom))^(1/s)-1)
  
  # Slope at the inflexion point
  B = (top - bottom)*log(10)*(scal)*(s/(s+1))^(s+1)	# + bottom 							# 5P	 	pour d=1, (d/(d+1))^(d+1) = 1/4
  A = Yflex  - B*(Xflex)																					# 5P		pour d=1, (d/(d+1))^d = 1/2
  
  # Visualization
  if(Plot){
    plot(y ~ x, type = 'n',...)
    pas = (max(x) - min(x))/20
    if(!is.null(sd))
      for(s in 1:length(sd)){
        segments(x0 = x[s], x1 = x[s], y0 = y[s]-sd[s], y1 = y[s]+sd[s], lty = 3)
        segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]-sd[s], y1 = y[s]-sd[s], lty = 3)
        segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]+sd[s], y1 = y[s]+sd[s], lty = 3)
      }
    if(addPoints) {points(x, y, pch=19, col = pCol,...); points(x, y, cex=1.2)}
    if(addLine) lines(newY ~ newX, col = lCol, lwd=3,...)
    if(addXinf) {abline(A, B, col = 'blue', lwd = 2); points(-A/B, 0, pch = 19, col = 'purple')}
    if(wcoef==0) Sub = "Non weighted 5P logistic regr."
    title (sub = Sub)
  }
  
  return(list(pars=best$estimate, goodness=perf$goodness, err=perf$err, p=perf$p,
              fit=yfit, newX=newX, newY=newY,
              xflex=Xflex, yflex=Yflex, x50=X50, y50=Y50,
              slope=B, yIntercept=A, xInfl=-A/B))
  
}
