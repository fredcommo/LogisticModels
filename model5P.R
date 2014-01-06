# model5P <- function(x, y, sd=NULL, W.coef=0.25, Plot=FALSE, Title="", Xlab="", Ylab="", pCol = "indianred2", lCol = "royalblue3", Sub = "Weighted 5P logistic regr.", ...){
 model5P <- function(x, y, sd=NULL, W.coef=0.25, Plot=TRUE, addLine = TRUE, addPoints = TRUE, addXinf = TRUE, pCol = "indianred2", lCol = "royalblue3", Sub = "Weighted 5P logistic regr.", ...){

# Avec nlm
	# Fonction logistique 5P
		nlm.fit.5P <- function(bottom, top, xmid, scal, s,  X){
			Y.fit <-bottom+(top-bottom)/(1+10^((xmid-X)*scal))^s
			return(Y.fit)
			}

	# Fonction sce : weighted sum squared error
		sce.5P <- function(param, X, yobs, Weights, W.coef) {
			bottom <- param[1]
			top = param[2]
			xmid = param[3]
			scal = param[4]
			s = param[5]
			ytheo <- nlm.fit.5P(bottom, top, xmid, scal, s, X)
			residus <- yobs - ytheo
			Weights <- (1/(residus^2))^(W.coef)
			return(sum(Weights*(yobs - ytheo)^2))
			}

	# initialisation des valeurs
		z = NULL
		w.coef = W.coef
		bottom.ini = min(y)
		top.ini = max(y)
		xmid.ini = (max(x)+min(x))/2
		z <- (y - bottom.ini)/(top.ini - bottom.ini)
		z[z==0] <- 0.01
		z[z==1] <- 0.99
		scal.ini = coef(lm(x~log(z/(1-z))))[2]
		scal.ini <- as.numeric(scal.ini)
		s.ini = 1
		# weights <- (z*(1-z))^(w.coef)
		# weights <- weights/sum(weights)
		weights <- rep(1, length(y))		

		init <- c(bottom.ini, top.ini, xmid.ini, scal.ini, s.ini)
		best<-nlm(f = sce.5P, p = init, X = x, yobs = y, Weights=weights, W.coef=w.coef)						# estimation of parameters

	# Best estimates
		best.bottom <- best$estimate[1]
		best.top <- best$estimate[2]
		best.xmid<-best$estimate[3]
		best.scal <- best$estimate[4]
		best.s <- best$estimate[5]

	# Estimation des valeurs
		newX <- seq(min(x)*0.75, max(x)*1.2, length=500)						
		Y.best <- nlm.fit.5P(best.bottom, best.top, best.xmid, best.scal, best.s, newX)
		y.fit <- nlm.fit.5P(best.bottom, best.top, best.xmid, best.scal, best.s, x)
		lm.test <- lm.rob(y, y.fit)
		p.fit <- coef(summary(lm.test$model))[2,4]

			# coordinates of inflexion point
				Xflex = best.xmid + (1/best.scal)*log10(best.s)
				Yflex = best.bottom + (best.top - best.bottom)*(best.s/(best.s+1))^best.s

			# coordinates of rep = 0.5
				Y50 = 0.5*(best.top + best.bottom)
				X50 = best.xmid - 1/best.scal*log10(((best.top - best.bottom)/(Y50 - best.bottom))^(1/best.s)-1)

			# Slope at the inflexion point
				B = (best.top - best.bottom)*log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1)	# + best.bottom 							# 5P	 	pour d=1, (d/(d+1))^(d+1) = 1/4
				A = Yflex  - B*(Xflex)																					# 5P		pour d=1, (d/(d+1))^d = 1/2

	# Visualization
		if(Plot){
			plot(y ~ x, type = 'n',...)
			#xlim = range(min(newX), max(newX)), ylim = range(min(y, Y.best),max(y, Y.best)),...)
			pas = (max(x) - min(x))/20
			if(!is.null(sd))
				for(s in 1:length(sd)){
					segments(x0 = x[s], x1 = x[s], y0 = y[s]-sd[s], y1 = y[s]+sd[s], lty = 3)
					segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]-sd[s], y1 = y[s]-sd[s], lty = 3)
					segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]+sd[s], y1 = y[s]+sd[s], lty = 3)
					}
				if(W.coef==0) Sub = "Non weighted 5P logistic regr."
				title (sub = Sub)
				}
			if(addPoints) points(x, y, col = pCol,...)
			if(addLine) lines(Y.best~newX, col = lCol,...)
			if(addXinf) {abline(A, B, col = 'blue', lwd = 2); points(-A/B, 0, pch = 19, col = 'purple')}

			# legend("bottomleft", legend = paste("IC50 :", round(10^X50, 2),unit), bty="n")

	# return(list(lo.fit = best.bottom, hi.fit = best.top, fit.values = y.fit, p.value = p.fit))
	return(list(model = best, Xfit = newX, fitted.values = Y.best, First.fit = y.fit[1], Last.fit = y.fit[length(x)], p.value = p.fit,
						xflex = Xflex, yflex = Yflex, x50 = X50, y50 = Y50, slopeInfl = B, yIntercept = A, xInfl = -A/B))
}
