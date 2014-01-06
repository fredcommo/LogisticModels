 modelp5P <- function(x, y, sd=NULL, W.coef=0, Plot=FALSE, Title="", Xlab="", Ylab="", ptCol = "indianred2", lCol = "royalblue3"){

# Avec nlm
	# Fonction logistique 5P
		nlm.fit.p5P <- function(top, xmid, scal, s,  X){
			Y.fit <- top/(1+10^((xmid-X)*scal))^s
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus) avec pondérations
		sce.p5P <- function(param, X, yobs, Weights, W.coef) {
			top <- param[1]
			xmid <- param[2]
			scal <- param[3]
			s <- param[4]
			ytheo <- nlm.fit.p5P(top, xmid, scal, s, X)
			residus <- yobs - ytheo
			Weights <- (1/(residus^2))^(W.coef)
			return(sum(Weights*(yobs - ytheo)^2))
			}

	# initialisation des valeurs
		z = NULL
		w.coef = W.coef
		top.ini = max(y)
		xmid.ini = (max(x)+min(x))/2
		z <- y /top.ini
		z[z==0] <- 0.01
		z[z==1] <- 0.99
		scal.ini = coef(lm(x~log(z/(1-z))))[2]
		scal.ini <- as.numeric(scal.ini); scal.ini
		s.ini = 1
		# weights <- (z*(1-z))^(w.coef)
		# weights <- weights/sum(weights)
		weights <- rep(1, length(y))		

		init <- c(top.ini, xmid.ini, scal.ini, s.ini)
		best<-nlm(f = sce.p5P, p = init, X = x, yobs = y, Weights=weights, W.coef=w.coef)						# calcul des paramètres

	# Récupération des paramètres
		best.top <- best$estimate[1]
		best.xmid<-best$estimate[2]
		best.scal <- best$estimate[3]
		best.s <- best$estimate[4]

	# Estimation des valeurs
		newX <- seq(min(x)*0.75, max(x)*1.2, length=500)						
		Y.best <- nlm.fit.p5P(best.top, best.xmid, best.scal, best.s, newX)
		y.fit <- nlm.fit.p5P(best.top, best.xmid, best.scal, best.s, x)
		lm.test <- lm.rob(y, y.fit)
		p.fit <- coef(summary(lm.test))[2,4]

			# coordonnées du pt d'inflexion
			#	Xflex = best.xmid + (1/best.scal)*log10(best.s)
			#	Yflex = best.bottom + (best.top - best.bottom)*(best.s/(best.s+1))^best.s

			# coordonnées du pt rep = 0.5
			#	Y50 = 0.5
			#	X50 = best.xmid - 1/best.scal*log10(((best.top - best.bottom)/(Y50 - best.bottom))^(1/best.s)-1)

			# pente au pt d'inflexion
			#	# B = best.bottom + (best.top - best.bottom)*log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1); B	# + best.bottom 							# 5P	 	pour d=1, (d/(d+1))^(d+1) = 1/4
			#	B = log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1); B								# top = 1, bottom = 0
			#	A = Yflex  - B*(Xflex); A													# 5P		pour d=1, (d/(d+1))^d = 1/2
			#	print(c(bottom = best.bottom, top = best.top,
			#			xmid = best.xmid, scal = best.scal, s = best.s, 
			#			Xflex = Xflex, Yflex = Yflex, X50 = X50, D50 = round(10^X50, 2)))

	# Représentations graphiques
		if(Plot){
			plot(y~x, pch = 20, cex = 2.5, col = ptCol, xlim = range(min(newX), max(newX)), xlab = Xlab,
				ylim = range(min(y, Y.best),max(y, Y.best)), ylab = Ylab)
			pas = (max(x) - min(x))/20
			if(!is.null(sd))
				for(s in 1:length(sd)){
					segments(x0 = x[s], x1 = x[s], y0 = y[s]-sd[s], y1 = y[s]+sd[s], lty = 3)
					segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]-sd[s], y1 = y[s]-sd[s], lty = 3)
					segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]+sd[s], y1 = y[s]+sd[s], lty = 3)
					}
			lines(Y.best~newX, col = lCol, lwd = 3)
			Sub = "Weighted p5P logistic regr."
			if(W.coef==0) Sub = "Non weighted p5P logistic regr."
			title (main = Title, sub = Sub)
			# legend("bottomleft", legend = paste("IC50 :", round(10^X50, 2),unit), bty="n")
		}
	# return(list(lo.fit = best.bottom, hi.fit = best.top, fit.values = y.fit, p.value = p.fit))
	# return(list(First.fit = y.fit[1], Last.fit = y.fit[length(x)], p.value = p.fit))
	return(list(model = best, fitted.values = y.fit))
}
