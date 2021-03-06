# clouds(d ~ t | x & x_2 ,data)
# N.B.: plot and lines have a transf parameter to change y variable!

clouds <- function(formula,data,measure.varname="measure",system.varname="system") {
  ## simple transf 
  obj <- new.env()
  obj$formula <- formula
  class(obj) <- "clouds"
  if(missing(data)) {
   init.clouds(obj) 
  } else { # for measures and systems identification
  	init.clouds(obj,data=data)
  	obj$data <- data
  	init.system.clouds(obj,measure.varname,system.varname)
  }
  obj
}

init.clouds <- function(obj,...) {
  ## varnames 
  varnames<-list(
    y=all.vars(obj$formula[[2]]),
    t=all.vars(obj$formula[[length(obj$formula)]][[2]]),
    x=all.vars(obj$formula[[length(obj$formula)]][[3]])
  )
  ## environment for the first regression
  obj$varnames <- varnames  # stressed varnames
  obj$model<-model.frame(obj,...)
  ## xx = unique identifier for the level of stress
  obj$model$xx <- as.factor(apply(obj$model[obj$varnames$x],1,function(l) paste(paste(names(l),l,sep="="),collapse="&")))

  obj #return the obj itself
}

init.system.clouds <- function(obj,measure.varname="measure",system.varname="system") {
	# detect system variable
	if(is.character(system.varname)) {
		if(system.varname %in% names(obj$data)) obj$varnames$system <- system.varname
		else system.varname <- 1
	}
	if(is.numeric(system.varname)) obj$varnames$system <- names(data)[[system.varname]]
	if(obj$varnames$system %in% all.vars(obj$formula)) stop("system variable name matches to some variable of the model")
	
	# detect measure variable
	if(is.character(measure.varname)) {
		if(measure.varname %in% names(obj$data)) obj$varnames$measure <- measure.varname
		else measure.varname <- 2
	}
	if(is.numeric(measure.varname)) obj$varnames$measure <- names(data)[[measure.varname]]
	if(obj$varnames$measure %in% all.vars(obj$formula)) stop("measure variable name matches to some variable of the model")
	
	if(!is.numeric(obj$data[[obj$varnames$measure]])) stop("measure variable is supposed to be numeric!")
	obj$syst_id <- paste(obj$data[[obj$varnames$system]],apply(obj$data[obj$varnames$x],1,paste,collapse="&"),sep="&")
	obj$measure <- obj$data[[obj$varnames$measure]]
}

model.frame.clouds <- function(obj,...) {
  if(is.null(obj$model)) {
    data <- list(...)$data # NULL if data not provided inside ...
    fformula <- obj$formula
    fformula[[length(fformula)]][[1]] <- as.name("+")
    fformula[[length(fformula)]][[3]] <- parse(text=paste(obj$varnames$x,collapse=" + "))[[1]]
    if(is.null(data)) model.frame(fformula) else model.frame(fformula,data=data)
  } else obj$model
}


plot.clouds <- function(obj,only=NA,method=c("default","measures"),xlim=NULL,main=NULL,transf=NULL,arrows.only=NA,...) {

	method <- match.arg(method)
	if(is.null(transf)) transf <- function(x) x
	xx <- obj$model$xx
	xx.uniq <- unique(xx)

	args <- list(...)
	#default value
	if(is.null(args$lwd)) args$lwd <- 2 
	if(is.null(args$lty)) args$lty <- 1
	if(is.null(args$col)) args$col <- seq(xx.uniq)+1 
	if(is.null(args$pch)) args$pch <- seq(xx.uniq)
	## complete the vector to the proper length in order to extract index
	for(v in c("lty","lwd","col","pch")) args[[v]] <- rep(args[[v]],length.out=length(xx.uniq)) 

	do.legend <- function(args,...) {
	default <- list(...)
	for(v in names(default)) if(is.null(args[[v]])) args[[v]] <- default[[v]]
	if(!is.null(args$pos)) {
	  pos <- args$pos
	  args$pos <- NULL
	}
	legend <- list(x=pos[1])
	if(length(pos)==2) legend$y <- pos[2]
	legend <- c(legend,args)
	if(is.null(legend$legend)) legend$legend <- xx.uniq
	do.call("legend",legend)
	}

	if(length(only)>1 || !is.na(only)) {
		args$col[setdiff(seq(xx.uniq),only)] <- 0
	}  

	symbols <- c(letters,1:100) # TODO

	plot(obj$model[[2]],transf(obj$model[[1]]),xlab=obj$varnames$t,ylab=obj$varnames$y,main=if(is.null(main)) "clouds" else main,col=0,xlim=xlim)
	switch(method,
	measures={
		for(i in seq(xx.uniq)) {
			cpt <- 0
		  	sub <- xx==xx.uniq[i]
		  	systs <- unique(obj$syst_id[sub])
		  	for(k in systs) {
		  		pch <- symbols[cpt <- cpt+1]
		  		syst <- obj$syst_id == k
		  		x<-obj$model[syst,2];y<-transf(obj$model[syst,1])
		  		if(length(arrows.only)>1 && (cpt %in% arrows.only[[i]])) {
		  			arrows(x[-length(x)],y[-length(y)],x[-1],y[-1],col=args$col[i])
		  		}
		  		points(x,y,col=args$col[i],pch=pch,lwd=args$lwd[i])
				
			}
		}
		args$pch <- NULL
		do.legend(args,pos="bottomleft",cex=.8)
	},{
		for(i in seq(xx.uniq)) {
		  	sub <- xx==xx.uniq[i]
		  	points(obj$model[sub,2],transf(obj$model[sub,1]),col=args$col[i],pch=args$pch[i],lwd=args$lwd[i])
		}
		do.legend(args,pos="bottomleft",cex=.8)
	})

} 

lines.clouds <- function(obj,only=NA,method=c("default","same.intercept"),ic=NULL,transf=NULL,...) {
	if(is.null(transf)) transf <- function(x) x
	xx <- obj$model$xx
	xx.uniq <- unique(xx)

	method <- match.arg(method)

	args <- list(...)
	#default value
	if(is.null(args$alpha.col)) args$alpha.col <- .1 
	if(is.null(args$lwd)) args$lwd <- 2 
	if(is.null(args$lty)) args$lty <- 1
	if(is.null(args$col)) args$col <- seq(xx.uniq)+1 
	if(is.null(args$pch)) args$pch <- seq(xx.uniq)
	## complete the vector to the proper length in order to extract index
	for(v in c("lty","lwd","col","pch","alpha.col")) args[[v]] <- rep(args[[v]],length.out=length(xx.uniq)) 
  	
	if(length(only)>1 || !is.na(only)) {
		args$col[setdiff(seq(xx.uniq),only)] <- 0
	} 

	formula1 <- obj$formula
	formula1[[3]] <- formula1[[3]][[2]]

	switch(method,
		same.intercept={
			## ajout de la variable xx (niveau des stress)
		    formula1[[3]] <- parse(text=paste("xx",deparse(formula1[[3]]),sep=":"))[[1]]
		    data <- obj$model
		    data[[1]] <- transf(data[[1]])
		    lm.1 <-lm(formula1,data=data)
		    lm.1.coef <- if(attr(terms(formula1),"intercept")) lm.1$coef else c(0,lm.1$coef)
		    if(!is.null(ic)) {
			    lm.1.xlim <- range(c(0,obj$model[[2]])) 			#include intercept in range
			    lm.1.xlim <- lm.1.xlim + diff(lm.1.xlim)*c(-.8,1.2) #increase range
			    lm.1.new <- data.frame(seq(lm.1.xlim[1],lm.1.xlim[2],length=100))
			    names(lm.1.new) <- obj$varnames$t
			}
		    for(i in seq(xx.uniq)) {
			  	abline(a=lm.1.coef[1],b=lm.1.coef[1+i],col=args$col[i],pch=args$pch[i],lwd=args$lwd[i],lty=args$lty[i])
				if(!is.null(ic)) {
					lm.1.newdata <- cbind(lm.1.new,data.frame(xx=obj$model[which(xx==xx.uniq[i])[1],"xx"]))
		    		lm.1.pred <- predict(lm.1,newdata=lm.1.newdata,interval="confidence",level=1-ic)
					if(args$col[i]>0) polygon(c(lm.1.newdata[[obj$varnames$t]],rev(lm.1.newdata[[obj$varnames$t]])),c(lm.1.pred[,2],rev(lm.1.pred[,3])),col=do.call("rgb",as.list(c(col2rgb(args$col[i])/255,args$alpha.col[i]))),border=args$col[i],lwd=1,lty=args$lty[i])	
				}
			}

		#},
		# same.intercept={
		# 	## ajout de la variable xx (niveau des stress)
		#     formula1[[3]] <- parse(text=paste(as.character(formula1[[3]]),"xx",sep=":"))[[1]]
		#     lm.1 <-lm(formula1,data=obj$model)
		#     lm.1.coef <- lm.1$coef
		#     if(!is.null(ic)) {
		#     	lm.1.se <- (tmp<-summary(lm.1))$coef[-1,2] + tmp$coef[1,2]
		#     	lm.1.side <- c(-1,-1,1,1)*qnorm(1-ic/2)
		# 	    lm.1.xlim <- range(c(0,obj$model[[2]]))
		# 	    lm.1.xlim <- (lm.1.xlim + diff(lm.1.xlim)*c(-.5,1.5))[c(1,2,2,1)] 
		# 	}
		#     for(i in seq(xx.uniq)) {
		# 	  	abline(a=lm.1.coef[1],b=lm.1.coef[1+i],col=args$col[i],pch=args$pch[i],lwd=args$lwd[i],lty=args$lty[i])
		# 		if(!is.null(ic)) {
		# 			lm.1.ylim <- lm.1.coef[1]+lm.1.coef[1+i]*lm.1.xlim + lm.1.side*lm.1.se[i]
		# 			polygon(lm.1.xlim,lm.1.ylim,col=do.call("rgb",as.list(c(col2rgb(args$col[i])/255,.1))),border=args$col[i],lwd=1,lty=args$lty[i])	
		# 		}
		# 	}

		},{
		  	
		  	for(i in seq(xx.uniq)) {
			  	sub <- xx==xx.uniq[i]
			  	abline(lm(formula1,data=obj$model[sub,]),col=args$col[i],pch=args$pch[i],lwd=args$lwd[i])
			}
		}
	)
}

