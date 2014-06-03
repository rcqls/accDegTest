# data is required here!
# measure==0 and measure==1 are supposed to be inside data
# measure in 1 position or named measure as a default

addt0 <- function(formula,xref=1,transf=log,data,measure.varname="measure",system.varname="system") {
	if(missing(data)) stop("data is missing!")

	## simple transf 
	if(length(transf)==1) {
		if(identical(transf,log)) transf <- c(log,exp)
		else if(identical(transf,exp)) transf <- c(exp,log)
		else if(identical(transf,identity)) transf <- c(identity,identity) #no transform
		else stop("transf argument generally needs to be of length 2.")
	}

	obj <- new.env()
	obj$data <- data
	obj$formula <- formula
	# remove 
	if(!attr(terms(as.formula(strsplit(deparse(formula),"\\|")[[1]][1])),"intercept")) stop("Intercept required since it is supposed to be a random effect.")
	obj$xref<-xref
	obj$transf<-transf

	class(obj) <- c("addt0","addt") # addt for model.frame!

	init.addt(obj,data=data)

	# reuse because already defined for general clouds plotting
	init.system.clouds(obj,measure.varname,system.varname)
	
	prepare.estim.addt0(obj) 
	
	obj
}

prepare.estim.addt0 <- function(obj) {
	# measure 0
	y0 <- obj$transf[[1]](obj$data[obj$measure==0,obj$varnames$y])
	#print(y0)
	obj$coef.y0 <- mean(y0)
	obj$sigma.y0 <- sd(y0)

	# measure 1
	data1 <- cbind(obj$data[obj$measure==1,], data.frame(dy=sapply(unique(obj$syst_id),function(s) {
			obj$transf[[1]](obj$data[obj$syst_id == s & obj$measure == 1, obj$varnames$y])
			- obj$transf[[1]](obj$data[obj$syst_id == s & obj$measure == 0, obj$varnames$y])
		}
	)))

	#print(dim(data1))
	form <- update(eval(as.call(c(obj$formula[[1]],obj$formula[[2]],obj$formula[[3]][[2]]))),.~0+.)
	formula <- obj$formula
	formula[[3]][[2]] <- form[[3]]
	obj$addt.dy <- addt(formula,obj$xref,transf=identity,data=data1)
	obj$coef.dy <- coef(obj$addt.dy)
	obj$sigma.dy <- summary(obj$addt.dy)$sigma

	## model coefficients
	obj$coef <- c(obj$coef.0,obj$coef.dy)
	names(obj$coef)[[1]] <- "intercept"
	obj$sigma <- obj$sigma.dy/sqrt(2)
	obj$sigma.0 <- sqrt(obj$sigma.y0^2 - obj$sigma^2) # random intercept effect noise 
}


coef.addt0 <- function(obj,type=c("default","acceleration","af","AF")) {
  type <- match.arg(type)
  switch(type,
    af=,AF=,acceleration=coef(obj$addt.dy,"AF"),
    obj$coef
  )
}

summary.addt0 <- function(obj) list(coefficients=coef(obj),sigma0=obj$sigma.0,sigma=obj$sigma,r.squared=summary(obj$addt.dy)$r.squared)

residuals.addt0 <- function(obj) {
    resid <- obj$model[[1]]
    aa1 <- obj$addt.dy$a1(obj$par)*coef(obj,"AF")
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      resid[sub] <- obj$transf[[1]](obj$model[[1]][sub])-(obj$coef.y0 + aa1[i]*obj$model[[2]][sub]) 
    }
    resid
  # }
}

plot.addt0 <- function(obj,type="all degradations",with.layout=TRUE,fit=TRUE,only=NA,fitFreeAccel=FALSE,xlim=NULL,ylim=NULL,...) {

  if(is.numeric(type)) type <- switch(type,"degradation","g-degradation","all degradations","interval stress plot","time resid","stress resid","normal","resid","points clouds")
  else type <- match.arg(type,c("degradation","g-degradation","all degradations","interval stress plot","time resid","stress resid","normal","resid","points clouds"))

  if(type %in% c("degradation","g-degradation","all degradations","stress plot","points clouds") && is.null(xlim) ) 
    xlim <- range(c(0,obj$model[[2]]))
  
  xx.uniq <- obj$model$xx[obj$rank$start]

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

  aa1<-obj$addt.dy$a1(obj$addt.dy$par)*coef(obj$addt.dy,"AF")
  #cat("aa1->");print(aa1)

  if(type=="points clouds") {
    plot.clouds(obj,xlim=xlim,main="points clouds",only=only,...)
  } else if(type=="degradation") {
    plot(obj$model[[2]],obj$model[[1]],xlab=obj$varnames$t,ylab=obj$varnames$y,main="degradation vs time",xlim=xlim,ylim=ylim)
    
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      points(obj$model[sub,2],obj$model[sub,1],col=args$col[i],pch=args$pch[i],lwd=args$lwd[i])
      if(fit) {
        #if(step %in% c(1,1.2,12,2.1,21)) 
        curve(obj$transf[[2]](obj$coef.y0+aa1[i]*x),col=args$col[i],lty=args$lty[i],lwd=args$lwd[i],add=TRUE)
        #if(step %in% c(2,1.2,12,2.1,21)) curve(obj$transf[[2]](mean(obj$intercept.1)+obj$slope.2[i]*x),col=args$col[i],lty=args$lty[i],lwd=args$lwd[i],add=TRUE)
        if(fitFreeAccel) {
          argsFree <- args
          argsFree$lty <- argsFree$lty*2
          curvesFreeAccelModel.addt(obj,argsFree)
        }
      }
    }
    do.legend(args,pos="bottomleft",cex=.8)
  } else if(type=="g-degradation") {
    plot(obj$model[[2]],obj$transf[[1]](obj$model[[1]]),xlab=obj$varnames$t,ylab=paste("g(",obj$varnames$y,")",sep=""),main="g-degradation vs time",xlim=xlim,ylim=ylim )
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      points(obj$model[sub,2],obj$transf[[1]](obj$model[sub,1]),col=args$col[i],pch=args$pch[i],lwd=args$lwd[i])
      # if(step %in% c(1,1.2,12,2.1,21)) abline(a=mean(obj$intercept.1),b=obj$slope.1[i],col=args$col[i],lty=args$lty[i]*lty.factor,lwd=args$lwd[i])
      # if(step %in% c(2,1.2,12,2.1,21)) abline(a=mean(obj$intercept.1),b=obj$slope.2[i],col=args$col[i],lty=args$lty[i],lwd=args$lwd[i])
      abline(a=obj$coef.y0,b=aa1[i],col=args$col[i],lty=args$lty[i],lwd=args$lwd[i])
      if(fitFreeAccel) {
          argsFree <- args
          argsFree$lty <- argsFree$lty*2
          linesFreeAccelModel.addt(obj,argsFree)
      }
    }
    do.legend(args,pos="bottomleft",cex=.8)
  } else if(type=="all degradations") {
    layout(matrix(c(1,2), 1,2, byrow = TRUE))
    plot(obj,"degradation",xlim=xlim,...)
    plot(obj,"g-degradation",xlim=xlim,...)
    layout(1)
  } else if(type=="interval stress plot") {
     plot.clouds(obj,xlim=xlim,main="stress plot",only=only,transf=obj$transf[[1]],...)
     lines.clouds(obj,method="same.intercept",ic=0.05,lty=2,only=only,transf=obj$transf[[1]])
     lines(obj,only=only,ic=NULL)
  } else if(type=="time resid") {
    resid <- residuals(obj)
    # uncomment if standardized residuals
    #resid <- resid/sqrt(mean(resid^2)/(length(resid)-2)*length(resid))
    plot(obj$model[[2]],resid,xlab=obj$varnames$t,main="residuals vs time")
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      points(obj$model[sub,2],resid[sub],col=args$col[i],pch=args$pch[i],lwd=args$lwd[i])
      abline(h=0,lwd=3)
    }
    do.legend(args,pos="top",cex=.8)
  } else if(type=="stress resid") {
    resid <- residuals(obj)
    # uncomment if standardized residuals
    #resid <- resid/sqrt(mean(resid^2)/(length(resid)-2)*length(resid))
    nx <- length(obj$varnames$x)
    if(with.layout) layout(matrix(c(rep(1,nx),2:(nx+1)), nx,2, byrow = FALSE))
    plot(obj$model$xx,resid,main="residuals vs stress")
    if(with.layout) for(v in obj$varnames$x) {
      plot(obj$model[[v]],resid,xlab=v,main=paste("residuals vs",v))
    }
    if(with.layout) layout(1)
  } else if(type=="normal") {
    resid <- residuals(obj)
    # uncomment if standardized residuals
    #resid <- resid/sqrt(mean(resid^2)/(length(resid)-2)*length(resid))
    qqnorm(resid)
  } else if(type=="resid") {
    layout(matrix(c(1,1,2,3), 2,2, byrow = TRUE))
    plot(obj,"time resid",...)
    plot(obj,"stress resid",with.layout=FALSE)
    plot(obj,"normal")
    layout(1)
  } 

}
