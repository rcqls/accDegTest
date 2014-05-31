# model: g(D_k)=a_0+a_1*Af(x_k)*t_k+U_k
# with Af(x)=exp(b^T %*% (z(x_ref)-z(x)))
# input=data.frame D(égradation)_k,x(Stress)_k,t(emps)_k
# given: g and z, x_ref
# param: a_0,a_1,b,s_U 


# Version (New) based on the maximum-likekihood equivalent to the least square method!!!
# addt(d ~ t | 1/x & 1/x_2 ,data,xref=?,transf=log)
addt <- function(formula,xref=1,transf=log,data) {
  ## simple transf 
  if(length(transf)==1) {
    if(identical(transf,log)) transf <- c(log,exp)
    else if(identical(transf,exp)) transf <- c(exp,log)
    else if(identical(transf,identity)) transf <- c(identity,identity) #no transform
    else stop("transf argument generally needs to be of length 2.")
  }

  obj <- new.env()
  obj$formula <- formula
  # remove 
  obj$with_intercept <- attr(terms(as.formula(strsplit(deparse(formula),"\\|")[[1]][1])),"intercept")
  obj$xref<-xref
  obj$transf<-transf

  class(obj) <- "addt" 

  if(missing(data)) init.addt(obj) else init.addt(obj,data=data)

  prepare.estim.addt(obj)

  # first estimation initialized with the addtTwoPass method
  addt2 <- addtTwoPass(formula=formula,xref=xref,transf=transf,data=data)
  run(obj,coef(addt2)[-(1:2)])
  # save coef for free acceleration model 
  obj$coef4freeAccelModel <- list(intercept=addt2$intercept.1,slope=addt2$slope.1)

  obj
}

names.addt <- function(obj) {
  c("optim","par","formula","xref","transf")
} 

run.addt <- function(obj,par0,method=c("BFGS","CG","L-BFGS-B","Nelder-Mead","SANN"),mode=c("two.steps","one.step","best"),verbose=FALSE,...) {
  method <- match.arg(method)
  mode <- match.arg(mode)
  if(is.null(obj$par)) obj$par<-rep(0,ncol(obj$zmodel))
  if(missing(par0)) par0 <- obj$par
  if(mode=="two.steps") {
    if(length(par0)==1) obj$optim <- list(par=par0) else obj$optim<-optim(par0,obj$ols)
    obj$optim<-optim(obj$optim$par,obj$ols,obj$ols.gr,method=method,...)
  } 
  else if(mode=="one.step") obj$optim<-optim(par0,obj$ols,obj$ols.gr,method=method,...)
  
  else if(mode=="best") {
    i <-0
    res <- obj$par
    repeat {
      i<-i+1
      old <- res
      res<-run(obj,mode="two.steps",verbose=FALSE)
      if(identical(all.equal(old,res),TRUE)) break
    }
    cat(i,"steps done!\n")
  }

  if(verbose) print(obj$optim)
  obj$par <- obj$optim$par
  coef(obj)
}


coef.addt <- function(obj,type=c("default","acceleration","af","AF")) {
  type <- match.arg(type)
  switch(type,
    af=,AF=,acceleration=as.vector(exp(as.matrix(obj$zmodel[obj$rank$start,]) %*% obj$par)),
    {
      if(!is.null(obj$optim) &&  obj$optim$convergence!=0) {
        warning("Non convergence of optimisation method!")
        print(obj$optim)
      }
      if(obj$with_intercept) {
        res<-c(obj$a0(obj$par),obj$a1(obj$par),obj$par)
        names(res)[1:2] <- c("(intercept)",obj$varnames$t)
      } else {
        res<-c(obj$a1(obj$par),obj$par)
        names(res)[1] <- c(obj$varnames$t)
      } 
      res
    }
  )
}

prepare.estim.addt <- function(obj) {
  # to get expression more explicit
  y <- obj$transf[[1]](obj$model[[1]])
  t <- obj$model[[2]]
  zz<- as.matrix(obj$zmodel)

  # Cache to theoretically increase computation
  #obj$x <- NULL;obj$dx <- NULL

  # Vector x with element AF(x_i;b)*t_i
  x <-function(b) {#if(is.null(obj$x)) (obj$x <<- 
    as.vector(exp(zz%*%b)*t)
    #) else obj$x
    }
  # slope estimation depending on b
  if(obj$with_intercept) {
    obj$a1<-function(b) {
      as.vector(cov(y,x(b))/var(x(b)))
    }
    # contrast to minimize
    obj$ols<-function(b) {
      ##obj$x <<- NULL # to be updated when first called 
      sum((y-mean(y)-obj$a1(b)*(x(b)-mean(x(b))))^2)
    }

    # function to export results
    obj$a0<-function(b) {
      #obj$x <<- NULL 
      mean(y)-obj$a1(b)*mean(x(b))
    }
  } else {
    obj$a1<-function(b) {
      as.vector(mean(y*x(b))/mean(x(b)^2))
    }
    # contrast to minimize
    obj$ols<-function(b) {
      ##obj$x <<- NULL # to be updated when first called 
      sum((y-obj$a1(b)*(x(b)))^2)
    }
    # function to export results
    obj$a0<-function(b) 0
  }

  # gradients
  x.gr <- function(b) #if(is.null(obj$dx)) (obj$dx<<-
            matrix(as.vector(zz)*x(b),nrow=nrow(zz))
          #) else obj$dx

  if(obj$with_intercept) {
    a1.gr<-function(b) apply(x.gr(b),2,function(d) (cov(y,d*var(x(b))-cov(y,x(b))*2*cov(x(b),d)))/var(x(b))^2)
   
    obj$ols.gr<-function(b) {
      #obj$x <<- NULL;obj$dx<<-NULL # to be updated when first called    
      dx <- x.gr(b)
      da1<- a1.gr(b)
      eps<-y-mean(y)-obj$a1(b)*(x(b)-mean(x(b)))
      sapply(1:length(da1),function(j) -2*sum(eps*(da1[j]*(x(b)-mean(x(b)))+obj$a1(b)*(as.vector(dx[,j])-mean(as.vector(dx[,j]))))))
    }

  } else {

    a1.gr<-function(b) apply(x.gr(b),2,function(d) (mean(y*d*mean(x(b)^2)-mean(y*x(b))*2*mean(x(b)*d)))/mean(x(b)^2)^2)
   
    obj$ols.gr<-function(b) {
      #obj$x <<- NULL;obj$dx<<-NULL # to be updated when first called    
      dx <- x.gr(b)
      da1<- a1.gr(b)
      eps<-y-obj$a1(b)*x(b)
      sapply(1:length(da1),function(j) -2*sum(eps*(da1[j]*(x(b))+obj$a1(b)*as.vector(dx[,j]))))
    }
  }
  
  obj
}


init.addt <- function(obj,...) {
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
  ## sort model by xx
  obj$model <- obj$model[order(as.character(obj$model$xx)),]
  ## rank for class
  obj$rank <- list(start=match(unique(obj$model$xx),obj$model$xx))
  obj$rank$end <- c(obj$rank$start[-1]+1,length(obj$model$xx))

  ## xref = usual stress transformed as a named list used as an environment
  xref <- obj$xref
  if(!is.list(xref)) {
    if(is.null(names(xref))) {
      if(length(xref)<length(varnames$x)) xref <- rep(xref,length=length(varnames$x))
      xref <- as.list(xref)
      names(xref) <- varnames$x
    } else if(all(sort(names(xref))==sort(varnames$x))) xref <- as.list(xref)
    else stop("Vector xref needs to be completed!")
  } else if(any(sort(names(xref))!=sort(varnames$x))) {
    stop("List xref needs to be completed!")
  }

  ## list used as environment for reference stress values
  obj$xref <- xref
  ## zexprs = expressions of the terms in the acceleration
  obj$zexprs <- strsplit(deparse(obj$formula[[length(obj$formula)]][[3]]),"&")[[1]]

  ## z = the values of the zcodes applied to the data as a data frame
  obj$zmodel <- stress.expression.addt(obj,obj$model)

  obj #return the obj itself
}

stress.expression.addt <- function(obj,df) {
  zdf <- as.data.frame(sapply(obj$zexprs,function(code) {
    zexpr <- parse(text=code)[[1]]
     eval(zexpr,envir=obj$xref)-eval(zexpr,envir=df)
  },simplify=FALSE))
  colnames(zdf) <- paste("z",1:ncol(zdf),sep="")
  zdf
}

model.frame.addt <- function(obj,...) {
  if(is.null(obj$model)) {
    data <- list(...)$data # NULL if data not provided inside ...
    fformula <- obj$formula
    fformula[[length(fformula)]][[1]] <- as.name("+")
    fformula[[length(fformula)]][[3]] <- parse(text=paste(obj$varnames$x,collapse=" + "))[[1]]
    if(is.null(data)) model.frame(fformula) else model.frame(fformula,data=data)
    } else obj$model
}

predict.addt <- function(obj,to.predict=c("degradation","acceleration","af"),...,data,plot=FALSE,transf=TRUE) {
  to.predict <- match.arg(to.predict)
  # data.frame of values to predict
  df <- if(missing(data)) {
     if(length(list(...))>0) data.frame(...)[intersect(c(obj$varnames$t,obj$varnames$x),names(list(...)))]
     else obj$model
  } else data

  # if(is.null(df$xx)) {
  #   df$xx <- as.factor(apply(df[obj$varnames$x],1,function(l) paste(paste(names(l),l,sep="="),collapse="&")))
  # }

  coeff <- coef(obj)

  #print(df)
  zdf <- stress.expression.addt(obj,df)
  #print(zdf)
  acc <- exp(as.vector(as.matrix(zdf) %*% coeff[-(1:2),drop=FALSE]))
  if(to.predict %in% c("acceleration","degradation")) acc <- acc*coeff[2]
  if(to.predict=="degradation") {
    res <- coeff[1]+acc*df[[obj$varnames$t]]
    if(transf) res <- obj$transf[[2]](res)
    if(plot==T) points(obj$model[[obj$varnames$t]],res)
    res
  } else {
    if(plot && length(df)==1) {
      plot(df[[1]],acc,type="l",main=paste("acceleration vs ",obj$varnames$x,sep=""),xlab=obj$varnames$x,ylab="acceleration")
      invisible(acc)
    } else acc
  }
}


residuals.addt <- function(obj) {
  # if(!missing(coeff)) {
  #   resid <- obj$model[[1]]
  #   for(i in seq(obj$rank$start)) {
  #     sub <- obj$rank$start[i]:obj$rank$end[i]
  #     resid[sub] <- obj$model.1[[1]][sub]-(coeff[1]+ coeff[1+i]*obj$model.1[[2]][sub]) 
  #   }
  #   resid
  # } else if(step==1) {
  #   if(obj$method==1) {
  #     resid <- obj$model[[1]]
  #     for(i in seq(obj$rank$start)) {
  #       sub <- obj$rank$start[i]:obj$rank$end[i]
  #       resid[sub] <- obj$transf[[1]](obj$model[[1]][sub])-(obj$a0(obj$par)+obj$slope.1[i]*obj$model.1[[2]][sub]) 
  #     }
  #     resid
  #   } else residuals(obj$lm.1)
  # } else {
    resid <- obj$model[[1]]
    aa1 <- obj$a1(obj$par)*coef(obj,"AF")
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      resid[sub] <- obj$transf[[1]](obj$model[[1]][sub])-(obj$a0(obj$par)+ aa1[i]*obj$model[[2]][sub]) 
    }
    resid
  # }
}

ftq.addt <- function(obj,threshold,p=.95,tail=FALSE,...,plot=FALSE) { #... is for the 
  df <- data.frame(...)[obj$varnames$x]
  sigma <- sqrt(sum(residuals(obj)^2)/(nrow(obj$model)-1))
  #print(list(obj$model,sigma,coef(obj)[1],qnorm(if(tail) 1-p else p)*sigma,obj$transf[[1]](threshold)))
  ftq <- (obj$transf[[1]](threshold)-coef(obj)[1]-qnorm(if(tail) 1-p else p)*sigma)/predict(obj,"accel",data=df)
  names(ftq) <- NULL
  if(plot && length(df)==1) {
    plot(df[[1]],ftq,type="l",main=paste("failure time quantile (threshold=",threshold,",p=",p*100,"%)",sep=""),xlab=obj$varnames$x,ylab=obj$varnames$t)
    invisible(ftq)
  } else ftq

}

ftm.addt <- function(obj,threshold,...,plot=FALSE) { #... is for the 
  df <- data.frame(...)[obj$varnames$x]
  ftm <- (obj$transf[[1]](threshold)-coef(obj)[1])/predict(obj,"accel",data=df)
  names(ftm) <- NULL
  if(plot && length(df)==1) {
    plot(df[[1]],ftm,type="l",main=paste("failure time mean (threshold=",threshold,")",sep=""),xlab=obj$varnames$x,ylab=obj$varnames$t)
    invisible(ftm)
  } else ftm

}

print.addt <- function(obj, ... ) {
  print(coefficients(obj))
}

summary.addt <- function(obj,...) {
  r2 <- if(obj$with_intercept) 1-var(resid(obj))/var(obj$transf[[1]](obj$model[[1]]))
        else 1-mean(resid(obj)^2)/mean((obj$transf[[1]](obj$model[[1]]))^2)
  list(
    coefficients=coef(obj),
    sigma=sqrt(mean((resid<-residuals(obj))^2)/(length(resid)-2)*length(resid)),
    r.squared=r2
  )
}

# inverse <- function (f, lower = -100, upper = 100) {
#    function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
# }

linesFreeAccelModel.addt <- function(obj,args) {
  for(i in seq(obj$rank$start))
    abline(a=mean(obj$coef4freeAccelModel$intercept),b=obj$coef4freeAccelModel$slope[i],col=args$col[i],lty=args$lty[i],lwd=args$lwd[i])    
}

curvesFreeAccelModel.addt <- function(obj,args) {
  for(i in seq(obj$rank$start))
    curve(obj$transf[[2]](mean(obj$coef4freeAccelModel$intercept)+obj$coef4freeAccelModel$slope[i]*x),col=args$col[i],lty=args$lty[i],lwd=args$lwd[i],add=TRUE)      
}

plot.addt <- function(obj,type="all degradations",with.layout=TRUE,fit=TRUE,only=NA,fitFreeAccel=FALSE,xlim=NULL,ylim=NULL,...) {

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

  aa1<-obj$a1(obj$par)*coef(obj,"AF")
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
        curve(obj$transf[[2]](obj$a0(obj$par)+aa1[i]*x),col=args$col[i],lty=args$lty[i],lwd=args$lwd[i],add=TRUE)
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
      abline(a=obj$a0(obj$par),b=aa1[i],col=args$col[i],lty=args$lty[i],lwd=args$lwd[i])
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

lines.addt <- function(obj,only=NA,method=c("default","free.accel"),ic=NULL,...) {
  method <- match.arg(method)
   

  xx.uniq <- obj$model$xx[obj$rank$start]
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

  aa1<-obj$a1(obj$par)*coef(obj,"AF")

  switch(method,
    free.accel={
      #TODO! pb with a0 what is the best estimation method?
    },{
      for(i in seq(obj$rank$start)) {
        abline(a=obj$a0(obj$par),b=aa1[i],col=args$col[i],lty=args$lty[i],lwd=args$lwd[i])
      }
    }
  )

}



#TODO: random effect model dealt in two steps with Y_i,1-Y_i,0 (sigma) and Y_i,0 (sigma^(0))
# extension; Y_i,j-hat{a}_i and hat{a}^(0)_i when j>=1 (if j=1 hat{a}^(0)_i=Y_i,0)


