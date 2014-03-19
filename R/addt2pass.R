# model: g(D_k)=a_0+a_1*Af(x_k)*t_k+U_k
# with Af(x)=exp(b^T %*% (z(x_ref)-z(x)))
# input=data.frame D(égradation)_k,x(Stress)_k,t(emps)_k
# given: g and z, x_ref
# param: a_0,a_1,b,s_U 


# Version (OLD) with two passes! Mainly used as an initialization of addt.
# addtTwoPass(d ~ t | 1/x & 1/x_2 ,data,xref=?,transf=log)
addtTwoPass <- function(formula,xref=1,transf=c(log,exp),data,method=4,weight.method3=function(x) 1/sqrt(sum(x^2)),more.step=0) {
  ## simple transf 
  if(length(transf)==1) {
    if(identical(transf,log)) transf <- c(log,exp)
    else if(identical(transf,exp)) transf <- c(exp,log)
    else stop("transf argument generally needs to be of length 2.")
  }

  obj <- list(formula=formula,xref=xref,transf=transf,method=method)
  class(obj) <- c("addtTwoPass","addt")

  obj <- if(missing(data)) init.addt(obj) else init.addt(obj,data=data)

  ###### First regressions
  ## formula1 = d ~ t => first regression to apply for each stress level
  formula1 <- obj$formula
  formula1[[3]] <- formula1[[3]][[2]]
  model1 <- obj$model
  model1[[1]] <- obj$transf[[1]](model1[[1]])
  names(model1) <- names(obj$model)
  obj$model.1 <- model1
  ##print(model1)

  if(method %in% c(1,3,4)) {
    ## apply on the unique(xx) and compute the first regression for each level of stress
    lm1 <- list()
  
    res1 <- sapply(seq(obj$rank$start),function(i) {
      sub <-  sub <- obj$rank$start[i]:obj$rank$end[i]
      lm(formula1,data=obj$model.1[sub,]) -> tmp
      lm1[[i]] <<-tmp
      c(tmp$coeff,sum(residuals(tmp)^2)/tmp$df.residual)
    })
    
    obj$lm.1 <- lm1
    obj$intercept.1 <- res1[1,] #here, intercept is a vector
    obj$slope.1 <- res1[2,]
    obj$sigma2.1 <- res1[3,]

  } else if (method==2) {

    ## ajout de la variable xx (niveau des stress)
    formula1[[3]] <- parse(text=paste(as.character(formula1[[3]]),"xx",sep=":"))[[1]]
    obj$lm.1 <-lm(formula1,data=model1)
    obj$intercept.1 <- obj$lm.1$coef[1]
    obj$slope.1 <- obj$lm.1$coef[-1]
    obj$sigma2.1 <- sum(residuals(obj$lm.1)^2)/obj$lm.1$df.residual
  }
  # further step for method3
  if (method %in% c(3,4)) {
    if(method==3) {
      xmodel <- obj$model[obj$rank$start,obj$varnames$x] 
      weight <- apply(sapply(seq(obj$xref),function(i) xmodel[[i]]-xref[[i]]),1,function(x) weight.method3(x))
      obj$weight <- weight/sum(weight)
      obj$intercept.1 <- sum(obj$intercept.1*obj$weight)
    } else if(method==4) {
      if(all(obj$slope.1>0))  obj$intercept.1 <- min(obj$intercept.1)
      else if(all(obj$slope.1<0))  obj$intercept.1 <- max(obj$intercept.1)
      else warning("Error is coming because of method 1!")
    }
    ## ajout de la variable xx (niveau des stress)
    model1.3 <- model1
    model1.3[[1]] <- model1.3[[1]] - obj$intercept.1
    formula1 <- as.formula(paste(deparse(formula1),":xx + 0 ",sep=""))
    obj$lm.1 <-lm(formula1,data=model1.3)
    obj$slope.1 <- obj$lm.1$coef
    obj$sigma2.1 <- sum(residuals(obj$lm.1)^2)/obj$lm.1$df.residual
  }

  obj$formula.1 <-  formula1
  obj$sign <- if(all(obj$slope.1>0)) 1
  else if(all(obj$slope.1<0)) -1
  else {
    print(c(obj$intercept.1,obj$slope.1))
    stop("addtTwoPass model error: first regressions do not have same slope sign!")
  }
  ###### Second regression## obj$model2 = data.frame for the second regression
  model2 <- cbind(data.frame(y2=log(as.vector(obj$sign*obj$slope.1))),obj$zmodel[obj$rank$start,,drop=FALSE])
  rownames(model2) <- model1$xx[obj$rank$start]
  formula2 <- as.formula(paste("y2 ~ ",paste(colnames(obj$zmodel),collapse=" + "),sep=""))
  obj$lm.2 <-lm(formula2,data=model2)
  res2 <-obj$lm.2$coef
  obj$formula.2 <- formula2
  obj$model.2 <- model2
  ## final obj
  obj$coeff <- c(mean(obj$intercept.1),obj$sign*exp(res2[1]),res2[-1]) #TODO: need to be weighet????
  names(obj$coeff) <- c("(intercept)",obj$varnames$t,obj$zexprs)
  ## corrected slope after second step 
  obj$slope.2 <- exp(as.matrix(obj$zmodel[obj$rank$start,]) %*% obj$coeff[-(1:2)])*obj$coeff[2]
  # final step
  tmp <- more.a.step.addtTwoPass(obj)
  obj$slope.2 <- tmp$slope
  obj$coeff <- tmp$coefficients

  if(more.step>0) {
    ## another step one with fixed intercept
    obj$formula.1.i <- obj$formula
    obj$formula.1.i[[3]] <- obj$formula.1.i[[3]][[2]]
    obj$formula.1.i[[3]] <- parse(text=paste(as.character(obj$formula.1.i[[3]]),":xx+0",sep=""))[[1]]
    obj$formula.2.i <- as.formula(paste("y2 ~ 0 + ",paste(colnames(obj$zmodel),collapse=" + "),sep=""))
    
    ## 
    for(i in 1:more.step) {
      tmp <- more.b.step.addtTwoPass(obj)
      obj$lm.2 <- tmp$lm.2.i
      obj$coeff <- tmp$coefficients 
      tmp <- more.a.step.addtTwoPass(obj)
      obj$slope.2 <- tmp$slope
      obj$coeff <- tmp$coefficients
    }
  }
  obj
}

more.a.step.addtTwoPass <- function(obj,plot=FALSE) {
  t2 <- predict(obj,"af")*obj$model[[obj$varnames$t]]
  coeff <- c(lm(obj$model.1[[obj$varnames$y]]~t2)$coef,coef(obj)[-(1:2)])
  names(coeff) <- c("(intercept)",obj$varnames$t,obj$zexprs)
  slope <- coef(obj,"af")*coeff[2]
  names(slope) <- as.character(obj$model$xx[obj$rank$start])
  if(plot) sapply(slope,function(b) abline(a=coeff[1],b=b))
  list(coefficients=coeff,slope=slope)
}

more.b.step.addtTwoPass <- function(obj,plot=FALSE) {
  # with fixed intercept and common slope
  model.1.i <- obj$model.1
  model.1.i[[1]] <- model.1.i[[1]] - obj$coeff[1]
  slope.1.i <- lm(obj$formula.1.i,data=model.1.i)$coef

  if(!all(slope.1.i>0) && !all(slope.1.i<0)) {
    print(slope.1.i)
    stop("addtTwoPass model error: first regressions do not have same slope sign!")
  }

  model.2.i <- cbind(data.frame(y2=log(as.vector(slope.1.i)/obj$coeff[2])),obj$zmodel[obj$rank$start,,drop=FALSE])
  rownames(model.2.i) <- model.1.i$xx[obj$rank$start]
  
  lm.2.i<-lm(obj$formula.2.i,data=model.2.i)
  obj$coeff[-(1:2)]<-lm.2.i$coef
  names(obj$coeff)[-(1:2)] <- obj$zexprs
  list(coefficients=obj$coeff,lm.2.i=lm.2.i)
}

predict.addtTwoPass <- function(obj,to.predict=c("degradation","acceleration","af"),...,data,plot=FALSE,transf=TRUE) {
  to.predict <- match.arg(to.predict)
  # data.frame of values to predict
  df <- if(missing(data)) {
     if(length(list(...))>0) data.frame(...)[intersect(c(obj$varnames$t,obj$varnames$x),names(list(...)))]
     else obj$model
  } else data

  # if(is.null(df$xx)) {
  #   df$xx <- as.factor(apply(df[obj$varnames$x],1,function(l) paste(paste(names(l),l,sep="="),collapse="&")))
  # }

  #print(df)
  zdf <- stress.expression.addt(obj,df)
  acc <- exp(as.vector(as.matrix(zdf) %*% obj$coeff[-(1:2),drop=FALSE]))
  if(to.predict %in% c("acceleration","degradation")) acc <- acc*obj$coeff[2] 
  if(to.predict=="degradation") {
    res <- obj$coeff[1]+acc*df[[obj$varnames$t]]
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

coef.addtTwoPass <- function(obj,type=c("default","acceleration","af")) {
  type <- match.arg(type)
  switch(type,af=,acceleration=as.vector(exp(as.matrix(obj$zmodel[obj$rank$start,]) %*% obj$coeff[-(1:2)])),obj$coeff)
}

residuals.addtTwoPass <- function(obj,step=2,coeff) {#residuals after the <step>th regression
  if(!missing(coeff)) {
    resid <- obj$model[[1]]
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      resid[sub] <- obj$model.1[[1]][sub]-(coeff[1]+ coeff[1+i]*obj$model.1[[2]][sub]) 
    }
    resid
  } else if(step==1) {
    if(obj$method==1) {
      resid <- obj$model[[1]]
      for(i in seq(obj$rank$start)) {
        sub <- obj$rank$start[i]:obj$rank$end[i]
        resid[sub] <- obj$model.1[[1]][sub]-(mean(obj$intercept.1)+obj$slope.1[i]*obj$model.1[[2]][sub]) 
      }
      resid
    } else residuals(obj$lm.1)
  } else {
    resid <- obj$model[[1]]
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      resid[sub] <- obj$model.1[[1]][sub]-(obj$coeff[1]+ obj$slope.2[i]*obj$model.1[[2]][sub]) 
    }
    resid
  }
}

ftq.addtTwoPass <- function(obj,threshold,p=.95,...,plot=FALSE) { #... is for the 
  df <- data.frame(...)[obj$varnames$x]
  sigma <- sqrt(sum(residuals(obj)^2)/(nrow(obj$model.1)-1))
  ftq <- (obj$transf[[1]](threshold)-obj$coeff[1]-qnorm(if(obj$sign<0) 1-p else p)*sigma)/predict(obj,"accel",data=df)
  names(ftq) <- NULL
  if(plot && length(df)==1) {
    plot(df[[1]],ftq,type="l",main=paste("failure time quantile (threshold=",threshold,",p=",p*100,"%)",sep=""),xlab=obj$varnames$x,ylab=obj$varnames$t)
    invisible(ftq)
  } else ftq

}

summary.addtTwoPass <- function(obj,...) {
  list(coefficients=coef(obj),r.squared=1-var(resid(obj))/var(obj$model.1[[1]]),slope.2=obj$slope.2,r.squared.2=summary(obj$lm.2)$r.sq)
}

# inverse <- function (f, lower = -100, upper = 100) {
#    function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
# }


plot.addtTwoPass <- function(obj,type="degradations",with.layout=TRUE,fit=TRUE,step=2,only=NA,...) {

  if(is.numeric(type)) type <- switch(type,"degradation","g-degradation","degradations","stress plot","time resid","stress resid","normal","resid","points cloud")
  if(type=="points cloud") {
    type <- "degradation"
    step <- 0
  }
  
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

  if(step==1) lty.factor <- 1 else lty.factor <- 2

  if(!is.na(only)) {
    args$col[setdiff(seq(xx.uniq),only)] <- 0
  }  

  if(type=="degradation") {
    plot(obj$model[[2]],obj$model[[1]],xlab=obj$varnames$t,ylab=obj$varnames$y,main="degradation vs time")
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      points(obj$model[sub,2],obj$model[sub,1],col=args$col[i],pch=args$pch[i],lwd=args$lwd[i])
      if(fit) {
        if(step %in% c(1,1.2,12,2.1,21)) curve(obj$transf[[2]](mean(obj$intercept.1)+obj$slope.1[i]*x),col=args$col[i],lty=args$lty[i]*lty.factor,lwd=args$lwd[i],add=TRUE)
        if(step %in% c(2,1.2,12,2.1,21)) curve(obj$transf[[2]](mean(obj$intercept.1)+obj$slope.2[i]*x),col=args$col[i],lty=args$lty[i],lwd=args$lwd[i],add=TRUE)
      }
    }
    do.legend(args,pos="bottomleft",cex=.8)
  } else if(type=="g-degradation") {
    plot(obj$model.1[[2]],obj$model.1[[1]],xlab=obj$varnames$t,ylab=paste("g(",obj$varnames$y,")",sep=""),,main="g-degradation vs time" )
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      points(obj$model.1[sub,2],obj$model.1[sub,1],col=args$col[i],pch=args$pch[i],lwd=args$lwd[i])
      if(step %in% c(1,1.2,12,2.1,21)) abline(a=mean(obj$intercept.1),b=obj$slope.1[i],col=args$col[i],lty=args$lty[i]*lty.factor,lwd=args$lwd[i])
      if(step %in% c(2,1.2,12,2.1,21)) abline(a=mean(obj$intercept.1),b=obj$slope.2[i],col=args$col[i],lty=args$lty[i],lwd=args$lwd[i])
    }
    do.legend(args,pos="bottomleft",cex=.8)
  } else if(type=="degradations") {
    layout(matrix(c(1,2), 1,2, byrow = TRUE))
    plot(obj,"degradation",step=step,...)
    plot(obj,"g-degradation",step=step,...)
    layout(1)
  } else if(type=="stress plot") {
    nx <- length(obj$varnames$x)
    if(with.layout) layout(matrix(1:nx,1,nx, byrow = TRUE))
    if(with.layout) for(i in 1:(ncol(obj$model.2)-1)) {
      plot(obj$model.2[[1+i]],obj$model.2$y2,xlab=colnames(obj$model.2)[1+i],ylab="log(slope)",main="stress plot")
      #print(coef(obj$lm.2))
      #abline(a=coef(obj$lm.2)[1],b=coef(obj$lm.2)[1+i])
      abline(a=log(abs(coef(obj)[2])),b=coef(obj)[2+i])
    }
    if(with.layout) layout(1)
  } else if(type=="time resid") {
    resid <- residuals(obj)
    resid <- resid/mean(resid^2)
    plot(obj$model[[2]],resid,xlab=obj$varnames$t,main="residuals vs time")
    for(i in seq(obj$rank$start)) {
      sub <- obj$rank$start[i]:obj$rank$end[i]
      points(obj$model[sub,2],resid[sub],col=args$col[i],pch=args$pch[i],lwd=args$lwd[i])
      abline(h=0,lwd=3)
    }
    do.legend(args,pos="top",cex=.8)
  } else if(type=="stress resid") {
    resid <- residuals(obj)
    resid <- resid/mean(resid^2)
    nx <- length(obj$varnames$x)
    if(with.layout) layout(matrix(c(rep(1,nx),2:(nx+1)), nx,2, byrow = FALSE))
    plot(obj$model.1$xx,resid,main="residuals vs stress")
    if(with.layout) for(v in obj$varnames$x) {
      plot(obj$model[[v]],resid,xlab=v,main=paste("residuals vs",v))
    }
    if(with.layout) layout(1)
  } else if(type=="normal") {
    resid <- residuals(obj)
    qqnorm(resid)
  } else if(type=="resid") {
    layout(matrix(c(1,1,2,3), 2,2, byrow = TRUE))
    plot(obj,"time resid",...)
    plot(obj,"stress resid",with.layout=FALSE)
    plot(obj,"normal")
    layout(1)
  } 

} 


