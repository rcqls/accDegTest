
### addt simulator
sim.addt <-function(formula,param,noise.level=.1,xref=1,transf=c(log,exp),data) {
  if(length(transf)==1) {
    if(identical(transf,log)) transf <- c(log,exp)
    else if(identical(transf,exp)) transf <- c(exp,log)
    else stop("transf argument generally needs to be of length 2.")
  }
  obj<-list(formula=formula,param=param,noise.level=noise.level,xref=xref,transf=transf)
  class(obj) <- "sim.addt"

  obj <- if(missing(data)) init.addt(obj) else init.addt(obj,data=data)

  obj
}

sim.addt.demo <- function(a=c(2,3),b=2,nX=6,nSys=10,...) {
  x<-rep(2:nX,if(length(nSys)==1) rep(nSys,nX-1) else nSys)
  sim.addt(
    ~ t | 1/x,
    param=c(a,b),
    data=data.frame(x=x,t=runif(length(x),1,4)),
    ...
  )
}

model.frame.sim.addt <- function(obj,...) model.frame.addt(obj,...)

## simulate a new realisation
simulate.sim.addt <- function(sim,nsim=1,seed=NULL,as.addt=FALSE,y) { #or "y"
  if(!missing(y)) {
    y<-as.character(substitute(y))
    as.addt <- TRUE
  } else y <- "y"
  Af <- exp(as.matrix(sim$zmodel) %*% as.matrix(sim$param[-(1:2)]))           #version de z unidimensionnelle en premier
  res <- sim$transf[[2]](sim$param[1]+sim$param[2]*Af*sim$model[[1]]+rnorm(length(Af),sd=sim$noise))     #rnorm can be changed
  if(as.addt) {
    sim$model <- cbind(y=res,sim$model)
    names(sim$model)[1] <- y
    sim$formula <- as.formula(paste(y,"~",deparse(sim$formula[[2]])))
    res <- addt(sim$formula,data=sim$model)
    #attach metadata to addt object to return
    attr(res,"param") <- sim$param              
    attr(res,"noise.level") <- sim$noise.level  
  }
  res 
}

print.sim.addt <- function(obj,...) {
  cat("Simulator of ADDT model:\n")
  cat("Formula:",as.character(formula(obj)),"\n")
  cat("parameters values:",obj$param,"\n")
}

