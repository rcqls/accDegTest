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

	class(obj) <- "addt0"

	# varnames 
	obj$varnames<-list(
		y=all.vars(obj$formula[[2]]),
		t=all.vars(obj$formula[[length(obj$formula)]][[2]]),
    	x=all.vars(obj$formula[[length(obj$formula)]][[3]])
	)
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
	
	# 
	prepare.estim.addt0(obj) 
	
	obj
}


prepare.estim.addt0 <- function(obj) {
	# system are uniquely determined with cond and system name 
	syst <- paste(obj$data[[obj$varnames$system]],apply(obj$data[obj$varnames$x],1,paste,collapse="&"),sep="&")
	measure <- obj$data[[obj$varnames$measure]]
	# measure 0
	y0 <- obj$transf[[1]](obj$data[measure==0,obj$varnames$y])
	#print(y0)

	# measure 1
	data1 <- cbind(obj$data[measure==1,], data.frame(dy=sapply(unique(syst),function(s) {
			obj$transf[[1]](obj$data[syst == s & measure == 1, obj$varnames$y])
			- obj$transf[[1]](obj$data[syst == s & measure == 0, obj$varnames$y])
		}
	)))

	#print(dim(data1))
	form <- update(eval(as.call(c(obj$formula[[1]],obj$formula[[2]],obj$formula[[3]][[2]]))),.~0+.)
	formula <- obj$formula
	formula[[3]][[2]] <- form[[3]]
	obj$addt <- addt(formula,obj$xref,obj$transf,data=data1)

}