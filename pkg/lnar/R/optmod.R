optmod <- function(cout,nthetas,
                     mydata,
                     maxiter=300,
                     syssize=sum(mydata[1,-1]),
                     tcrit=.0001,
                     relerr=1e-9,
                     abserr=1e-9,
                     hessianh=1e-4,
                     method=1,
                     usebfgs=0,
                     hess=1,
                     dfunction
                     )
  {
    if(!is.numeric(nthetas))
      stop("nthetas must be numeric.")
    if(is.data.frame(mydata)) mydata <- as.matrix(mydata)
    if(!is.matrix(mydata) )
       stop("mydata must be a matrix or data frame.")
    if(!is.numeric(maxiter))
       stop("maxiter must be numeric.")
    if(!is.numeric(maxiter))
       stop("maxiter must be numeric.")
    if(!is.numeric(syssize))
       stop("syssize must be numeric.")

    if(!is.numeric(tcrit))   stop("Termination criterion must be numeric.")
    if(!is.numeric(relerr))   stop("relerr must be numeric.")
    if(!is.numeric(abserr))   stop("abserr must be numeric.")
    if(!is.numeric(hessianh))   stop("hessianh must be numeric.")
    if(!is.numeric(method))   stop("method must be numeric.")
    if(!isGeneric(dfunction@generic[1]))
       stop("dfunction must be a generic compiled method.")
  ans <- with(cout,{
  dimdata <- dim(mydata)
  odesize<-length(cspecies)+length(Cov)+length(Means)
  mydata <- as.vector(t(mydata))
  if(is.na(match(method,c(0,1,3)))) stop("Please specify a valid method.")
  if(is.na(match(usebfgs,c(0,1))))
    stop("The argument usebfgs can be either 0 or 1 .")
  if(is.na(match(hess,c(0,1))))
    stop("The argument hess can be either 0 or 1 .")
  ##By calling the compiled code, we force R to load it.
  dfunction(odesize,mydata[1],mydata[2:(1+odesize)],rep(0,odesize),nthetas)
  
  storage.mode(mydata)<-"double"
  .Call("runmodel",
              derivs=getNativeSymbolInfo(dfunction@generic[1])$address,
              NPAPAMS = as.integer(length(nthetas)), 
              NODES=as.integer(odesize),
              maxiter=as.integer(maxiter),
              nmols=as.integer(syssize),
              tcrit=as.double(tcrit), 
              hessianh= as.double(hessianh),
              relerr=as.double(relerr),
              abserr=as.double(abserr),
              dataset= as.double(mydata), #format by row
              nrows=as.integer(dimdata[1]),
              ncols=as.integer(dimdata[2]),
              ord=as.integer(Orders),
              thetas=as.double(nthetas),
              usebfgs=as.integer(usebfgs),
              method=as.integer(method),
              hess=as.integer(hess),
              PACKAGE="lnar")
   })
    #Unload the compile code
    #dyn.unload(
    #  unlist(getNativeSymbolInfo(dfunction@generic[1])$package["path"]))
    return(ans)
}
