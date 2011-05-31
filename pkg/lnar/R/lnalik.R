lnalik <- function(cout,nthetas,
                     mydata,
                     syssize=NA,
                     relerr=1e-9,
                     abserr=1e-9,
                     method=0,
                     dfunction
                     )
  {
    if(!is.numeric(nthetas))
      stop("nthetas must be numeric.")
    if(is.data.frame(mydata)) mydata <- as.matrix(mydata)
    if(!is.matrix(mydata) )
       stop("mydata must be a matrix or data frame.")
    if(is.na(syssize)) 
    {
      syssize<-sum(mydata[1,-1])
      cat("Using the initial state to calculate the system's size: ",
          syssize,"\n")
    }
    if(!is.numeric(syssize))
       stop("syssize must be numeric.")

    if(!is.numeric(relerr))   stop("relerr must be numeric.")
    if(!is.numeric(abserr))   stop("abserr must be numeric.")
    if(!is.numeric(method))   stop("method must be numeric.")
    if(!isGeneric(dfunction@generic[1]))
       stop("dfunction must be a \"generic\" compiled method.")
    ans <- with(cout,{
      dimdata <- dim(mydata)
      odesize<-length(cspecies)+length(Cov)+length(Means)
      mydata <- as.vector(t(mydata))
      dfunction(as.integer(odesize), as.double(mydata[1]), 
            as.double(mydata[2:(1 + length(cspecies))]), 
            as.double(rep(0, odesize)), as.double(nthetas))
      if(is.na(match(method,c(0,1,3)))) stop("Please specify a valid method.")
      storage.mode(mydata)<-"double"
      .Call("calclik",
            ptrf=getNativeSymbolInfo(dfunction@generic[1])$address,
            NPARAMS = as.integer(length(nthetas)), 
            NODES=as.integer(odesize),
            nmols=as.integer(syssize),
            relerr=as.double(relerr),
            abserr=as.double(abserr),
            dataset= mydata, #format by row
            nrows=as.integer(dimdata[1]),
            ncols=as.integer(dimdata[2]),
            ord=as.integer(Orders),
            thetas=as.double(nthetas), 
            method= as.integer(method),PACKAGE="lnar")
    })
    #Unload the compile code
    #dyn.unload(
    #  unlist(getNativeSymbolInfo(dfunction@generic[1])$package["path"]))

    return(ans)
  }
