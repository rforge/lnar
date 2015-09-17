calcdens <- function(initdata,
                   edata=NA,
                   tstart=0,
                   tend,
                   initode=NA,
                   initmean=rep(0,length(initdata)),
                   initvar=rep(0,length(initdata)*(length(initdata)+1)/2),
                   thetas,                   
                   relerr=1e-8,
                   abserr=1e-8,
                   syssize, 
                   logprob=FALSE,
                   dfunction
                   )
  {
    if(!is.numeric(thetas))
      stop("Thetas must be specified as numeric.")
    if(!is.numeric(initdata))
       stop("initdata must be numeric.")
    if(!is.numeric(syssize))
       stop("syssize must be numeric.")
    if(!is.numeric(tstart))
       stop("tstart must be numeric.")
    if(!is.numeric(tend))
       stop("tend must be numeric.")
    
    if(!is.numeric(relerr))   stop("relerr must be numeric.")
    if(!is.numeric(abserr))   stop("abserr must be numeric.")

        if(is.na(initode)) initode <- initdata/syssize
    
    if(is.matrix(initvar)){
     
      if(!isSymmetric(initvar)) stop("The covariance matrix is not symmetric.")
      ecomp <- eigen(initvar, symmetric = TRUE)
      if(!all(ecomp$values>0))
        stop("The covariance matrix is not positive-definite")
      ## Return the upper diagonal by row.
      initvar <- t(initvar)[rev(upper.tri(initvar,diag=TRUE))]
    }
    
    odesize<-length(initdata)+length(initvar)+length(initmean)

    #Test the dfunction argument
    if(is.character(dfunction))
       {
         if(!class(getNativeSymbolInfo(dfunction)$address) == "NativeSymbol")
           {
             stop("The character string supplied in dfunction does not correspond to a symbol name.")
           }
         else
           {
             ptr = getNativeSymbolInfo(dfunction)$address
           }
      }   else
      {
        if(isGeneric(dfunction@generic[1]))
          {
            dfunction(odesize,tstart,c(initdata,initvar,initmean),rep(0,odesize),thetas)
            ptr = getNativeSymbolInfo(dfunction@generic[1])$address            
          }
        else{
          stop("dfunction must be either a generic compiled method or a character string giving a symbol name.")
          }
      }
    

    ans<-.Call("calcdens",
          initdata= as.double(initode),
          tstart= as.double(tstart),
          tend= as.double(tend),
          initmean= as.double(initmean),
          initvar= as.double(initvar),
          thetas=as.double(thetas),
          sodesize=as.integer(odesize),
          relerr=as.double(relerr),
          abserr=as.double(abserr),            
          ptrf=ptr,
          PACKAGE="lnar")

    if(any(is.na(edata))) {
      #ans$prob <- NA
    } else
    {
      if(is.matrix(edata)){
        lendat<- dim(edata)[1]
        if(lendat==length(tend)){
          for(i in 1:lendat){ 
            if(!all(is.numeric(edata[i,])))
              stop("data must be numeric.")
            if(sum(initmean)==0)
              {
                x <- edata[i,] - ans[[i]]$ODE*syssize
                y <- ans[[i]]$VAR*syssize
              } else
            {
              x <- edata[i,] - ans[[i]]$ODE*syssize -
                ans[[i]]$MEAN * sqrt(syssize)
              y <- ans[[i]]$VAR*syssize
            }
            ans[[i]]$prob <- dmnorm(x,y)
          }
        }        
      } else
      {
        if(!all(is.numeric(edata)))
          stop("data must be numeric.")
        if(sum(initmean)==0)
          {
            x <- edata - ans[[1]]$ODE*syssize
            y <- ans[[1]]$VAR*syssize
          } else
        {
          x <- edata- ans[[1]]$ODE*syssize -
            ans[[1]]$MEAN * sqrt(syssize)
          y <- ans[[1]]$VAR*syssize
        }
        ans[[1]]$prob <- dmvnorm(x,sigma=y,log=logprob) 
      }
    }
    return(ans)
    ##dyn.unload(unlist(getNativeSymbolInfo("derivs")$package["path"]))
  }
