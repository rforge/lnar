Mult <- function(A,B){
    dA=dim(A)
    inspar <- function(x) paste0("(",x,")")
    if(is.vector(B)){        
        stopifnot(dA[2]==length(B))
        res <- character(dA[1])
        for(i in 1:dA[1]){
            for(j in 1:dA[2]){
                tmp=paste(inspar(A[i,j]),inspar(B[j]),sep = "*")
                tmp=Deriv::Simplify(tmp)
                res[i] <- ifelse(j==1,tmp,paste(res[i],tmp,sep = "+"))
            }
            res[i]=Deriv::Simplify(res[i])
        }
    }else{
        dB=dim(B)
        stopifnot(dA[2]==dB[1])
        res <- matrix(" ",
                      nrow = dA[1], ncol = dB[2])
        for(i in seq.int(dA[1]))
            for(j in seq.int(dB[2])){
                for(k in seq.int(dB[1])){
                    tmp=paste(inspar(A[i,k]),inspar(B[k,j]),sep="*")
                    tmp=Deriv::Simplify(tmp)
                    res[i,j] <- ifelse(k==1,tmp,paste(res[i,j],tmp,sep="+"))                    
                }
                res[i,j]=Deriv::Simplify(res[i,j])
            }                
    }
    res
}
Add <- function(A,B){
    dA <- dim(A)
    dB <- dim(B)
    stopifnot(dA[1]==dB[1] && dA[2]==dB[2])
    res <- matrix(NA,dA[1],dA[2])
    res[] <- mapply(paste,A,B,sep="+")
    res[] <- vapply(res,Deriv::Simplify,as.character(0))
    res
}
eqpars <- function(x,subs){
    expr1 = parse(text=x)
    envsubs = sapply(subs,as.name)
    tmp1=lapply(expr1, function(x)
        do.call(what=substitute, list(expr=x,env=envsubs)))
    vapply(X=tmp1,FUN=deparse,FUN.VALUE = as.character(0),width.cutoff=500,backtick = FALSE)
}

parsemod<-function(stoich,rfun,thetas,species,constants=NA)
{
    ## C language conversions
    cexpr <- c("pow")
    names(cexpr) <- c("^")
    ######## Set up Variables/Names
    ## Thetas
    cthetas = sapply(seq_len(length(thetas))-1,function(x) sprintf("vthetas[%i]",x))
    names(cthetas) = thetas
    nspecies = length(species)
    ncov = nspecies*(nspecies +1)/2 
    odesize = 2*nspecies + ncov
    # ODE Species
    yode = sapply(0:(nspecies-1),function(x) sprintf("y[%i]",x))
    names(yode) = species    
    ## Covariance
    ## Here we have different names => symbols 
    covnames1<-outer(species,species,FUN="paste",sep=",")
    ## Default: Cov(x1,x2) => y[2], human friendly-explicit
    ccov = sapply(seq(ncov)+nspecies-1,function(x) sprintf("y[%i]",x))
    names(ccov) <- paste("Cov(",covnames1[upper.tri(covnames1,TRUE)],")",sep="")
    ccovs <- sapply(seq(ncov),function(x) sprintf("s%i",x))
    ## Math notation, s[1, 1] => "y[2]" -- Not used
    if(F){
        ccov2=ccov # S-names with y notation
        NamesMat <- outer(paste0("s[",seq(nspecies)),paste0(seq(nspecies),"]"),paste,sep=",")
        NamesVec <- NamesMat[upper.tri(NamesMat,T)]
        covnames <- vapply(
            NamesVec, Deriv::Simplify,
            as.character(0), USE.NAMES = F)
        names(ccov2) <- covnames
    }
    ## s1 => y[2] used for symbolic calcs
    ccov3=ccov
    names(ccov3) = ccovs 
    Smat = outer(paste0("s[",seq(nspecies)),paste0(seq(nspecies),"]"),paste,sep=",")
    ##Smat = matrix(paste0("s",seq(nspecies^2)),nspecies,nspecies,byrow = T)
    ## Crate Smat
    Smat[lower.tri(Smat,T)] = ccovs 
    Smat[upper.tri(Smat)]=t(Smat)[upper.tri(Smat)]
    ## Residual Process
    mt = sapply(seq(nspecies)-1 + ncov +nspecies,function(x) sprintf("y[%i]",x))
    names(mt) = paste0("m",seq(nspecies))
    c(yode,ccov3,mt)
    ##-------- Start the symbolic calcs
    ## dim(A')=(Nspecies,Nreactions), stoich = A'
    ## dy = A' h(y)
    (dY=Mult(stoich,rfun))
    subs=c(species,cthetas,yode,ccov3,mt,cexpr)
    ## Generate Fmat = A'*F'  where: F[i,j] = d(h_j(y)) / dm_i
    lenrf = length(rfun)
    Fmat=matrix(NA,nspecies,lenrf,byrow = T)
    for(i in seq(nspecies)){
        for(j in seq(lenrf)){
            Fmat[i,j] = Deriv::Deriv(rfun[j],species[i])
        }
    }
    Fmat=Mult(stoich,t(Fmat))
    ## dS = F S + S F' + A' diag(h(y)) A
    s1=Mult(Fmat,Smat)
    s2=Mult(Smat,t(Fmat))
    s3 = Add(s1,s2)
    Mh <- diag(lenrf)
    diag(Mh) <- rfun
    s4=Mult(stoich,Mult(Mh,t(stoich)))
    s5 = Add(s4,s3)
    dS <-  s5[lower.tri(s5,T)]
    ## dm = F m
    m = paste0("m",seq(nspecies))
    dm = Mult(Fmat,m)
    ## Output C code
    eqs = c(eqpars(dY,subs),eqpars(dS,subs),eqpars(dm,subs))
    cout=vapply(seq(odesize),function(x) paste0("fout[",x-1,"]= ",eqs[x],";"),as.character(0))
    cout=paste(cout,collapse = "\n")
    yrfun <- eqpars(rfun,subs)
    names(mt) <- paste("Mean",species)
    return(list(ccode=cout,
                cspecies=yode,
                cthetas= cthetas,
                Cov=ccov,
                Means=mt,
                ##Finds the reaction orders
                Orders= ( sapply(strsplit(yrfun,"y\\["),length)-1 )
                )
           )
}
