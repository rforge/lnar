parsemod<-function(y,rfun,thetas,species,
                     constants=NA)
{
    lenthetas <- length(thetas)
    lenspecies <- length(species)
    
    dimy<- dim(y)
    txt<-paste(sapply(seq(dimy[1]),
	    function(x) paste("{",paste(y[x,],sep="",collapse=","),"}") ),
	    collapse=",",sep=",")
    
    stoich<-paste("A:={",txt,"}",sep="") #Stoichiometry Matrix
    yspecies <- paste("y[",seq(lenspecies),"]",sep="")
    names(yspecies) <- species
    ythetas <- paste("theta[",seq(lenthetas),"]",sep="")
    names(ythetas) <- thetas
    yrfun <- rfun

    makesubs <- function(orig,repl,string)
      {
        ##Regular Expressions between math symbols
        regpre <- "([[:punct:][:space:]]+)"
        regpost <- "([[:punct:][:space:]]+)"
        string <- gsub(paste(regpre,orig,regpost,sep=""),
                    paste("\\1",repl,"\\2",sep=""),
                    string)
        string <- gsub(paste("^",orig,regpost,sep=""),
                       paste(repl,"\\1",sep=""),
                       string)
        string <- gsub(paste(regpre,orig,"$",sep=""),
                    paste("\\1",repl,sep=""),
                    string)
	string <- gsub(paste(regpre,orig,regpost,sep=""),
                    paste("\\1",repl,"\\2",sep=""),
                    string)

        return(string)
      }

    for (x in seq(lenspecies)) {
      yrfun <- makesubs(species[x],yspecies[x],yrfun)
    }

    for (x in seq(lenthetas)) {
      yrfun <- makesubs(thetas[x],ythetas[x],yrfun)
    }
    if(!is.na(constants))
      {
        #Substitute the constants in the equation
        lencon <- length(constants)
        for (x in seq(lencon)) {
          yrfun <- makesubs(names(constants)[x],constants[x],yrfun)
        }
      }    
    yrates <- paste("rates:={",paste(yrfun,collapse=","),"}")

    #paste(find.package("lnar"),"/extra/test.ys",sep="") #yacas template
    require(Ryacas)
    #options(yacas.method = 'system')
    yacfile <- paste(find.package("lnar"),"/extra/test.ys",sep="")

    if(package_version(packageDescription("Ryacas",fields="Version"))
       <"0.2.11") {
      yacas(stoich)
      yacas(yrates)
      cout <- yacas(paste("Load(\"",yacfile,"\")",sep=""))
      yacasStop()
    }else{
      cout <- yacas(
                    paste("[",
                          stoich,";",
                          yrates,";",
                          paste("Load(\"",yacfile,"\")",sep=""),";]",
                          sep=""),
                    method="system")
    }
    #print(cout)
    cout$PrettyForm <- grep("^[fd]",cout$PrettyForm,value=TRUE)
    cout <- paste( paste(cout$PrettyForm,collapse="\n") ,collapse="")
    cout<-sub("True\n","",as.character(cout))
    cspecies<- paste("y[",seq(lenspecies)-1,"]",sep="")
    cthetas <- paste("vthetas[",seq(lenthetas)-1,"]",sep="")
    names(cspecies) <- species
    names(cthetas)<-thetas
    lens <- (lenspecies+1)*(lenspecies)/2
    Covariances<-paste("y[",seq(lenspecies, lenspecies+lens -1),"]",sep="")
    xx<-outer(species,species,FUN="paste",sep=",")
    names(Covariances) <- paste("Cov(",xx[upper.tri(xx,TRUE)],")",sep="")
    mspecies <- paste("y[",seq(lens+lenspecies, lens+2*lenspecies -1),"]",
                      sep="")
    names(mspecies) <- paste("Mean",species)
    return(list(ccode=cout,
                cspecies=cspecies,
                cthetas= cthetas,
                Cov=Covariances,
                Means=mspecies,
                Orders= ( sapply(strsplit(yrfun,"y\\["),length)-1 ),
                yac=paste("[",
                          stoich,";",
                          yrates,";",
                          paste("Load(\"",yacfile,"\")",sep=""),";]",
                          sep="")
		##Finds the reaction orders
                )
           )
}
