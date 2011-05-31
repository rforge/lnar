require(lnar)
tt <- matrix(c(1,-1,0,0,1,-1),nrow=2,ncol=3,byrow=TRUE)
rfun <- c("con1 * Prey","con2 * Prey * Predator","con3 * Predator")
thetas <- paste("con",1:3,sep="")
species <- c("Prey","Predator")
cout <- parsemod(tt,rfun,thetas,species)
mydata<-c(0.0, 5000.0, 3000, 1, 5989, 2992, 2, 7165, 3107, 3, 8534, 3306, 
	4, 10041, 3709, 5, 11624, 4265, 6, 13306, 5181, 7, 14741, 6492, 
	8, 15867, 8337, 9, 16025, 10981)
mydata2 <- matrix(mydata,10,3,byrow=TRUE) # Example dataset

compmod(cout,"derivs") # Compile the model
nthetas<-c(.4,.1,0.4)  # The initial parameter values

#Find the Maximum Likelihood Estimates and Wald CIs
(run1<-optmod(cout,nthetas=nthetas, mydata=mydata2, method=1,
              maxiter=300, tcrit=1e-5, relerr=1e-9,
              abserr=1e-9, hessianh=1e-4,
              dfunction=derivs))

##Calculate the transition density's parameters at t=1
calcdens(mydata2[1,-1],tend=1,thetas=run1$ES,syssize=80,dfunction=derivs)

##Evaluate the log-likelihood at the mles
(l1<-lnalik(cout,nthetas=run1$ES, mydata=mydata2, method=1,
                   relerr=1e-9, abserr=1e-9,
                   dfunction=derivs) )
