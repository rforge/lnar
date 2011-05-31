library(lnar)
##Number of Species: 4
species=c('RNA','P','P2','DNA')
#Number of Parameters: 8
params=c('k1','r1','k2','k3','k4','r4','k5','k6')
stoich=matrix(c( 0,0,1,0,0,0,-1,0,  
                 0,0,0,1,-2,2,0,-1,
                 -1,1,0,0,1,-1,0,0,
                 -1,1,0,0,0,0,0,0
  ),4,8,byrow=TRUE)
#Number of Reactions: 8
 reac=c(
    'k1*DNA*P2',
    'r1*(0.2941176-DNA)',
    'k2*DNA',
    'k3*RNA',
    'k4*0.5*P*P',
    'r4*P2',
    'k5*RNA',
    'k6*P'
   )
##generate the model and c code
model1<-parsemod(stoich,reac,params,species)
compmod(model1,"tder") #compile model, loads it as "tder"

##load the data
data(ardata)

##We set all c's to .2 for our initial values
nthetas <- rep(.2,8)
nthetas[1]=nthetas[1]*34 # corresponds to a 2nd order reaction
nthetas[5]=nthetas[5]*34 # corresponds to a 2nd order reaction

##Optimize with Nelder-Mead
(model1opt<-optmod(model1,nthetas=nthetas, mydata=ardata, method=0,
              maxiter=1800, tcrit=1e-5, relerr=1e-12,
              abserr=1e-12, hessianh=1e-4,syssize=34,
              dfunction=tder))
##Continue with BFGS
(model2opt<-optmod(model1,nthetas=model1opt$ES, mydata=ardata, method=0,
              maxiter=25, tcrit=1, relerr=1e-12,
              abserr=1e-12, hessianh=1e-4,syssize=34,usebfgs=1,
              dfunction=tder))
##Calculate the transition density's parameters at t=1
calcdens(as.numeric(ardata[1, -1]),tend=1,
         thetas=model2opt$ES,syssize=34,dfunction=tder)

##Evaluate the log-likelihood at the mles

c.par <- model2opt$ES
c.par[model1$Order>1] <- c.par[model1$Order>1]/34
(l2<-lnalik(model1,nthetas=c.par, mydata=ardata, method=0,
                   relerr=1e-9, abserr=1e-9,
                   dfunction=tder,syssize=34) )

c.par2 <- model2opt$ES
c.par2[model1$Order>1] <- c.par2[model1$Order>1]*34
(l1<-lnalik(model1,nthetas=model2opt$ES, mydata=ardata, method=1,
                   relerr=1e-9, abserr=1e-9,syssize=34,
                   dfunction=tder) )

(l1<-lnalik(model1,nthetas=c.par2, mydata=ardata, method=1,
                   relerr=1e-9, abserr=1e-9,syssize=34,
                   dfunction=tder) )
