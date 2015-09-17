\name{calcdens}
\alias{calcdens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimates the transition density.
}
\description{
  The system of ODEs is solved subject to initial conditions and the
  estimates of the mean, the variance, the macroscopic equations and
  the transition density are returned.
}
\usage{
 calcdens(initdata, edata=NA, tstart=0, tend,
   initode=NA,initmean=rep(0,length(initdata)),
   initvar=rep(0,length(initdata)*(length(initdata)+1)/2),
   thetas,relerr=1e-8, abserr=1e-8, syssize,dfunction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{initdata}{
    A numerical vector indicating the initial point. It is unscaled,
    e.g. expressed as number of molecules.}
    \item{edata}{
      Optional, a numerical vector indicating the ending point. It also
      is unscaled, e.g. expressed as number of molecules.}
    \item{tstart}{
      The starting time, defaults to 0.}
    \item{tend}{
      Either a vector or a scalar with the time-points to be estimated.}
    \item{initode}{
      Optional, the initial values of the macroscopic ODEs, defaults to
      the scaled \code{initdata}, e.g. the concentration of the species.
    }
    \item{initmean}{
      A numerical vector indicating the initial values for the means.
      Defaults to a vector of zeroes, otherwise it is expected to be
      scaled by the inverse of the square root of the system size.}
    \item{initvar}{
      Either a matrix indicating the initial Variance-Covariance matrix or
      a vector representing the upper diagonal (including the main
      diagonal) following a row orientation. Defaults to a matrix of
      zeroes and is expected to be on the scale of macroscopic ODEs.}
    \item{thetas}{
      A numerical vector with the parameter values.
    }
    \item{relerr}{
      Numerical, the relative error for the numerical ordinary
      differential equations (ODEs) solver.
    }
    \item{abserr}{
      Numerical, the absolute error for the numerical ordinary
      differential equations (ODEs) solver.
    }
    \item{logprob}{
      Boolean, indicate if the log of the probability should be returned.
    }
    \item{syssize}{
      Numerical, indicating the system size.}
    \item{dfunction}{
      The compiled function, given as a loaded dynamic library
      in R or as a character string of the symbol's name (similar to the
      \code{\link[base]{is.loaded}} argument).
    }
  }
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list of the following components, estimated at each \code{tend}
  time-point:
  \item{Time}{The time instance of the estimates.}
  \item{ODE}{The value of the ODE equation (the macroscopic model),
    expressed in concentration.}
  \item{MEAN}{The mean of the SDE process, expressed in mesoscopic units
    (multiply by sqrt(system size) to convert to original units).}
  \item{VAR}{The covariance of the SDE process (multiply by system size to
      convert to original units).}
  \item{prob}{Optional, expresses the estimated transition probability
    density, available only if\code{edata} is defined.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Vasileios Giagos
}
\note{
  All densities are conditioned on the initial time-point
  \code{tstart}.  The \code{MEAN} and \code{VAR} elements are not
  at the same scale but they depend on the scale of the initial
  values. We assume that the initial values are given as number of
  molecules.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  The model parsing is described in \code{\link[lnar]{parsemod}}.
}
\examples{
\dontrun{
require(lnar)
##We consider the Lotka-Volterra Model
tt <- matrix(c(1,-1,0,0,1,-1),nrow=2,ncol=3,byrow=TRUE)
rfun <- c("con1 * Prey","con2 * Prey * Predator","con3 * Predator")
thetas <- paste("con",1:3,sep="")
species <- c("Prey","Predator")
cout <- parsemod(tt,rfun,thetas,species) #Parse the model

##Inputs a dataset
initdata<-c(50.0, 30)
data2<-c(51, 28)

compmod(cout,"derivs") #Compile the model

##Test that derivs is working.
derivs(1,c(initdata[1],initdata[2],c(0,0,0,0,0)),
       rep(0,7),c(.1,.0001,.1)) 

syssize=sum(initdata)
thetas <- c(0.25,0.20,0.125) #scaled kinetic constants
##Calculate Transition Density
(calc1<-calcdens(initdata,data2,tstar=0,tend=.1,
                thetas=thetas,
                syssize=syssize,
                dfunction=derivs))
##Test:
log(calc1[[1]]$prob) # -4.835931

##Calculate the trans. dens. parameters in some time points

(calc2<-calcdens(initdata,tstar=0,tend=c(.1,.5,3),
                thetas=thetas,
                syssize=syssize,
                dfunction=derivs))
}
}

