\name{optmod}
\alias{optmod}
\title{
MLEs for a given model.
}
\description{
  Fits the compiled model to a given datased using a numerical Maximul
  Likelihood Estimation procedure.
}
\usage{
optmod(cout,nthetas, mydata, maxiter=300,
       syssize=sum(mydata[1,-1]), tcrit=.0001,
       relerr=1e-9, abserr=1e-9, hessianh=1e-4,
       method=1, usebfgs=0, dfunction)
}
\arguments{
  \item{cout}{
    The parsed model containing the C code and the name
    relations.
  }
  \item{nthetas}{
    A numerical vector with the initial values for the
    parameters to be optimized.
  }
  \item{mydata}{
    A \code{data.frame} or a matrix with the data.
    The first column must indicate the \emph{time} of the observations.
  }
  \item{maxiter}{
    Numerical, indicated the maximum number of iterations for the
    optimization algorithm. 
}
  \item{syssize}{
    Numerical, the system size defaults to the initial population.
}
  \item{tcrit}{
    Numerical, the convergence criterion for the optimization
    algorithm.
}
  \item{relerr}{
    Numerical, the relative error for the numerical ordinary
    differential equations (ODEs) solver.
}
  \item{abserr}{
    Numerical, the absolute error for the numerical ordinary
    differential equations (ODEs) solver.
}
  \item{hessianh}{
    Numerical, indicates the approximation step for the central
    differences calculations of the Hessian matrix.
}
  \item{method}{
    Numerical which takes the following integer values:
    \itemize{
      \item{0: Restarting Method using number of molecules.}
      \item{1: Restarting Method using concentrations.}
      \item{3: Non-Restarting Method using concentrations.}
    }
    a discussion on the Restarting and the Non Restarting method can be
    found in Giagos (2010). 
} 

  \item{usebfgs}{
    Specify whether to use the BFGS algorithm (1), or the default
    Nelder-Mead    simplex algorithm (0)
  }
  \item{dfunction}{
    The compiled function, given as a loaded dynamic library in R.
}
}
\details{
  By default the \code{L-BFGS-B} optimization procedure is employed
  (see \code{\link[stats]{optim}} for more details).
}
\value{
  A list with the following elements:
  \item{UP}{The upper confidence bound.}
  \item{ES}{The MLEs}
  \item{LO}{The lower confidence bound.}
}
\references{
  Giagos, V.: 2010, \emph{Inference for auto-regulatory genetic networks
  using diffusion   process approximations}, Thesis, Lancaster University,
  2010.
}
\author{
Vasileios Giagos
}
\note{
  Note that the \code{tcrit} has a different interpretation for the
  Nelder-Mead algorithm and different for BFGS
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  The model parsing is described in \code{\link[lnar]{parsemod}}.
}
\examples{

\dontrun{
##We consider the Lotka-Volterra Model
tt <- matrix(c(1,-1,0,0,1,-1),nrow=2,ncol=3,byrow=TRUE)
rfun <- c("con1 * Prey","con2 * Prey * Predator","con3 * Predator")
thetas <- paste("con",1:3,sep="")
species <- c("Prey","Predator")
cout <- parsemod(tt,rfun,thetas,species) #Parse the model

##Inputs a dataset
mydata<-c(0.0, 5000.0, 3000, 1, 5989, 2992, 2, 7165, 3107, 3, 8534, 3306, 
	4, 10041, 3709, 5, 11624, 4265, 6, 13306, 5181, 7, 14741, 6492, 
	8, 15867, 8337, 9, 16025, 10981)
mydata2<-t(matrix(mydata,3,10)) # Put the data into a matrix

compmod(cout,"derivs") #Compile the model

##Test that derivs is working.
derivs(mydata[1],c(mydata[2],mydata[3],c(0,0,0,0,0)),
       rep(0,7),c(.1,.0001,.1)) 

##Run model
run1<-optmod(cout,nthetas=c(.5,.1,.2), mydata=mydata2,
             maxiter=300, simplexsize=1e-5, relerr=1e-9,
             abserr=1e-9, hessianh=1e-4,dfunction=derivs)

##Try a random starting value. WARNING! May crash!
(nthetas<-exp(runif(3,-3,1)))
run2<-optmod(cout,nthetas=nthetas, mydata=mydata2,
             maxiter=500, simplexsize=1e-6, relerr=1e-9,
             abserr=1e-9, hessianh=1e-4,dfunction=derivs)
}
}

