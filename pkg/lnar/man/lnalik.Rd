\name{lnalik}
\alias{lnalik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear Noise Approximation log-likelihood.
}
\description{
Estimates the log-likelihood of the LNA approximation.
}
\usage{
lnalik(cout, nthetas, mydata, syssize = NA,
       relerr = 1e-09, abserr = 1e-09, method = 0, dfunction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cout}{
    The parsed model.
}
  \item{nthetas}{
    The vector of the parameters.
}
  \item{mydata}{
    Either a matrix or a data frame of the data to be evaluated. The
    first column is assumed to correspond to the time of each
    observation.
}
  \item{syssize}{
    Optional, a scalar indicating the system size.
}

  \item{relerr}{
    Optional, a scalar indicating the relative error for the ODE solver.
}
  \item{abserr}{
    Optional, a scalar indicating the absolute error for the ODE solver.
}
  \item{method}{
    Optional, a scalar with possible options:
    \itemize{
      \item 0: Restarting method using concentrations. The
    parameters are assumed to be scaled, i.e. \eqn{\theta}{thetas}.
      \item 1: Restarting method using number of molecules. The
    parameters are assumed to be un-scaled, i.e. \eqn{c}{c}.
      \item 3: Non-Restarting method using concentrations. The
    parameters are assumed to be scaled as well.
    }
}
  \item{dfunction}{
    The compiled model.
}
}
\details{
  See Giagos (2010) for a discussion on the Restarting and the Non
  Restarting method.  
}
\value{
  Returns the estimated log-likelihood.
}
\references{
  Giagos, V.: 2010, \emph{Inference for auto-regulatory genetic networks
  using diffusion   process approximations}, Thesis, Lancaster University,
  2010.
}
\author{
  Vasileios Giagos
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
require(lnar)
tt <- matrix(c(1,-1,0,0,1,-1),nrow=2,ncol=3,byrow=TRUE)
rfun <- c("con1 * Prey","con2 * Prey * Predator","con3 * Predator")
thetas <- paste("con",1:3,sep="")
species <- c("Prey","Predator")
cout <- parsemod(tt,rfun,thetas,species)

mydata<-c(0.0, 5000.0, 3000, 1, 5989, 2992, 2, 7165, 3107, 3, 8534,
          3306,4, 10041, 3709, 5, 11624, 4265, 6, 13306, 5181, 7,
          14741, 6492,8, 15867, 8337, 9, 16025, 10981)

mydata2 <- matrix(mydata,10,3,byrow=TRUE)#Example dataset

compmod(cout,"derivs")

#Our initial values
nthetas<-c(.25,.20,0.125)

print(derivs(mydata[1],c(mydata[2],mydata[3],
                         c(0,0,0,0,0)),rep(0,7),nthetas))

(l1<-lnalik(cout,nthetas=nthetas, mydata=mydata2, method=1,
              relerr=1e-9, abserr=1e-9,
              dfunction=derivs) )
nthetas2<-c(.25,.20/8000,0.125)
(l2<-lnalik(cout,nthetas=nthetas2, mydata=mydata2, method=0,
              relerr=1e-9, abserr=1e-9,
              dfunction=derivs) )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
