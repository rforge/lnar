\name{parsemod}
\alias{parsemod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Parses a biological model to C code.
}
\description{
  Given as input the reaction rates, stochastic constants and model
  constants, the function outputs the C code of the underlying ODEs. The
  system of odes express the macroscopic approximation (ODE) as well as
  the estimates for the instantaneous mean and covariance of the linear
  SDE from the LNA approximation.
}
\usage{
parsemod(y,rfun,thetas,species,constants=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    The stoichiometry matrix, note that the dimensions are assumed
    to be: (# Species)\eqn{\times}{*}(# Reactions).
  }
  \item{rfun}{
    A character vector with each elements expressing the reaction
    rates.  
  }
  \item{thetas}{
    A character vector which denotes the stochastic constant names
    associated with each reaction.
  }
  \item{species}{
    A character vector which denotes the species names.
  }
  \item{constants}{
    Not supported yet.
    %A character vector which denotes the model constants. The model
    %constants are substituted in the C code by their numeric value.
  }
}

\value{
  Returns a list with the following elements:
  \item{ccode}{
    The actual C code as text.
  }
  \item{cspecies}{
    A character vector of the species names in the C code and
    their corresponding model names are given in the \code{names}
    attribute.
  }
  \item{cthetas}{
    A character vector of the names of the stochastic constants
    in the C code and their corresponding model names are given
    in the \code{names} attribute.
  }
  \item{Cov}{
    A character vector of the functions names that corresponds the
    upper triangular matrix (given as the names attribute) of the
    instantaneous variance-covariance matrix.
  }
  \item{Means}{
    A character vector of the functions names that corresponds to
    the instantaneous means (given as the names attribute).
  }
  \item{Orders}{
    A numerical vector indicating the order of each reaction.
  }

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
  Vasileios Giagos
}
%\note{
%%  ~~further notes~~
%}
\details{
  The derivation of the ODEs  is accomplished with the help of
  Deriv package and some symbolic algebra routines written in R. The
  derived ODEs are joined together, i.e. treated as a single system, and
  are expressed as C code. In particular, they are expressed as a C
  function of the following form:
  \preformatted{
double * name (double * t, double * y, double * fout, double * vthetas)}
  }
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  After parsing you need to compile the model with
  \code{\link{compmod}}.
}
\examples{
##Parses the Lotka-Volterra Model.

#Stoichiometry matrix
tt <- matrix(c(1,-1,0,0,1,-1),nrow=2,ncol=3,byrow=TRUE)
#The reaction rates
rfun <- c("con1 * Prey","con2 * Prey * Predator","con3 * Predator")
#Parameters
thetas <- paste("con",1:3,sep="")
#Species
species <- c("Prey","Predator")
\dontrun{
cout <- parsemod(tt,rfun,thetas,species)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
