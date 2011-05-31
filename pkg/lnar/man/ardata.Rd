\name{ardata}
\alias{ardata}
\docType{data}
\title{
SSA simulated data from the prokaryotic transcription model.
}
\description{
  The SSA algorithm was used to simulate the prokaryotic transcription
  model of Golightly and Wilkinson (2005) for a time interval of 49 
  time units. The data contain observations of the model's state at
  every time unit without any observational error.
}
\usage{data(ardata)}
\format{
  A data frame with 50 observations on the following 5 numeric variables.
  \describe{
    \item{\code{Time}}{The time instance the observation was recorded.}
    \item{\code{RNA}}{The number of molecules for the \bold{RNA} species.}
    \item{\code{P}}{The number of molecules for the \bold{P} species.}
    \item{\code{P2}}{The number of molecules for the \bold{P2} species.}
    \item{\code{DNA}}{The number of molecules for the \bold{DNA} species.}
  }
}
\details{
  The prokaryotic transcription model contains 5 species but due to
  the linear dependence between \bold{DNA} and \bold{DNA.P2}
  (\eqn{ {DNA} +  {DNA.P2} = 10}) , the latter was omitted.
}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
\references{
Golightly, A., Wilkinson, D. J.: 2005, \emph{Bayesian inference for
  stochastic kinetic models using a diffusion approximation}, 
  Biometrics, \bold{61}(3):781-788.
}
\examples{
data(ardata)
}
\keyword{datasets}
