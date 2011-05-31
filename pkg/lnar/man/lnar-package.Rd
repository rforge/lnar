\name{lnar-package}
\alias{lnar-package}
\alias{lnar}
\docType{package}
\title{
  Inference for stochastic kinetic genetic networks using the Linear
  Noise Approximation.}
\description{
  A collection of functions that implement the Linear Noise Approximation
  for stochastic kinetic models with emphasis on genetic auto-regulatory
  networks.
}
\details{
 \tabular{ll}{
  Package: \tab lnar\cr
  Type: \tab Package\cr
  Version: \tab 0.0.5\cr
  Date: \tab 2011-05-31\cr
  License: \tab GPL (>=2.0) \cr
  LazyLoad: \tab yes\cr
 }
 %~~ An overview of how to use the package, including the most important ~~
 %~~ functions ~~
 The \pkg{lnar} package provides inferential tools for a class of genetic
 auto-regulatory networks based on the Linear Noise Approximation (Kurtz
 1972). Two LNA-based estimation methods are provided: the
 \emph{Restarting} and the \emph{Non Restarting} method, see (Giagos 2010)
 for more details. Such networks, are specified as a system of biochemical
 reactions in the \code{\link{parsemod}} method which, in turn, outputs
 the underlying Linear Noise Approximation as C code to be compiled with
 the \code{\link{compmod}} method.  The compiled model can be fitted to
 a dataset using \code{\link{optmod}}, a Maximum Likelihood Estimation
 procedure. Try \code{demo(lv)} for an example implementing the
 Lotka-Voltera model and \code{demo(autoreg)} for the implementation of
 a prokaryotic transcription model (Golightly and Wilkinson 2005). 
} 
\note{
  All methods expect the parameters to be expressed in terms of
  \eqn{\theta}{thetas}, i.e.  scaled according to their order. Normally, in
  a biological model, e.g. a SBML file, the parameters (\eqn{c}{c})
  correspond to kinetics equations based on the number of molecules.

}
\section{Warning}{
  This is an experimental and unstable package. Most of the C code has
  been ported from an earlier version implemented in C using the Gnu
  Scientific Library (GSL).
}
\author{
 Vasileios Giagos \email{v.giagos@kent.ac.uk}
}
\references{
  Kurtz, T. G.: 1972, \emph{The relationship between stochastic and
    deterministic models for chemical reactions}, The Journal of Chemical
  Physics \bold{57}(7), 2976-2978.\cr
  
  Golightly, A. and Wilkinson D. J.: 2005, \emph{Bayesian Inference for
   Stochastic Kinetic Models Using a Diffusion Approximation},
   Biometrics \bold{61}(3), 781-788.

  Giagos, V.: 2010, \emph{Inference for auto-regulatory genetic networks
  using diffusion   process approximations}, Thesis, Lancaster University,
  2010.
}
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%ls()
%}
