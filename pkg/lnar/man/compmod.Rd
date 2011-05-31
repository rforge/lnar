\name{compmod}
\alias{compmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compile parsed models.
}
\description{
Compiles the generated C source code of a parsed model using
the \pkg{inline} package.
}
\usage{
compmod(cout, name = "derivs")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cout}{
    A parsed model.
}
  \item{name}{
    A string indicating the name of the compiled function, defaults
    to "derivs".
}
}
\details{
  Uses the \pkg{inline} package to compile the generated c code.
}
\value{
  Returns a compiled function named as the \code{name} argument  and
  is included in R's environment \bold{automatically}, i.e. \emph{no
    assignment is needed}. The compiled function
  has the following form in C:
  \preformatted{
double * name ( double * t, double * y, 
                double * fout, double * vthetas) }
}
%\references{
%% ~put references to the literature/web site here ~
%}
\note{
  We note again that there is no need to assign the return value to a
  variable. The \code{compmod} function loads the compiled function
  directly to the R's environment.}
\author{
  Vasileios Giagos
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  The function \code{\link{parsemod}} generates the c code for the
  models. The transition density can be calculated with
  \code{\link{calcdens}} and MLEs are obtained with
  \code{\link{optmod}}.
}
\examples{
\dontrun{
##After generating a model "cout"
compmod(cout,"derivs") #no need to assign the return value
}
}
