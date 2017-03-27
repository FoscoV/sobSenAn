\name{lwrDens}
\alias{lwrDens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lwrDens(parVal, shapeA1, shapeA2, shapeB1, shapeB2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parVal}{
%%     ~~Describe \code{parVal} here~~
}
  \item{shapeA1}{
%%     ~~Describe \code{shapeA1} here~~
}
  \item{shapeA2}{
%%     ~~Describe \code{shapeA2} here~~
}
  \item{shapeB1}{
%%     ~~Describe \code{shapeB1} here~~
}
  \item{shapeB2}{
%%     ~~Describe \code{shapeB2} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (parVal, shapeA1, shapeA2, shapeB1, shapeB2) 
{
    denA <- dnorm(parVal, shapeA1, shapeA2)
    denB <- dnorm(parVal, shapeB1, shapeB2)
    pmin(denA, denB)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
