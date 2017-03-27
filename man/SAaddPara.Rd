\name{SAaddPara}
\alias{SAaddPara}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SAaddPara()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    cat(c("Type the name of the parameter which sensitivity you want to analyse: \n"), 
        fill = TRUE)
    namePara <- scan(, what = "text")
    cat(c("Write down the values ", namePara, " may assume.\n (return blank when done)\n"), 
        fill = TRUE)
    fndPara <- scan()
    fndDist <- fitdistr(fndPara, "normal")
    GoFfndDistr <- ks.test(fndPara, "pnorm", fndDist$estimate[1], 
        fndDist$estimate[2])
    scrPara <- combn(fndPara, (length(fndPara) - 1))
    singParEff <- unlist(lapply(X = seq(1, length(fndPara)), 
        FUN = deltDens))
    cat(c(namePara, "seems normally distributed with a probability aestimed at ", 
        GoF2fndDistr$p.value, "\% \n each one of the values provided affect the distribution for ", 
        1 - mean(singParEff), "\% (mean)\n", 1 - max(singParEff), 
        "\% (max) \n do you want to save it for further exploration? \n (y|n)"), 
        fill = T)
    promptGo <- scan(what = "text", nmax = 1)
    while (promptGo != "y" & promptGo != "n") {
        cat("answer y or n \n ")
        promptGo <- scan(, what = "text", nmax = 1)
    }
    if (promptGo == "y") {
        npDist <- data.frame(param = namePara, dist = "norm", 
            P1 = fndDist$estimate[1], P2 = fndDist$estimate[2])
        if (any(ls(SAsobEN) == "parDists")) {
            rbind(SAsobEN$parDists, npDist)
        }
        else {
            SAsobEN$parDists <- npDist
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
