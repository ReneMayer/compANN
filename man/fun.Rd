\name{fun}
\alias{fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
insert here the title
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fun(wtmin, wtmex, hist, stimu, clampT, t, d, a, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wtmin}{
%%     ~~Describe \code{wtmin} here~~
}
  \item{wtmex}{
%%     ~~Describe \code{wtmex} here~~
}
  \item{hist}{
%%     ~~Describe \code{hist} here~~
}
  \item{stimu}{
%%     ~~Describe \code{stimu} here~~
}
  \item{clampT}{
%%     ~~Describe \code{clampT} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function (wtmin, wtmex, hist, stimu, clampT, t, d, a, n) 
.Call("filee8c1c7bd55f", wtmin, wtmex, hist, stimu, clampT, t, 
    d, a, n, PACKAGE = "compANN")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
