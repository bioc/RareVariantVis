\name{movingAverage}
\alias{movingAverage}

\title{
Computation of moving average
}
\description{
Function calculates moving average from a vector of numeric values.
}
\usage{
movingAverage(x, n, centered)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{a vector of numeric values for which moving average is computed}
    \item{n}{numeric value giving the frame length for moving average}
    \item{centered}{logic variable indicating if moving average should be centered (default = FALSE)}
}

\value{
    \item{comp1}{function returns vector of moving average values}
}

\author{
Winston Chang
}

\examples{
movingAverage(1:20, n=3, centered=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~moving average }

