\name{randomBeta}
\alias{randomBeta}

\title{
Generating selection gradients/vectors in random directions
}

\description{
\code{randomBeta} generates unit length vectors (selection gradients) uniformly distributed in a k-dimensional hypersphere.
}

\usage{
randomBeta(n = 1, k = 2)
}

\arguments{
  \item{n}{
Number of selection gradients/vectors
}
  \item{k}{
Number of dimensions 
}
}

\details{
\code{randomBeta} exploits the spherical symmetry of a multidimensional Gaussian density function. Each element of each vector is randomly sampled from a univariate Gaussian distribution with zero mean and unit variance. The vector is then divided by its norm to standardize it to unit length.
}

\value{
\code{randomBeta} returns a matrix where the vectors are stacked column wise. 
}

\author{
Geir H. Bolstad
}

\examples{
randomBeta(n = 5, k = 3)
}

\keyword{array}
\keyword{algebra}
