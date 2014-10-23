\name{matrix2long}
\alias{matrix2long}

\title{Convert a quadratic matrix to long format}
\description{
This function converts data from a quadratic round robin matrix into long format.
}

\usage{
matrix2long(M, new.ids=TRUE, var.id="value")
}

\arguments{
	\item{M}{A matrix with actors in rows and partners in columns)}
	\item{new.ids}{Should new ids for actors and partners be defined? (If new.ids=FALSE the row and column names are taken. In that case, you have to make sure, that rows and columns have the same set of names.)}
	\item{var.id}{The name of the column with the measured variable in the returned data frame}
}

\value{
A data frame in long format
}

\seealso{
\code{\link{long2matrix}}
}

\examples{
#The example data are taken from the "Mainz Freshman Study" and consist 
# of ratings of liking as well as ratings of the metaperception of 
# liking at zero-acquaintance using a Round-Robin group of 54 participants 
# (Back, Schmukle, & Egloff, in pres)

# load a data set in matrix style
data("liking_a")

str(liking_a)
long <- matrix2long(liking_a)

str(long)
 
}