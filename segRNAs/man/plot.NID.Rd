% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seg_criteria.R
\name{plot.NID}
\alias{plot.NID}
\title{Plotting NID object}
\usage{
\method{plot}{NID}(nids)
}
\arguments{
\item{nid}{: the NID objet returned by the function nid.penalties}
}
\value{
the graph according to the parameters entered in the function
}
\description{
Plotting NID object
}
\examples{
nid.beta <- nid.penalties(rna.data1 = dataset1 , rna.data2 = dataset2, penalty_range = c(15,75) )
plot(nid.beta)

}
