% Generated by roxygen2 (4.0.2): do not edit by hand
\name{MLPAiter}
\alias{MLPAiter}
\title{Iterate the MLPA model}
\usage{
MLPAiter(params, num_times, ivec, N0 = 1)
}
\arguments{
\item{params}{A list with all the biological and biophysical parameters
needed to run the model. Will normally be generated by makeMLPAparams()}

\item{num_times}{The number of time steps to iterate}

\item{ivec}{An integer vector with length num_times, with each element
specifying the index of the dispersal kernel (third dimension of params$kernel)
to use in each simulation year. Will normally be generated with sample();
included so that the same sequence of years can be anaylzed across model variants.}

\item{N0}{Initial population abundance in each patch, in units of number of
 individuals per unit of habitat. All age classes in all habitats are intialized
 to the same abundance}
}
\value{
A list with two components. Nt is a three-dimensional array with the
 abundance in each age class at each location at each time. The dimensions are
 location, age, and time respecively. It includes time zero, so the third
 dimension has length num_times+1.

 The second component is Wt, which is a two-dimensional array of spwaning
    biomass, arranged by space and time.
}
\description{
Iterate the MLPA model
}

