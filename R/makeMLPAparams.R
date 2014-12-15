#' Create parameter list for size structured model
#' 
#' @param m Instantaneous natural mortality rate
#' @param Linf Asymptotic size for von Bertelanffy growth model
#' @param k Initial growth rate for von Bertelanffy growth model
#' @param x0 Initial size for von Bertelanffy growth model (may be negative)
#' @param c Coefficent for transforming size into biomass
#' @param d Exponent for tranforming size (length) into biomass
#' @param xM Age at first reproduction
#' @param xmax Maximum age
#' @param kernel Collection of square dispersal kernels, arrange in a three-dimensional array with source locations in columns, destination locations in rows, and years along the third dimension of the array
#' @param CR The compensation ratio (sets intensity of density dependence; for some reason, it is four for all species!)
#' @param scale Equilibrium spawning stock biomass in a deterministic, nonspatial model with unit area
#' @param habitat The amount of habitat in each location. May be a scalar (all locations the same) or a vector with length matching the first dimension of kernel
#' 
#' @return A list containing all the parameters needed to run MLPAiter

makeMLPAparams <- function(m, Linf, k, x0, c, d, xM, xmax, kernel, CR=4, scale=100, 
                           habitat=1) {
  # Annual survival
  p <- exp(-m)
  
  # Cumulative survival to age x, conditioned on surviving ot age 1.
  cum_p <- cumprod(rep(p,xmax))/p 
  
  # Size at age
  x <- 1:xmax
  Lx <- Linf*(1 - exp(-k*(x-x0)))
  wx <- c*Lx^d
  
  # alpha
  alpha <- CR/drop(cum_p[-(1:(xM-1))] %*% wx[-(1:(xM-1))])
  
  # Adjust alpha to take into accout the dispersal matrix
  Kmean <- apply(kernel, 1:2, mean)
  lambda <- Re(eigen(Kmean)$values[1])
  alpha <- alpha/lambda
  
  # Set beta. If habitat is a scalar, then beta will be a scalar and will be applied in 
  #  all patches. If habitat is a vector, it desctibes the relative amount of recruitment 
  #  habitat in each patch, and beta will be a vector. The input paramter scale sets the 
  #  desired equilibrium spawning biomass in a nonspatial population with habitat=1.
  beta <- (CR-1)/(scale/habitat)
  beta <- beta/lambda
  
  return(list(
    alpha=alpha,
    beta=beta,
    p=p,
    wx=wx,
    xM=xM,
    xmax=xmax,
    kernel=kernel))
}
