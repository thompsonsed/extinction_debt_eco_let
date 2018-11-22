# Preston function and related solutions
library(Bessel)

# The Preston Function

#' Preston function for approximating species richness from a spatially explicit neutral model under
#' the provided parameters
#'
#' @param A_on_sigma_sq area (number of individuals) divided by the variance of dispersal
#' @param nu the speciation rate
#'
#' @return the number of species
Preston <- function(A_on_sigma_sq,nu)
{
  a = sqrt(A_on_sigma_sq/pi)
  nu_eff = nu*log(1/nu)/(1-nu)
  num = 2*pi*a*(1-nu_eff)*BesselI(a,1,expon.scaled=T)
  denom = (1/sqrt(nu_eff))*BesselI(a,1,expon.scaled=T)*
    (BesselK(sqrt(nu_eff)*a,0,expon.scaled=T)/
    BesselK(sqrt(nu_eff)*a,1,expon.scaled=T))+BesselI(a,0,expon.scaled=T)
  return(nu_eff*A_on_sigma_sq+num/denom)
}

# Calculate the contiguous result
S_contig<- function(A,nu,sigma_sq)
{
  sigma_sq*Preston(A/sigma_sq,nu)
}


#' Species richness from the Preston function immediately after clearing a contiguous landscape
#' in a random pattern.
#'
#' @param A_max the area (number of individuals) in the previous, contiguous habitat
#' @param A the area (number of individuals) in the current, randomly fragmented habitat
#' @param nu the speciation rate
#' @param sigma_sq the variance of the dispersal kernel
#'
#' @return the number of species
S_random <- function(A_max,A,nu,sigma_sq)
{
  sigma_sq*Preston(A/sigma_sq,1-(1-nu)^(A_max/A))
}

#' Species richness from the Preston function at equilibrium after clearing a contiguous landscape
#' in a random pattern.
#'
#' @param A_max the area (number of individuals) in the previous, contiguous habitat
#' @param A the area (number of individuals) in the current, randomly fragmented habitat
#' @param nu the speciation rate
#' @param sigma_sq the variance of the dispersal kernel
#'
#' @return the number of species
S_random_equilibrium <- function(A_max, A, nu, sigma_sq)
{
  S_contig(A, nu, sigma_sq * (A/A_max))
}

#' Estimation of the sigma parameter of dispersal from a distance and the number of steps it took
#' to get there.
#' 
#' The sigma parameter is the standard deviation of a normal distribution (note that sigma^2 is 
#' used in other functions, but sigma is returned here).
#'
#' @param distance the distance travelled in n_steps
#' @param n_steps the number of steps
#'
#' @return the sigma parameter estimation
estimate_sigma_rayleigh <- function(distance, n_steps)
{
  return((distance) * (2/(pi* n_steps))^0.5)
}

#' Predits the number of species remaining from a simple power-law species-area curve
#' 
power_law_estimation <- function(A_max, A, z)
{
  return((A/A_max)^z)
}
