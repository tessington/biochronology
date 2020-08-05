require(rstan)
getci <- function(mcmc, per.2.use) {
  require(KernSmooth)
  out <- rep(NA, times = 3)
  smoothed <- bkde(mcmc,
                   kernel = "normal",
                   canonical = FALSE)
  delta <- smoothed$x[2] - smoothed$x[1]
  # generate cumulative probability
  cumprob <- rep(NA, length (smoothed$x))
  for (j in 1:length(smoothed$x)) cumprob[j] <- sum(smoothed$y[1:j])* delta
  
  # find lower percentile
  per <- c( (1- per.2.use) / 2, 0.5, 0.5 + 0.5 * per.2.use)
  out <- approx(cumprob, smoothed$x, 
                xout = per,
                method = "linear")$y
  
  return(out)
  
}

species <- "YFS"
# full model
filename <- paste("Outputs/",species, "_result.Rdata", sep = "")
load(file = filename)


# find wich elements of output are for "beta_t"
output <- extract(model.stan)

posterior <- getci(output$rho, 0.8) # this is the variance of the anomaly time series
print(posterior)


# repeat for LMMER version
rm(model.stan)
rm(output)

filename <- paste("Outputs/", species, "LMMER.Rdata", sep = "")
load(file = filename)

output <- extract(model.stan)

posterior <- getci(output$rho, 0.8)
print(posterior)
