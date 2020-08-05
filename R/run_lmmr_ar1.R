#' run_lmmr_ar1
#' @param species the species name, either YFS (yellowfin sole) or POP (pacific ocean perch)
#' @param niters the number of NUTS Hamiltonion iterations to use
#' @return Stan model object
#' 
#' @export speciesLMMR.Rdata the stan model object, sent to Outputs subdirectory
#'
#' @examples
#' 
run_lmmr_ar1 <- function(species, niters) {
require(rstan)
require(shinystan)
require(lme4)
n_chains = 3



if (species == "YFS")  filename <- "data/YFS_all.csv"
if (species == "POP") filename <- "data/POP_meas.csv"
thedata <- read.csv(file = filename, header = T)
thedata$FishID <- as.factor(thedata$FishID)
thedata$Year <- as.factor(thedata$Year)



# set up model matrixes
x1 <- model.matrix(~Increment_age , data = thedata)
x2 <- model.matrix(~-1 + Year, data = thedata)
x <- cbind(x1, x2)

z <- model.matrix(~-1+FishID, data = thedata)
nfish <- ncol(z)
nyear <-ncol(x2)

nbeta <- 2
ngamma <- ncol(z)
y <- log(thedata$Increment) - mean(log(thedata$Increment)) # log, so that growth anomalies are on same scale as vonbertalanffy model
ndata <- length(y)

stan.data <- list(y= y,
                  x = x,
                  z = z,
                  ngamma = ngamma,
                  nbeta = nbeta,
                  nyear = nyear,
                  ndata = ndata)

initf2 <- function(chain_id = 1, nbeta, ngamma) {
  list(
    beta = as.array(runif(nbeta, -1.0, 1.0)),
    gammaraw = rnorm(ngamma, 0, 1),
    logit_rho = runif(1, -1, 1),
    sigma = runif(n = 1, 0.1, 1),
    sigma_eta = runif(1, 0.1, 1),
    sigma_fish = runif(n = 1, 0.1, 1),
    year_epsilon_raw = rnorm(nyear,0, 1)

  )
}

params = c(
  "beta",
  "yearbeta",
  "sigma",
  "sigma_eta",
  "sigma_fish",
  "sigma_year",
  "rho",
  "logit_rho",
  "gamma",
  'y_pred_check'
)

model.name <- "Stan/LinearMixedEffect.stan"
# get initial values for mcmc
init_ll <- lapply(1:n_chains, function(id)
  initf2(chain_id = id, nbeta = nbeta, ngamma = ngamma))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# run STAN model
model.stan = stan(
  file = model.name,
  data = stan.data,
  iter = niters,
  pars = params,
  warmup = floor(niters / 2),
  chains = n_chains,
  thin = 1,
  algorithm = 'NUTS',
  init = init_ll,
  verbose = FALSE,
  control = list(adapt_engaged = TRUE, adapt_delta = 0.99, max_treedepth = 10)
)

savefilename <- paste("Outputs/",species,"LMMER.Rdata", sep = "")
save(file = savefilename, model.stan)
return(model.stan)
}

     