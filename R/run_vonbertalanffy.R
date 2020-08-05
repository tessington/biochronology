#' run_vonbertalanffy
#' @param species the species name, either YFS (yellowfin sole) or POP (pacific ocean perch)
#' @return Stan model object
#' 
#' @export species_result.Rdata the stan model object, sent to Outputs subdirectory
#'
#' @examples model.stan <- run_vonbertalanffy("YFS")
#' 
#' 
run_vonbertlanffy <- function(species){

if (species == "POP") filename <- "data/POP_meas.csv"
if (species == "YFS") filename <- "data/YFS_all.csv"
savefilename <- paste(species, "_result.Rdata",sep = "")


#savefilename <- "YFS_result_withstartprior.Rdata"
thedata <- read.csv(file = filename, header = T)
n.data <- nrow(thedata)
min.first.age <- min(thedata$Increment_age)
minage.fn <- function(x, fishids, fish.2.use) min(x[fishids==fish.2.use])

# turn fish id into 1...n
thedata$FishIDNum <- as.numeric(as.factor(thedata$FishID))
thedata$CohortNum <- as.numeric(as.factor(thedata$Cohort))
n.fish <- max(thedata$FishIDNum)
thedata$InitAge <- rep(NA, nrow(thedata))
# create a new column that shows initial age for each fish
for (i in 1:n.fish) {
  fish.index <- which(thedata$FishIDNum==i)
  thedata$InitAge[fish.index] <- rep(min(thedata$Increment_age[fish.index]), length(fish.index))
}

# make a vector of calendar years with corresponding year index
initial.year <- min(thedata$Year - thedata$InitAge + min.first.age) # this is the earliest possible year.
final.year <- max(thedata$Year) # this is the latest possible year

real.year.list <- initial.year : final.year


# function to map calendar year onto year index (year index is 1, 2,.... last year)
year.to.index <- function(cyear, iyear) return(which(iyear %in% cyear))
# function that uses year.to.index, to look up a fishid, extract the years, and translate into year index
make.year.index <- function(x, fishids, fish.2.use, iyear) return(year.to.index(x[fishids==fish.2.use], iyear))


minage.fn <- function(x, fishids, fish.2.use) min(x[fishids==fish.2.use])
startyears.fn <- function(x, fishids, fish.2.use) x[fishids==fish.2.use][1]



# turn calendar years into year index
thedata$iyear <- rep(NA, n.data)
for (i in 1:n.fish) thedata$iyear[thedata$FishIDNum == i] <- 
  make.year.index(thedata$Year, thedata$FishIDNum, i, real.year.list)

# get vector of startyear indexes and start ages
startyears <- startages <- rep(NA,n.fish)
for (i in 1:n.fish) startyears[i] <- startyears.fn(thedata$iyear,thedata$FishIDNum,i)
for (i in 1:n.fish) startages[i] <- minage.fn(thedata$Increment_age,thedata$FishIDNum,i)

# initialize other Stan data
n.years <- length(real.year.list)
min.first.age <- min(startages)
max.first.age <- max(startages)
firstyear <- startyears - startages + min.first.age + 1
# create cumulative increment for each fish
thedata$cinc <- rep(NA, n.data)
for (i in 1:n.fish) {
  winc <- thedata$Increment[thedata$FishIDNum == i]
  cinc <- rep(NA, length(winc))
  cinc[1]<-0
  for (t in 2:length(winc)) cinc[t] <- cinc[t-1]+winc[t-1]
  thedata$cinc[thedata$FishIDNum == i] <- cinc
}
#### OK, below  runs the Stan routine
require(rstan)
require(shinystan)

stan.data <- list(
  ndata = n.data,
  nfish = n.fish,
  nyears = n.years,
  maxfirstage = max.first.age,
  minfirstage = min.first.age,
  winc = thedata$Increment,
  years = thedata$iyear,
  startyears = startyears,
  fishids = array(thedata$FishIDNum),
  startage = startages,
  cinc = thedata$cinc,
  firstyear = firstyear
)


initf2 <- function(chain_id = 1, ndata, nfish, nyears) {
  list(
    beta_t = runif(1,0,.2),
    rho = runif(1,0.0,1),
    psi = runif(1, 0, 0.1),
    k_base = runif(1, .15, 0.151),
    q_base = runif(1, 0.5, 2),
    k_sigma = runif(1, 0, 0.1),
    eps_k_raw = rnorm(nfish, 0, 1),
    inc_sigma = runif(1, 0.01, 0.2),
    omega = rnorm(nyears, 0, 1),
    w_start_p = runif(1, 0, 0.25),
    eps_qi_raw = rnorm(nfish, 0, 1),
    qi_sigma = runif(1,0, 0.01)
  )
}
n_chains <- 3
init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id, ndata = n.data, nfish = n.fish, nyears = n.years))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

model.name <- "Stan/vonbertalanffy.stan"
params <- c("k_base",
            "q_base",
            "inc_sigma",
            "k_sigma",
            "qi_sigma",
            "w_start",
            "rho",
            "psi",
            "beta_t",
            "eps_q",
            "eps_k",
            "eps_qi",
            "winf_base"
)
model.stan <- stan(file = model.name,
                   data = stan.data,
                   iter = 1500,
                   par = params,
                   warmup = floor(1500/2),
                   chains = 3,
                   thin = 15,
                   algorithm = 'NUTS',
                   init = init_ll,
                   verbose = FALSE,
                   control = list(adapt_engaged = TRUE, adapt_delta = 0.9, max_treedepth = 10))

save(output, model.stan, stan.data, file = savefilename)
return(model.stan)
}

