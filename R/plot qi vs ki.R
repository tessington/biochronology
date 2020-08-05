

require(rstan)
require(KernSmooth)
require(dplyr)


getci <- function(mcmc, per.2.use) {
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
make.qi.ki.plot <- function(spc, per.2.use) {
  if (spc =="POP") load("Outputs/POP_result.Rdata")

if (spc == "YFS") load("Outputs/YFS_result.Rdata")
output <- extract(model.stan)

kis <- output$eps_k
qis <- output$eps_qi
kbase <- output$k_base
qbase <- output$q_base
psi <- output$psi
psi.out <- getci(psi, per.2.use)
print(psi.out)

log.ki <- log.qi <- matrix(NA, nrow = nrow(kis), ncol = ncol(kis))

for (i in 1:ncol(kis)) log.ki[,i] <- kis[,i] + log(kbase)
for (i in 1:ncol(qis)) log.qi[,i] <- qis[,i] + log(qbase) + psi * kis[,i]

log.qi.out <- log.ki.out <- matrix(NA, nrow = ncol(kis), ncol = 3)


for (i in 1:ncol(kis)) {
  log.ki.out[i,]<- getci(log.ki[,i], per.2.use)
  log.qi.out[i,]<- getci(log.qi[,i], per.2.use)
  
}

# simple plot
plot(x = log.ki.out[,2],
     y = log.qi.out[,2],
     pch = 21,
     bg = "black", 
     xlab = expression(paste("log(k"["i"],")", "")),
     ylab = expression(paste("log(q"["i"],")", ""))
     )

qbase.out <- getci(qbase,per.2.use)
kbase.out <- getci(kbase,per.2.use)

min.k <- min(log.ki.out[,2])
max.k <- max(log.ki.out[,2])

min.q <- log(qbase.out[2]) + psi.out[2] *( min.k - log(kbase.out[2]))
max.q <- log(qbase.out[2]) + psi.out[2] *( max.k - log(kbase.out[2]))

lines(x = c(min.k, max.k),
     y = c(min.q, max.q),
     lwd= 2)
}

plotfilename <- "qivski.pdf"
pdf(file = plotfilename,
    height = 3.5,
    width = 6,
    useDingbats = F)

spc <- "YFS"
par(mfrow = c(1,2), las = 1)

per.2.use <- 0.8 # width of credibility interval

make.qi.ki.plot(spc,per.2.use)

spc <- "POP"
make.qi.ki.plot(spc,per.2.use)
dev.off()
system2("open", args = c("-a Skim.app", plotfilename))