
require(rstan)
require(shinystan)
require(KernSmooth)
library(RColorBrewer)
require(readxl)
require(dplyr)

spc <- "POP"

if (spc =="POP") load("Outputs/POP_result.Rdata")




if (spc == "YFS") load("Outputs/YFS_result.Rdata")
output <- extract(model.stan)
# load environmental time series
per.2.use <- 0.8 # width of credibility interval

eps.qs <- output$eps_q 
beta.t <- output$beta_t

# adjust for beta.t
for (i in 1:ncol(eps.qs)) eps.qs[,i] <- beta.t * eps.qs[,i]

n.years <- ncol(eps.qs)

# a place to save median and upper and lower bound
fitted.qs <- matrix(NA, nrow = n.years, ncol = 3)
colnames(fitted.qs) <- c("lower","median","upper")

# loop through columns, fit smoother, and report percentiles
for (i in 1:n.years) {
  q.year <- eps.qs[,i]
  q.smooth <- bkde(q.year,
                   kernel = "normal",
                   canonical = FALSE)
  delta.q <- q.smooth$x[2] - q.smooth$x[1]
  # generate cumulative probability
  cumprob <- rep(NA, length (q.smooth$x))
  for (j in 1:length(q.smooth$x)) cumprob[j] <- sum(q.smooth$y[1:j])* delta.q
  
  # find lower percentile
  per <- c( (1- per.2.use) / 2, 0.5, 0.5 + 0.5 * per.2.use)
  fitted.qs[i,] <- approx(cumprob, q.smooth$x, 
               xout = per,
               method = "linear")$y
}

# figure out year indices - important!
if (spc == "POP") filename <- "data/POP_meas.csv"
if (spc == "YFS") filename <- "data/YFS_all.csv"
thedata <- read.csv(file = filename, header = T)
min.first.age <- min(thedata$Increment_age)
thedata$FishIDNum <- as.numeric(as.factor(thedata$FishID))
thedata$CohortNum <- as.numeric(as.factor(thedata$Cohort))
n.fish <- max(thedata$FishIDNum)
thedata$InitAge <- rep(NA, nrow(thedata))
# create a new column that shows initial age for each fish
for (i in 1:n.fish) {
  fish.index <- which(thedata$FishIDNum==i)
  thedata$InitAge[fish.index] <- rep(min(thedata$Increment_age[fish.index]), length(fish.index))
}

initial.year <- min(thedata$Year - thedata$InitAge + min.first.age) # this is the earliest possible year.
final.year <- max(thedata$Year) # this is the latest possible year

real.year.list <- initial.year : final.year


# calculate polygones
bottom.poly<-fitted.qs[,1]
top.poly<-fitted.qs[,3]
span=.07
degree=2

smooth <- T

if (smooth) {
  family="gaussian"
  
  bottom.poly.smooth<-loess.smooth(real.year.list, bottom.poly, span = span, degree = degree,
                                   family = family, evaluation = length(real.year.list) * 4)
  top.poly.smooth<-loess.smooth(real.year.list, top.poly, span = span, degree = degree,
                                family = family, evaluation = length(real.year.list)*4)
  
  x.poly<-c(bottom.poly.smooth$x,rev(top.poly.smooth$x))
  y.poly<-c(bottom.poly.smooth$y,rev(top.poly.smooth$y))
} else {
  x.poly <- c(real.year.list, rev(real.year.list))
  y.poly <- c(fitted.qs[,1], rev(fitted.qs[,3]))
}


## get the number of samples in each year

if (spc == "POP") xlims <- c(1905, 2015)
if (spc == "YFS")  xlims <- c(1950, 2015)
samples.years <- xlims[1]:xlims[2]
n.per.year <- rep(NA, length(samples.years))
for (i in 1:length(samples.years)){
  n.per.year[i] <- length(which(thedata$Year==samples.years[i]))
}

#---------------------------------------------
# make plot


plotfilename <- paste("graphics/",spc,".pdf", sep = "")
plot.matrix <- c(1,1,2,2,2,2,2,3,3,3,3,3)

pdf(plotfilename, 
    height = 6,
    width = 4)
layout(plot.matrix)
par(mar = c(3,5,0,0), oma = c(3,3,3,3))
ylims <- c(-.6, 0.6)

# make sample size plot
plot(samples.years, n.per.year,
     type = "l",
     lwd = 2,
     xlab = "",
     ylab = "",
     axes = F,
     xlim = xlims,
     ylim = c(0, max(n.per.year)),
     xaxs = "i",
     yaxs = "i",
     xpd = NA,
     cex.axis = 1.5,
     cex.lab = 1.5)
axis(2, at = c(0, max(n.per.year)), las = 1, cex.axis = 1.5)

# make fitted q plot 
col.2.plot<-col2rgb("black", alpha = 0.2)
col.2.plot <- rgb(col.2.plot[1]/255, col.2.plot[2]/255, col.2.plot[3]/255, alpha = 0.2)
plot(real.year.list,fitted.qs[,2],
     type="n",
     lwd=2,
     main="",
     xlab="",
     ylab="",
     xlim=xlims,
     ylim = ylims,
     xaxs="i",
     yaxs="i",
     axes=T,
     las = 1,
     cex.axis = 1.5,
     cex.lab = 1.5)


polygon(x.poly,y.poly,col=col.2.plot,border=col.2.plot)
lines(real.year.list,fitted.qs[,2],lwd=2,col = "black")


# now, plot the LMER result
rm(output)
if (spc == "POP")  load("Outputs/POPLMMER.Rdata")
if (spc == "YFS") load("Outputs/YFSLMMER.Rdata")

output <- extract(model.stan)
real.year.list <- min(thedata$Year): max(thedata$Year)
beta.t <- output$yearbeta
n.years <- ncol(beta.t)
fitted.betas <- matrix(NA, nrow = n.years, ncol = 3)
for (i in 1:n.years) {
  beta.year <- beta.t[,i]
  beta.smooth <- bkde(beta.year,
                   kernel = "normal",
                   canonical = FALSE)
  delta.beta <- beta.smooth$x[2] - beta.smooth$x[1]
  # generate cumulative probability
  cumprob <- rep(NA, length (beta.smooth$x))
  for (j in 1:length(beta.smooth$x)) cumprob[j] <- sum(beta.smooth$y[1:j])* delta.beta
  
  # find lower percentile
  per <- c( (1- per.2.use) / 2, 0.5, 0.5 + 0.5 * per.2.use)
  fitted.betas[i,] <- approx(cumprob, beta.smooth$x, 
                          xout = per,
                          method = "linear")$y
}
fitted.betas <- fitted.betas - mean(fitted.betas[,2])

### make plotting polygons
# calculate polygones
bottom.poly<-fitted.betas[,1]
top.poly<-fitted.betas[,3]
span=.075
degree=2

smooth <- T

if (smooth) {
  family="gaussian"
  
  bottom.poly.smooth<-loess.smooth(real.year.list, bottom.poly, span = span, degree = degree,
                                   family = family, evaluation = (length(real.year.list)-1) * 4)
  top.poly.smooth<-loess.smooth(real.year.list, top.poly, span = span, degree = degree,
                                family = family, evaluation = length(real.year.list)*4)
  
  x.poly<-c(bottom.poly.smooth$x,rev(top.poly.smooth$x))
  y.poly<-c(bottom.poly.smooth$y,rev(top.poly.smooth$y))
} else {
  x.poly <- c(real.year.list, rev(real.year.list))
  y.poly <- c(fitted.betas[,1], rev(fitted.betas[,3]))
}


col.2.plot<-col2rgb("black", alpha = 0.2)
col.2.plot <- rgb(col.2.plot[1]/255, col.2.plot[2]/255, col.2.plot[3]/255, alpha = 0.2)
plot(real.year.list,fitted.betas[,2],
     type="n",
     lwd=2,
     main="",
     xlab="",
     ylab="",
     xlim=xlims,
     ylim = c(-.6, 0.6),
     xaxs="i",
     yaxs="i",
     axes=T,
     las = 1,
     cex.axis = 1.5,
     cex.lab = 1.5)


polygon(x.poly,y.poly,col=col.2.plot,border=col.2.plot)
lines(real.year.list,fitted.betas[,2],lwd=2,col = "black")
mtext(text = "Year Effect", side =2, outer = TRUE)
mtext(text = "Year", side =1, outer = TRUE)
dev.off()
system2("open", args = c("-a Skim.app", plotfilename))



