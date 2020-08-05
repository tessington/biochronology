# this code compares fitted model to the model that best fits conventional size at age estimation.  
require(rstan)
require(KernSmooth)


plotfilename <- "CompareBCtoVB.pdf"
pdf(file = plotfilename, height = 3.5, width = 7)

load("Outputs/YFS_result.Rdata")
output <- extract(model.stan)

Linf.male <- 33.7
k.male <- 0.151
tzero.male <- -0.111
# feMale
Linf.female <- 37.8 
k.female <- 0.237
tzero.female <- 0.112
ages <- 2:34
Lstart.male  <- Linf.male * (1- exp(-(ages - tzero.male)*k.male))
Lstart.female  <- Linf.female * (1- exp(-(ages - tzero.female)*k.female))

# get otolith biochronogy-based estimated
q.bar <- median(output$q_base)
k.bar <- median(output$k_base)
wstart.bar <- median(output$w_start)
Winf.bar <- q.bar / k.bar
wstart.width <- wstart.bar * Winf.bar

# get size at age of otliths
wstart <- rep(x = NA, times = length(ages))
wstart[1] <- wstart.width
for (j in 2:length(ages)) wstart[j] <-  wstart[j-1] + (1 - exp(-k.bar)) * (Winf.bar - wstart[j-1])

par(las = 1, mar = c(5,5,5,5), mfrow = c(1,2))
plot(x = ages,
     y = Lstart.male,
     type = "l",
     lwd = 2,
     col = "black",
     xlim = c(0,35),
     ylim = c(0, 40),
     ylab = "Length (cm)",
     xlab = "Age (years)")
lines(x = ages, y = Lstart.female, lwd = 2)

# scale the two so it makes sense.  scaling to size at age 30
h <- (Lstart.male[30]/wstart[30])
w.scaled <- wstart * (Lstart.male[30]/wstart[30])
lines(x = ages, y = w.scaled, lwd = 2, col = "gray")
w.ticks <- seq(from = 0, to = 1.5, by = 0.25)
axis( side = 4, at = w.ticks *h, labels = w.ticks )
mtext(side = 4, text = "Otolith width (cm)", las = 0, line = 3)

# repeat for pacicic ocean pearch
spc <- "POP"

  load("Outputs/POP_result.Rdata")
  output <- extract(model.stan)
  
Linf <- 41.55
k <-  0.14
tzero <- -1.317
ages <- 2:83
Lstart  <- Linf * (1- exp(-(ages - tzero)*k))
q.bar <- median(output$q_base)
k.bar <- median(output$k_base)
wstart.bar <- median(output$w_start)
Winf.bar <- q.bar / k.bar
wstart.width <- wstart.bar * Winf.bar

wstart <- rep(x = NA, times = length(ages))
wstart[1] <- wstart.width
for (j in 2:length(ages)) wstart[j] <-  wstart[j-1] + (1 - exp(-k.bar)) * (Winf.bar - wstart[j-1])


plot(x = ages,
     y = Lstart,
     type = "l",
     lwd = 2,
     col = "black",
     xlim = c(0,80),
     ylim = c(0, 80),
     ylab = "Length (cm)",
     xlab = "Age (years)")

# scale the two so it makes sense to plot on the same graph.  scaling to size at age 30
h <- (Lstart[30]/wstart[30])
w.scaled <- wstart * (Lstart[30]/wstart[30])
lines(x = ages, y = w.scaled, lwd = 2, col = "gray")
# make the second y axis, based on scaling
w.ticks <- seq(from = 0, to = 1.5, by = 0.25)
axis( side = 4, at = w.ticks *h, labels = w.ticks )
mtext(side = 4, text = "Otolith width (cm)", las = 0, line = 3)






# extract median qbar, kbar
dev.off()
system2("open", args = c("-a Skim.app", plotfilename))
