# this code compares fitted model to the model that best fits conventional size at age estimation.  
require(rstan)
require(KernSmooth)


plotfilename <- "Graphics/CompareBCtoVB.pdf"
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
ages <- 2:60
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
     ylim = c(0, 50),
     ylab = "Length (cm)",
     xlab = "Age (years)")
lines(x = ages, y = Lstart.female, lwd = 2)

par(new = TRUE)
h <- (Lstart.male[30]/wstart[30])
w.scaled <- wstart * (Lstart.male[30]/wstart[30])
plot(x = ages, y = wstart, type = "l",
     lwd = 2, col = "gray", xlab = "", ylab = "", xlim = c(0, 35),
     ylim = c(0, 1.5), axes = F, )

w.ticks <- seq(from = 0, to = 1.5, by = 0.25)
axis(side = 4, at = w.ticks, labels = w.ticks)
mtext(side = 4, text = "Otolith width (cm)", las = 0, line = 3)
# finally, save as dataframe for later use
df <- as.data.frame(cbind(ages, Lstart.male, wstart))
saveRDS(object = df, file = "Outputs/YFS_length_otolith.RDS")


par(new = FALSE)
# repeat for pacific ocean pearch
spc <- "POP"

  load("Outputs/POP_result.Rdata")
  output <- extract(model.stan)
  
Linf <- 41.55
k <-  0.14
tzero <- -1.317
ages <- 2:180
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
     ylim = c(0, 50),
     ylab = "Length (cm)",
     xlab = "Age (years)")

par(new = TRUE)
plot(x = ages, y = wstart, type = "l", lwd = 2, col = "gray", xlab = "", ylab = "", xlim = c(0, 80),
     ylim = c(0, 1.5), axes = F)
w.ticks <- seq(from = 0, to = 1.5, by = 0.25)
axis(side = 4, at = w.ticks, labels = w.ticks)
mtext(side = 4, text = "Otolith width (cm)", las = 0, line = 3)
df <- as.data.frame(cbind(ages, Lstart, wstart))
saveRDS(object = df, file = "Outputs/POP_length_otolith.RDS")
par(new = FALSE)







# extract median qbar, kbar
dev.off()
system2("open", args = c("-a Skim.app", plotfilename))
