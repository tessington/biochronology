require(rstan)
require(KernSmooth)
require(viridis)

####### Plotting Function #####
make.plot <-
  function(spc,
           output,
           thedata,
           min.year,
           max.year,
           small.data,
           ages.2.use,
           legend.text,
           legend.pos,
           years.2.trim = 10) {
 require(dplyr)
    
       output <- extract(model.stan)
    
    if (spc == "POP") {
      Linf <- 41.55
      k <-  0.14
      tzero <- -1.317
      ages <- 2:90
    }
    
    if (spc == "YFS") {
      Linf <- 33.7
      k <- 0.151
      tzero <- -0.111
      ages <- 2:35
      
    }
    
    # extract median qbar, kbar
    q.bar <- median(output$q_base)
    k.bar <- median(output$k_base)
    
    beta.t <- output$beta_t
    eps.q.raw <- output$eps_q
    eps.q <- matrix(NA, nrow= nrow(eps.q.raw), ncol = ncol(eps.q.raw))
    # adjust for beta.t
    for (i in 1:ncol(eps.q.raw)) eps.q[,i] <- beta.t * eps.q.raw[,i]
    eps.q <- apply(eps.q, MARGIN = 2, median)
    eps.q <- eps.q[-(1:years.2.trim)]
    
    rm(output)
    q.t <- q.bar * exp(eps.q)
    winf.bar <- q.bar / k.bar
    h <- Linf / winf.bar
    winf.t <- q.t / k.bar
    Linf.t <- winf.t * h
    
    
    # make a matrix of size at age, where rows are ages, columns are ages
    n.years <- length(Linf.t)
    n.ages <- length(ages)
    
    output <- matrix(data = NA,
                     nrow = n.ages,
                     ncol = n.years)
    l.t.start <- rep(x = NA, times = n.ages)
    l.t.start[1] <- Linf * (1 - exp(-(ages[1] - tzero) * k.bar))
    
    for (a in 2:n.ages)
      l.t.start[a] <-
      l.t.start[a - 1]  + (1 - exp (-k.bar)) * (Linf.t[1]  - l.t.start[a - 1])
    output[, 1] <- l.t.start
    for (i in 2:n.years) {
      output[1, i] <- l.t.start[1]
      for (a in 2:n.ages)
        output[a, i] <-
          output[a - 1, i - 1]  + (1 - exp (-k.bar)) * (Linf.t[i]  - output[a - 1, i -
                                                                              1])
    }
    
   # print(ages.2.use + 5 - 2)
  #  print(dim(output))
    minsize <- (apply(X=output[ages.2.use[1:4] + 5 - 2,], MARGIN = 1, FUN = min))
    maxsize <- (apply(X=output[ages.2.use[1:4] + 5 - 2,], MARGIN = 1, FUN = max))
    print(maxsize - minsize)
    for (i in 1:4) {
      if (i == 1)
        plot(
          1:n.years + min.year+years.2.trim,
          output[ages.2.use[i] + 5 - 2, ],
          type = "l",
          lwd = 2,
          ylim = c(10, 70),
          col = col[1],
          xlab = "Year",
          ylab = "Length (cm)",
          xlim = c(1+min.year+years.2.trim - 5, n.years + min.year + years.2.trim)
        )
      
      if (i > 1)
        lines(1:n.years + min.year + years.2.trim,
              output[ages.2.use[i] + 5 - 2, ],
              lwd = 2,
              col = col[i])
      
      tmp.data <- small.data %>%
        filter(Capture_age_chron < ages.2.use[i + 1] &
                 Capture_age_chron >= ages.2.use[i])
      points(tmp.data$Year,
             tmp.data$Length,
             pch = 21,
             bg = col[i])
    }
    
   
    legend(legend.pos,
           pch = 21,
           pt.bg = col,
           legend = legend.text,
           cex = 0.75,
           bty = "n")
    
  }
#######
plotfilename = "compare_fit_to_obs.pdf"
pdf(file = plotfilename,
    height = 3.5,
    width = 7)
par(mfrow = c(1,2), las = 1)

spc = "YFS"
  load("Outputs/YFS_result.Rdata")
  filename <- "data/YFS_all.csv"
  thedata <- read.csv(file = filename, header = T)
  min.year <- min(thedata$Year)
  max.year <- max(thedata$Year)
  # get final size-at-age
  unique.ids <- which(duplicated(as.character(thedata$FishID)))
  small.data <- thedata[-unique.ids,]
  ages.2.use <- c(15,20,25,30,35)
  legend.text <-  c("15-20", "20-25","25-30", "30-35")
  col <- plasma(n=16)[c(2,6,10,16)]
  make.plot(spc, output, thedata, min.year, max.year, small.data, ages.2.use, legend.text, legend.pos = "topright", years.2.trim = 5)


spc ="POP"
  load("Outputs/POP_result.Rdata")
  filename <- "data/POP_meas.csv"
  thedata <- read.csv(file = filename, header = T)
  thedata$Length <- thedata$Length/10
  min.year <- min(thedata$Year)
  max.year <- max(thedata$Year)
  # get final size-at-age
  unique.ids <- unique(thedata$FishID)
  small.data <- thedata[unique.ids,]
  ages.2.use <- c(50,60,70,80,90)
  legend.text <-  c("50-60", "60-70","70-80", "80-90")
  make.plot(spc, output, thedata, min.year, max.year, small.data, ages.2.use, legend.text, legend.pos = "topleft")

dev.off()
system2("open", args = c("-a Skim.app", plotfilename))


  