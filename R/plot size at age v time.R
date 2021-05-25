require(rstan)
require(KernSmooth)
require(viridis)

####### Plotting Function #####

interp.fun <- function(win, df) {
  df$h.at.age <- df$Lstart / df$wstart
  # find the right row
  index <- which(df$wstart>=win)[1]
  if (is.na(index)) h.predict <- df$h.at.age[nrow(df)]
  if (!is.na(index)) h.predict <- approx(x = df$wstart[c(index-1, index)], y = df$h.at.age[c(index-1, index)], xout = win)$y
  return(h.predict)
}

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
           years.2.trim = 10, 
           age.ref) {
 require(dplyr)
    
       output <- extract(model.stan)
    
    if (spc == "POP") {
      Linf <- 41.55
      k <-  0.14
      tzero <- -1.317
      ages <- 2:90
      df <- readRDS("Outputs/POP_length_otolith.RDS")
    }
    
    if (spc == "YFS") {
      Linf <- 33.7
      k <- 0.151
      tzero <- -0.111
      ages <- 2:35
      df <- readRDS("Outputs/YFS_length_otolith.RDS")
      
    }
    df$h.at.age <- df$Lstart / df$wstart
    
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
    
    
    h.at.age <- df$Lstart / df$wstart
    winf.t <- q.t / k.bar
  
    
    
    # make a matrix of size at age, where rows are ages, columns are ages
    n.years <- length(winf.t)
    n.ages <- length(ages)
    
    woutput <- matrix(data = NA,
                     nrow = n.ages,
                     ncol = n.years)
    w.t.start <- rep(x = NA, times = n.ages)
    w.t.start[1] <- winf.bar[1]*(1 - exp(-(ages[1] - tzero) * k.bar))
    
    for (a in 2:n.ages)
      w.t.start[a] <-
      w.t.start[a - 1]  + (1 - exp (-k.bar)) * (winf.t[1]  - w.t.start[a - 1])
    woutput[, 1] <- w.t.start
    for (i in 2:n.years) {
      woutput[1, i] <- w.t.start[1]
      for (a in 2:n.ages)
        woutput[a, i] <-
          woutput[a - 1, i - 1]  + (1 - exp (-k.bar)) * (winf.t[i]  - woutput[a - 1, i -
                                                                              1])
    }
    # now scale output based on h.at.age
    

    
   # print(ages.2.use + 5 - 2)
  #  print(dim(output))
    minsize <- (apply(X=woutput[ages.2.use[1:4] + 5 - 2,], MARGIN = 1, FUN = min))
    maxsize <- (apply(X=woutput[ages.2.use[1:4] + 5 - 2,], MARGIN = 1, FUN = max))
    print(maxsize - minsize)
    
    for (i in 1:4) {
      if (i == 1)
    
        plot(
          1:n.years + min.year+years.2.trim,
          df$h.at.age[df$ages == age.ref] * woutput[ages.2.use[i] + 5 - 2, ],
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
              df$h.at.age[df$ages == age.ref]* woutput[ages.2.use[i] + 5 - 2, ],
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
plotfilename = "Graphics/compare_fit_to_obs.pdf"
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
  make.plot(spc, output, thedata, min.year, max.year, small.data, ages.2.use, legend.text, legend.pos = "topright", years.2.trim = 5, age.ref = 20)


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
  age.ref <- 70
  make.plot(spc, output, thedata, min.year, max.year, small.data, ages.2.use, legend.text, legend.pos = "bottomright", age.ref = 70)

dev.off()
system2("open", args = c("-a Skim.app", plotfilename))


  