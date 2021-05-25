library(ggplot2)
library(gridExtra)
library(dplyr)

theme_set(theme_minimal())

species <- "YFS"

if (species == "POP") filename <- "data/POP_meas.csv"
if (species == "YFS") filename <- "data/YFS_all.csv"
thedata <- read.csv(file = filename, header = T)
thedata$FishID <- as.factor(thedata$FishID)

yfs.plot <- ggplot(thedata, aes(x = Increment_age, y = Increment, group = FishID, colour = FishID)) + geom_line() + geom_point(alpha = 0.5) +
  theme_bw() +
  scale_colour_viridis_d() +
  xlab("Age-at-Increment") + 
  xlim(0, 35) + 
  ylim(0, 0.125) +
  ylab("Increment (mm)") +
  theme(legend.position = "none") +
  theme(panel.grid.major= element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.text= element_text(size = 12)) +
geom_text(x = 0, y = 0.125, label = "A", colour = "black", size = 5)
        
# get min and max size by age

age.summary <- thedata %>%
  group_by(Increment_age) %>%
  summarise(min_inc = min(Increment), max_inc = max(Increment))
            
species <- "POP"

if (species == "POP") filename <- "data/POP_meas.csv"
if (species == "YFS") filename <- "data/YFS_all.csv"
thedata <- read.csv(file = filename, header = T)
thedata$FishID <- as.factor(thedata$FishID)

pop.plot <- ggplot(thedata, aes(x = Increment_age, y = Increment, group = FishID, colour = FishID)) + geom_line() + geom_point(alpha = 0.5) +
  theme_bw() +
  scale_colour_viridis_d() +
  xlab("Age-at-Increment") +
  xlim(0, 80) + 
  ylim(0, 0.04) + 
  ylab("Increment (mm)") +
  theme(legend.position = "none") +
  theme(panel.grid.major= element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.text= element_text(size = 12)) +
  geom_text(x = 0, y = 0.04, label = "B", colour = "black", size = 5)

plotfilename <- "Graphics/increment_vs_age.pdf"
pdf(file = plotfilename,
    width = 5,
    height = 2.75)

grid.arrange(yfs.plot, pop.plot, nrow = 1)
dev.off()
system2("open", args = c("-a Skim.app", plotfilename))


  
  
  
