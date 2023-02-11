library(ggplot2)
library(dplyr)
library(tidyverse)
#library(Hmisc)

timing <- read.csv("results/timesparticles.csv")
dir.create("figures")



#timing <- timing %>% filter(d %in% c(2.0))
#print(timing)

ds <- timing %>% distinct(d)

#####################
# Max particles     #
#####################
reglist <- list()

as <- c()
bs <- c()
astderrs <- c()
bstderrs <- c()
mylm = 1

label <- "tot #particles"
label <- "max particles"
for (dd in ds[["d"]]) {
  mylm <- lm(log(max_particles) ~ n + log(n), timing %>% filter(d == dd))

  stderrs <- sqrt(diag(vcov(mylm)))
  astderrs <- c(astderrs, stderrs[2])
  bstderrs <- c(bstderrs, stderrs[3])
  
  a <- coef(mylm)[["n"]]
  b <- coef(mylm)[["log(n)"]]
  as <- c(as, a)
  bs <- c(bs, b)
  reglist[toString(dd)] <- summary(mylm) 
}



## Reg coeffs with errorbars ########
plotdf <- data.frame(a=as,b=bs, aerr=astderrs, berr=bstderrs, d= ds[["d"]], coef="a")


ggplot(plotdf, aes(x=d, y=a)) + 
  ggtitle(paste0("for varying d: log(", label, ") ~  a*n + b*log(n) + c")) +
  xlab("graph (d)ensity") + ylab("a") +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=a-aerr, ymax=a+aerr), width=.2,
                position=position_dodge(0.05))

ggsave("figures/aplots.png")

ggplot(plotdf, aes(x=d, y=b)) + 
  ggtitle(paste0("for varying d: log(", label, ") ~  a*n + b*log(n) + c")) +
  xlab("graph (d)ensity") + ylab("b") +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=b-berr, ymax=b+berr), width=.2,
                position=position_dodge(0.05))

ggsave("figures/bplots.png")


## Reg coeffs ########
png("figures/reglines.png")
par(mfrow=c(2,1))

#errbar(as.factor(ds[["d"]]), as, y+error_values,y-error_values,type='b')
plot(as.factor(ds[["d"]]), as, pch = 1, xlab = "graph (d)ensity", ylab = "a")
title(paste0("for varying d: log(", label, ") ~  a*n + b*log(n) + c")) 

plot(as.factor(ds[["d"]]), bs, pch = 1, xlab = "graph (d)ensity", ylab = "b")
title(paste0("for varying d: log(", label, ") ~  a*n + b*log(n) + c")) 

dev.off()

## Boxplots
###########################################

#timing <- timing %>% filter(d %in% seq(1.2, 2.0, 0.2))
timing <- timing %>% filter(d %in% c(0.0, 0.5 , 1.0,1.1, 1.5,0.9, 2.0))
# Timings
graphics.off()
ggplot(timing, aes(x=as.factor(n), y=totaltime, col=as.factor(d))) + geom_boxplot()
ggsave("figures/boxplots_timing.png")

graphics.off()
ggplot(timing, aes(x=as.factor(n), y=log(totaltime), col=as.factor(d))) + geom_boxplot()
ggsave("figures/boxplots_timing_log.png")

# Max particles
graphics.off()
ggplot(timing, aes(x=as.factor(n), y=max_particles, col=as.factor(d))) + geom_boxplot() 
ggsave("figures/boxplot_maxparticles.png")

graphics.off()
ggplot(timing, aes(x=as.factor(n), y=log(max_particles), col=as.factor(d))) + geom_boxplot() + geom_smooth(method = "lm", se = FALSE)
ggsave("figures/boxplots_maxparticles_log.png")

# Tot particles
graphics.off()
ggplot(timing, aes(x=as.factor(n), y=tot_particles, col=as.factor(d))) + geom_boxplot()
ggsave("figures/boxplots_totparticles.png")

graphics.off()
ggplot(timing, aes(x=as.factor(n), y=log(tot_particles), col=as.factor(d))) + geom_boxplot()
ggsave("figures/boxplots_totparticles_log.png")

