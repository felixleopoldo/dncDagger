library(ggplot2)

timing <- read.csv("results/timesparticles.csv")
dir.create("results/figures")
# Timings
graphics.off()
ggplot(timing, aes(x=as.factor(n), y=totaltime, col=as.factor(d))) + geom_boxplot()
ggsave("results/figures/boxplots_timing.png")

graphics.off()
ggplot(timing, aes(x=as.factor(n), y=log(totaltime), col=as.factor(d))) + geom_boxplot()
ggsave("results/figures/boxplots_timing_log.png")

# Max particles
graphics.off()
ggplot(timing, aes(x=as.factor(n), y=max_particles, col=as.factor(d))) + geom_boxplot()
ggsave("results/figures/boxplot_maxparticles.png")

graphics.off()
ggplot(timing, aes(x=as.factor(n), y=log(max_particles), col=as.factor(d))) + geom_boxplot()
ggsave("results/figures/boxplots_maxparticles_log.png")

# Tot particles
graphics.off()
ggplot(timing, aes(x=as.factor(n), y=tot_particles, col=as.factor(d))) + geom_boxplot()
ggsave("results/figures/boxplots_totparticles.png")

graphics.off()
ggplot(timing, aes(x=as.factor(n), y=log(tot_particles), col=as.factor(d))) + geom_boxplot()
ggsave("results/figures/boxplots_totparticles_log.png")

