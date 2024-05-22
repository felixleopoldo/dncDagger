library(ggplot2)
library(dplyr)
# library(tidyverse)
library(latex2exp)
library(patchwork)

#timings <- read.csv("./paper_gobnilpjoint.csv") # this contains the runtimes used in the paper
timings <- read.csv("./paper_joint.csv") 
# dir.create("figures")

timings["d"] <- timings["d"]

algs <- c("dnc")
# algs <- c("dnc", "orderpruner")

for (alg in algs) {
    print(alg)
    
    timing <- timings %>% filter(alg %in% c(alg, "gobnilp"))

    ## Boxplots
    ##########################################


    timing <- timing %>% filter(d %in% c(0.0, 0.5, 1.0, 1.5, 2.0))

    # Timings
    graphics.off()
    ggplot(timing, aes(x = as.factor(d), y = totaltime, col = as.factor(alg))) +
        geom_boxplot() +
        xlab("p") +
        ylab("Time (s.)") +
        theme_bw() +
        labs(col = "d") +
        theme(text = element_text(size = 12), axis.title = element_text(size = 14))

    ggsave(paste0("figures/", alg, "_gobnilp_boxplots_timing.png"))

    graphics.off()
    ggplot(timing, aes(x = as.factor(d), y = log2(totaltime), col = as.factor(alg))) +
        geom_boxplot() +
        xlab("p") +
        ylab(TeX("$ \\log_2(t)$")) +
        theme_bw() +
        labs(col = "d") +
        theme(axis.title = element_text(size = 14))

    ggsave(paste0("figures/", alg, "_gobnilp_boxplots_timing_log.png"))

    graphics.off()

    df1 <- timing[timing$alg == "gobnilp", ]
    df2 <- timing[timing$alg == "dnc", ]
    # create a dataframe with differences in log times between dnc and gobnilp
    # df1$logtime <- log2(df1$totaltime)
    # df2$logtime <- log2(df2$totaltime)
    # df1$d <- df1$d
    # df2$d <- df2$d
    # df1$alg <- "gobnilp"
    # df2$alg <- "dnc"
    # df1 <- rbind(df1, df2)


    # use ggplot to plot log time of dnc against gobnilp log time    
    
    ggplot(df1, aes(x = log2(totaltime), y = log2(df2$totaltime), col = as.factor(d))) +
        geom_point() +
        xlab(TeX("$ \\log_2(t)$ GOBNILP")) +
        ylab(TeX("$ \\log_2(t)$ d&c pruner")) +
        xlim(c(-1, 3)) +
        ylim(c(-7, 12)) +
        labs(col="d") + 
        geom_abline(
            intercept = 0,
            slope = 1,
            col = "black",
            linetype = "dashed"
        ) +         
        theme_bw() +
      geom_abline(intercept =0, slope = 1, linetype = "dashed") +
        theme(axis.title = element_text(size = 14))

    ggsave(paste0("figures/", alg, "_gobnilp_boxplots_time_vs_log.png"), width = 7, height = 5)

    graphics.off()

    # plot the difference in log times of dnc and gobnilp, with respect to d
    ggplot(df1, aes(x = as.factor(d), y = log2(df2$totaltime) - log2(totaltime), col=as.factor(d))) +
        geom_boxplot() +
        xlab("d") +
        labs(col="d") + 
        ylab(TeX("$ \\log_2(t)$ d&c pruner - $ \\log_2(t)$ GOBNILP")) +
        geom_jitter( size=0.4, alpha=0.9) +
        theme_bw() +
        geom_hline(yintercept=0,  linetype = "dashed") +
        theme(axis.title = element_text(size = 14)) #+ coord_fixed(ratio = 0.1, clip="on")

    ggsave(paste0("figures/", alg, "_gobnilp_boxplots_time_diff.png"), width = 7, height = 5)

}

