library(ggplot2)
library(dplyr)
library(latex2exp)
library(patchwork)
#library(quantreg) . #not working on R 4.2.3  
timings <- read.csv("paper_joint.csv") # change this
#dir.create("figures")
timings["d"] <- timings["d"]

algs <- c("dnc")
#algs <- c("dnc", "orderpruner", "gobnilp")

for (alg in algs){
  print(alg)
  
  timing <- timings[timings["alg"]==alg,]
  
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
  
  label <- bquote(N[p])
  
  for (dd in ds[["d"]]) {
    #print(timing %>% filter(d == dd))
    mylm <- lm(log2(tot_particles) ~ n + log2(n), timing %>% filter(d == dd)) # Could change to max_particles
    
    stderrs <- sqrt(diag(vcov(mylm)))
    astderrs <- c(astderrs, stderrs[2])
    bstderrs <- c(bstderrs, stderrs[3])

    a <- coef(mylm)[["n"]]
    b <- coef(mylm)[["log2(n)"]] # gets NA for some reason
    as <- c(as, a)
    bs <- c(bs, b)
    reglist[toString(dd)] <- summary(mylm)
  }
  
 # Not working on R 4.2.3 
  # aqs <- c()
  # bqs <- c()
  # aqstderrs <- c()
  # bqstderrs <- c()
  
  # for (dd in ds[["d"]]) {
  #   myqm <- rq(log2(tot_particles) ~ n + log2(n), data = timing %>% filter(d == dd)) # Could change to max_particles
    
  #   stderrs <- sqrt(diag(vcov(mylm)))
  #   astderrs <- c(astderrs, stderrs[2])
  #   bstderrs <- c(bstderrs, stderrs[3])
    
  #   aq <- coef(myqm)[["n"]]
  #   bq <- coef(myqm)[["log2(n)"]] # gets NA for some reason
  #   aqs <- c(aqs, aq)
  #   bqs <- c(bqs, bq)
    
  #   temp <- rep(0, 4)
  #   try(temp <- summary(myqm)$coefficients[,2])
    
  #   aqstderrs <- c(aqstderrs, temp[2])
  #   bqstderrs <- c(bqstderrs, temp[3])
  # }
  
  ## Regression coefficients with errorbars 
  plotdf <- data.frame(a=as,b=bs, aerr=astderrs, berr=bstderrs, d= ds[["d"]], coef="a")
  
  p1 <- ggplot(plotdf, aes(x=d, y=a)) + 
    #xlab("graph (d)ensity") + ylab("a") +
    xlab("d") + ylab("a") +
    #ylim(-20,20) +
    geom_line(colour="blue") +
    geom_point(colour="blue")+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 14))+
    geom_errorbar(aes(ymin=a-aerr, ymax=a+aerr), width=.09, colour="blue")
  #                position=position_dodge(0.05))
  
  #ggsave("figures/aplots.png")
  #print(a)
  #print(b)
  p2 <- ggplot(plotdf, aes(x=d, y=b)) + 
    xlab("d") + ylab("b") +
    geom_line(colour="blue") +
    geom_point(colour="blue") +
    theme_bw() +
    #ylim(-300,300) +
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 14)) +
    geom_errorbar(aes(ymin=b-berr, ymax=b+berr), width=.09, colour="blue",
                  position=position_dodge(0.05))
  
  reg_plot <- p1 + p2 +  
    plot_annotation(title = TeX("$\\log_2(\\Sigma_p) \\sim  a*p + b*\\log_2(p) + c$")) &  theme(plot.title = element_text(hjust = 0.5)) 
  
  
  reg_plot 
  
  #ggsave("figures/bplots.png")
  ggsave(paste0("figures/",alg,"_regplots.png"), width = 14, height = 3.5)
  
  ### This requires the quantreg package which is not compatible with R 4.2.3 on the image
  ### However, it could be run locally.
  # plotdf <- data.frame(a=aqs,b=bqs, aerr=aqstderrs, berr=bqstderrs, d= ds[["d"]], coef="a")
  
  # p1 <- ggplot(plotdf, aes(x=d, y=a)) + 
  #   #xlab("graph (d)ensity") + ylab("a") +
  #   xlab("d") + ylab("a") +
  #   #ylim(-20,20) +
  #   geom_line(colour="blue") +
  #   geom_point(colour="blue")+
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 14))+
  #   geom_errorbar(aes(ymin=a-aerr, ymax=a+aerr), width=.09, colour="blue")
  # #                position=position_dodge(0.05))
  
  # #ggsave("figures/aplots.png")
  # #print(a)
  # #print(b)
  # p2 <- ggplot(plotdf, aes(x=d, y=b)) + 
  #   xlab("d") + ylab("b") +
  #   geom_line(colour="blue") +
  #   geom_point(colour="blue") +
  #   theme_bw() +
  #   #ylim(-300,300) +
  #   theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 14)) +
  #   geom_errorbar(aes(ymin=b-berr, ymax=b+berr), width=.09, colour="blue",
  #                 position=position_dodge(0.05))
  
  # reg_plot <- p1 + p2 +  
  #   plot_annotation(title = TeX("$\\log_2(\\Sigma_p) \\sim  a*p + b*\\log_2(p) + c$")) &  theme(plot.title = element_text(hjust = 0.5)) 
  
  # reg_plot 
  
  # #ggsave("figures/bplots.png")
  # ggsave(paste0("figures/",alg,"_regplots_med.png"), width = 14, height = 3.5)
 
  
  ## Reg coeffs ########
  png(paste0("figures/",alg,"_reglines.png"))
  par(mfrow=c(2,1))
  
  #errbar(as.factor(ds[["d"]]), as, y+error_values,y-error_values,type='b')
  plot(as.factor(ds[["d"]]), as, pch = 1, xlab = "graph (d)ensity", ylab = "a")
  title(paste0(TeX("For varying $d$: $\\log_2(N_p) \\sim  a*p + b*\\log_2(p) + c$"))) 
  #title(paste0("for varying d: log(N_p) ~  a*p + b*log(p) + c")) 
  
  plot(as.factor(ds[["d"]]), bs, pch = 1, xlab = "graph (d)ensity", ylab = "b")
  title(paste0(TeX("For varying $d$: $\\log_2(N_p) \\sim  a*p + b*\\log_2(p) + c$"))) 
  
  dev.off()
  
  ## Boxplots
  ##########################################
  
  timing_1 <- timing %>% filter(d %in% c(0.5, 1.0, 1.5, 2.0), n %% 2 == 0)
  # Timings
  graphics.off()
  #ggplot(timing, aes(x=as.factor(n), y=(totaltime/tot_particles), col=as.factor(d))) + geom_boxplot() #+ xlab("p") + ylab("Time (s.)") + labs(col = "d")
  ggplot(timing_1, aes(x=tot_particles, y=totaltime)) + geom_point() + xlab(label) + theme_bw() + ylab("Time (s.)") 
  ggsave(paste0("figures/",alg,"_boxplots_time_vs_totorders.png"), width = 7, height = 5)
  
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  
  # ratio
  graphics.off()
  ggplot(timing_1, aes(x=as.factor(n), y=(totaltime/tot_particles), col=as.factor(d))) + geom_boxplot() + xlab("p") + theme_bw() + ylab(TeX("$t / \\Sigma_N$")) + labs(col = "d") +
    theme(text = element_text(size = 12), axis.title = element_text(size = 14)) + 
    scale_color_manual(values=ggplotColours(5)[2:5])
  ggsave(paste0("figures/",alg,"_boxplots_n_vs_timeandtotorders.png"), width = 7, height = 5)
  
  # swrt ratio
  graphics.off()
  ggplot(timing_1, aes(x=as.factor(n), y=sqrt(totaltime/tot_particles), col=as.factor(d))) + geom_boxplot() + xlab("p") + theme_bw() + 
    #ylab("sqrt(Time (s.) / tot #orderings)") + 
    ylab(TeX("$(t / \\Sigma_N)^{1/2}$")) + 
    labs(col = "d") +
    theme(text = element_text(size = 12), axis.title = element_text(size = 14)) +
    scale_color_manual(values=ggplotColours(5)[2:5])
  ggsave(paste0("figures/",alg,"_boxplots_n_vs_sqrttimeandtotorders.png"), width = 7, height = 5)
  
  
  timing <- timing %>% filter(d %in% c(0.0, 0.5, 1, 1.5, 2.0))
  
  # Timings
  graphics.off()
  ggplot(timing, aes(x=as.factor(n), y=totaltime, col=as.factor(d))) + geom_boxplot() + xlab("p") + ylab("Time (s.)") + theme_bw() + labs(col = "d") +
    theme(text = element_text(size = 12), axis.title = element_text(size = 14)) 
  ggsave(paste0("figures/",alg,"_boxplots_timing.png"))
  
  graphics.off()
  ggplot(timing %>% filter(n %% 2 == 0), aes(x=as.factor(n), y=log2(totaltime), col=as.factor(d))) + geom_boxplot() + xlab("p") + ylab(TeX("$ \\log_2(t)$")) +theme_bw() + labs(col = "d") +
     theme(axis.title = element_text(size = 14)) 
  ggsave(paste0("figures/",alg,"_boxplots_timing_log.png"), width = 7, height = 5)

  graphics.off()
  ggplot(timing, aes(x=as.factor(n), y=log2(totaltime), col=as.factor(d))) + geom_boxplot() + xlab("p") + ylab(TeX("$ \\log_2(t)$")) +theme_bw() + labs(col = "d") +
    theme(axis.title = element_text(size = 14)) 
  ggsave(paste0("figures/",alg,"_boxplots_timing_log_wide.png"), width = 14, height = 7)
    
  # Max particles
  graphics.off()
  ggplot(timing, aes(x=as.factor(n), y=max_particles, col=as.factor(d))) + geom_boxplot() + xlab("p") + ylab("Max #orderings") +theme_bw() + labs(col = "d")
  ggsave(paste0("figures/",alg,"_boxplot_maxparticles.png"))
  
  graphics.off()
  ggplot(timing, aes(x=as.factor(n), y=log2(max_particles), col=as.factor(d))) + geom_boxplot() + geom_smooth(method = "lm", se = FALSE) +theme_bw() + xlab("p") + ylab("log(max #orderings)") + labs(col = "d")
  ggsave(paste0("figures/",alg,"_boxplots_maxparticles_log.png"))
  
  # Tot particles
  graphics.off()
  ggplot(timing, aes(x=as.factor(n), y=tot_particles, col=as.factor(d))) + geom_boxplot() + xlab("p") + ylab(TeX("$\\Sigma_N$")) +theme_bw() + labs(col = "d")
  ggsave(paste0("figures/",alg,"_boxplots_totparticles.png"))
  
  graphics.off()
  ggplot(timing %>% filter(n %% 2 == 0), aes(x=as.factor(n), y=log2(tot_particles), col=as.factor(d))) + geom_boxplot() + xlab("p") + ylab(TeX("$\\log_2(\\Sigma_N)$")) +theme_bw() + labs(col = "d") +
       theme(axis.title = element_text(size = 14))
  ggsave(paste0("figures/",alg,"_boxplots_totparticles_log.png"), width = 7, height = 5)
  
  graphics.off()
  ggplot(timing, aes(x=as.factor(n), y=log2(tot_particles), col=as.factor(d))) + geom_boxplot() + xlab("p") + ylab(TeX("$\\log_2(\\Sigma_N)$")) +theme_bw() + labs(col = "d") +
    theme(axis.title = element_text(size = 14))
  ggsave(paste0("figures/",alg,"_boxplots_totparticles_log_wide.png"), width = 14, height = 7)
}

