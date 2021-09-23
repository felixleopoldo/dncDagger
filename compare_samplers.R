# pgibbs_output <- readRDS("pgibbs_output.rds")
omcmc_output <- readRDS("omcmc_output.rds")
##mh_output <- readRDS("mh_output.rds")

oscores <-  omcmc_output$traceadd$orderscores
# pscores <- pgibbs_output$pgibbs_log_scores[seq(1,10000,10)]
# #mhscores <- mh_output$mh_log_scores[seq(1,300000000,300000)]
#mhscores <- mh_output$mh_log_scores[seq(1,10000000,1)]
#mh_output

max(oscores)
plot(oscores)
length(oscores)

# plot(oscores,ylim = range(c(-45350, -45150)), col="red", ylab="Log-score")
# points(pscores, col="blue")
# points(mhscores, col="green")
# title("Order log-scores BGE")
# legend(x="topright",legend=c("BiDAG order MCMC","PGibbs", "MH local swap move"), col=c("red", "blue", "green"),title="Algorithm", lty=1:4, cex=0.8)
# max(mhscores)
# max(oscores)
# max(pscores)



# acf(oscores, lag.max = 20)

