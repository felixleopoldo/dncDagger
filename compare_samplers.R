pgibbs_output <- readRDS("pgibbs_output.rds")
omcmc_output <- readRDS("omcmc_output.rds")

oscores <-  omcmc_output$traceadd$orderscores[1:1000]
pscores <- pgibbs_output$pgibbs_log_scores[seq(1,10000,10)]

plot(oscores,ylim = range(c(-45350, -45200)), col="red")
points(pscores, col="blue")

max(oscores)
min(pscores)