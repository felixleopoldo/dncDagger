
library(argparser)
source("helper_functions.R")
p <- arg_parser("A program for generating and saving score tables.")
p <- add_argument(p, "--filename", help = "Data filename")
argv <- parse_args(p)

filename <- argv$filename
ret <- get_scores(filename)

saveRDS(ret, file = paste(filename, "rds", sep = "."))
