## Clean up last session.
Sys.time()
unlink(".RData")
gc()

## Libraries.
library("findpython")
library("argparse")
sessionInfo()


## Set up parser to handle arguments of the simulation (nobs, nnoise, cor).
parser <- ArgumentParser(description = "Handling arguments for simulation")
parser$add_argument("--nobs", "-o", metavar = "nobs", type = "integer",
                    help = "number of observations to be used from the simulated data set")
parser$add_argument("--nnoise", "-n", metavar = "nnoise", type = "integer",
                    help = "number of noise variables")
parser$add_argument("--rho", "-r", metavar = "rho", type = "integer",
                    help = "correlation between variables")

## Parsing input arguments. 
args <- parser$parse_args()
if (is.null(args$nobs) | is.null(args$nnoise) | is.null(args$rho)) {
  parser$print_help()
  stop(9)
} else {
  cat("Starting simulation with:\n")
  cat(sprintf("    nobs:    %10d\n", args$nobs))
  cat(sprintf("    nnoise:  %10d\n", args$nnoise))
  cat(sprintf("    rho:     %10d\n", args$rho))
  
}

if (interactive()) {
  SGE_TASK_ID <- 1L
} else {
  SGE_TASK_ID <- as.integer(Sys.getenv("SGE_TASK_ID"))
}
stopifnot(is.integer(SGE_TASK_ID), length(SGE_TASK_ID) == 1L)

## Make sure nobs, nnoise, and rho integers.
stopifnot(is.numeric(args$nobs),   length(args$nobs) == 1,   args$nobs %% 1 == 0)
stopifnot(is.numeric(args$nnoise), length(args$nnoise) == 1, args$nnoise %% 1 == 0)
stopifnot(is.numeric(args$rho),    length(args$rho) == 1,    args$rho %% 1 == 0)

## Show the current working directory. 
cat("My current working directory is", getwd(), ".\n")
cat("\n")

## Create output directory if not existing. 
outdir <- "../simulation_results"
if (!dir.exists(outdir)) { 
  dir.create(outdir)
}

## Now define the output file name which contains ALL the param settings
outputfile <- file.path(outdir, sprintf("ZANBI_rep_ID_%04d_OBS_%06d_N_%03d_COR_%03d.rds", 
                                        SGE_TASK_ID, args$nobs, args$nnoise, args$rho))



## Run simulation and check if not already existing. 
if (file.exists(outputfile)) {
  cat("file", outputfile, "exists.\n\n Nothing to do exit. Check files.\n")
} else {
  cat("Working on", outputfile, "\n\n")
  ## Load sim_ZANBI.R 
  source("sim_ZANBI.R")
  set.seed(SGE_TASK_ID)
  ## Start simulation with nrep = SGE_TASK_ID, nobs = nobs, nnoise = noise. 
  res <- sim(
    nobs = args$nobs,
    nnoise = args$nnoise,
    rho = args$rho
  )
}

##Summary. 
ls()
print(res)

## Save results as .rds.
saveRDS(res, outputfile)

## Clean up.
gc(verbose = FALSE)
warnings()
rm(list = ls(all.names = TRUE)) 

