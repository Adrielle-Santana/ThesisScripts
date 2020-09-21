### This is the main script for computing and storing the DWT coefficients
### for all experiments.

### * Load needed local libraries
source ("paths.r")
source ("compute-dwt-coefs.r")
source ("experiments.r")

### * Load needed system libraries
load.pkgs ("parallel")

### * Compose the list of experiment features/types
### This will be used inside the parallelized procedure
cases <- list ()
count <- 1
for (ef in experiment.features)
    for (et in experiment.types) {
        cases [[count]] <- list (ef = ef, et = et)
        count <- count + 1
    }

### * Paralellized procedure
fun <- function (i) {
    case <- cases [[i]]
    compute.dwt.coefs (exp.feature = case$ef,
                       exp.type = case$et)
}

### * Main call

### Use the MC CORES environment variable for changing the amount of
### parallelization.  MC_CORES = 1 means serial execution.

##nb.cores <- Sys.getenv ("MC_CORES")
##if (nb.cores == "")
##    nb.cores <- detectCores ()
##return <- mclapply (seq (1, length (cases)),
##                    fun,
##                    mc.cores = nb.cores)

########################################
## Option to main call without cores

for (i in seq(1, length (cases))) fun(i)
########################################

### * Clean the progress meter
cat ("\n")
flush (stdout ())
