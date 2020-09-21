### This script generates the roldsis coefficients for all subjects,
### continua, tasks and for physical and psychophysical responses.

### * Load local libraries
source ("paths.r")
source ("dwt-lib.r")
source ("experiments.r")
source ("roldsis-bootstrap.r")

### * Load system package
load.pkgs ("wavelets")

### * Physical responses of each subject for each continuum
phy.resp <- list ()
for (ef in experiment.features) {
    dat <- read.csv (file = file.path (data.dir, sprintf ("id%s.csv", ef)))
    phy.resp [[ef]] <- cbind (dat [, seq (2, 10, by = 2)])
}

electrodes <- read.table (file.path (data.dir, "electrodes.tab")) [, 1]
selected.electrodes <- c ("Fz", "F7", "F8", "TP9", "TP10")
index.electrodes <- sapply (selected.electrodes,
                            function (x) which (electrodes == x))

resp <- c ("phy", "psy")

### Specify position of the DWT coefficients
nb.wavelets <- 2 * length (dwt (rep (0, dwt.length))@W [[dwt.start.level]])
idx.wavelets <- seq (dwt.length - nb.wavelets + 1, dwt.length)

### * Initialize output variables
dirVec.boot <- list ()
dirVec <- list ()

boot.rep <- 100

for (ef in experiment.features) {

    ## Initialize feature slot
    dirVec [[ef]] <- list ()
    dirVec.boot [[ef]] <- list ()

    for (et in experiment.types) {

        ## Initialize type slot
        dirVec [[ef]] [[et]] <- list ()
        dirVec.boot [[ef]] [[et]] <- list ()

        for (subj in cohort) {

            ## Defines which matrix of physical responses to use
            output <- list (psy = psycho.resp,
                            phy = phy.resp [[ef]] [subj, ])

            subj.name <- sprintf("subj%d",subj)

            dirVec [[ef]] [[et]] [[subj.name]] <- list ()
            dirVec.boot [[ef]] [[et]] [[subj.name]] <- list ()

            cat (sprintf ("\r%s - %s  subj.: %d", ef, et, subj))
            flush (stdout ())

            load (dwt.coefs.filename (ef, et, subj))

            dat <- data.frame ()

            ## Start stim and electrodes loops to compose the big matrix
            ## to be used in the bootstrap procedure.
            for (stim in seq (1, 5)) {

                r <- c()
                idx <- 0
                cf <- 0

                idx <- which (dwt.coefs$stimulus == stim)

                ## For selected stimulus, compose the big matrix with all electrodes data
                for (elt in index.electrodes) {

                    cf <- dwt.coefs$response [[elt]] [idx, ]
                  
                    r <- cbind(r, cf [, idx.wavelets])
                }## elt

                ## Combines the big matrices for the 5 stimuli in a data.frame
                dat <- c(dat, data.frame(I(r)))

            }##stim

            ## Distribute responses in the structure separating each electrode response
            for (re in resp) {

                ## Performs bootstrap in the data.frame of the big matrices, compute the
                ## roldsis for physical and psychophysical responses. Observe that roldsis
                ## is computed for all electrodes together
                b <- roldsis.boot (dat, output [[re]], boot.rep)

                dirVec [[ef]] [[et]] [[subj.name]] [[re]] <- list()
                dirVec.boot [[ef]] [[et]] [[subj.name]] [[re]]  <- list()

                c <- 1
                ini <- 1

                for (elt in selected.electrodes) {

                    id <- seq (ini, length (idx.wavelets) * c)
                    c <- c + 1
                    ini <- id [length (idx.wavelets)] + 1

                    dirVec [[ef]] [[et]] [[subj.name]] [[re]] [[elt]] <- b$t0 [id]
                    dirVec.boot [[ef]] [[et]] [[subj.name]] [[re]] [[elt]] <- b$t [, id]

                } # elt

            } # resp

        } # subj

    } # et

} # ef

save (dirVec, file = file.path (results.dir, "dirVec.dat"))
save (dirVec.boot, file = file.path (results.dir, "dirVec-boot.dat"))
