### * Load the local libraries
source ("paths.r")
source ("dwt-lib.r")

### * Load the necessary R.matlab library
load.pkgs ("R.matlab")
load.pkgs ("signal")

### * Function for computing the DWT transformation of the EEG signals

compute.dwt.coefs <- function (exp.feature,
                               exp.type,
                               .cohort = cohort,
                               .dwt.filter = dwt.filter,
                               .dwt.length = dwt.length,
                               .dwt.nb.levels = dwt.nb.levels,
                               .baseline.length = baseline.length,
                               verbose = TRUE) {

    ## ** Number of stimuli
    nb.stim <- 5

    ## ** Create the indexation vectors
    idx.baseline <- seq (1, .baseline.length)

    ## ** Loop over the subjects
    for (subj in .cohort) {

        ## *** Create the output variables
        stimulus <- nb.trials <- c ()
        response <- signal.mean <- signal.sd <- rep (list (c ()), eeg.channels)

        ## *** Loop over the stimuli
        for (stim in seq (1, nb.stim)) {

            ## **** Read the data
            subj.file <- file.path (exp.dir,
                                    sprintf ("Sujeito%d", subj),
                                    exp.feature,
                                    exp.type,
                                    sprintf ("Channels%d.mat", stim))
            dat <- readMat (subj.file)

            ## **** Store the stimulus number
            n <- dim (dat$channels) [3]
            nb.trials <- c (nb.trials, n)
            stimulus <- c (stimulus, rep (stim, n))

            for (channel in seq (1, eeg.channels)) {

                ## ***** Initialize matrix of coefficients
                coef.mat <- signal.mat <- c ()

                ## ***** Loop over the responses
                for (k in seq (1, n)) {

                    ## ****** Get the signal (entire signal including baseline)
                    signal <- dat$channels [channel, , k]
                    signal <- resample(signal, eeg.sampfreq, eeg.origfreq)

                    ## ****** Baseline correction
                    baseline <- dat$channels [channel, idx.baseline, k]
                    baseline <- resample(baseline, eeg.sampfreq, eeg.origfreq)
                    signal <- signal - mean (baseline)

                    ## ****** Cumulate signals
                    signal.mat <- rbind (signal.mat, signal)

                    ## ****** Compute the wavelet transform
                    wt <- dwt (signal, filter = .dwt.filter,
                               n.levels = .dwt.nb.levels,
                               boundary = "periodic")
                    namesW <- names (wt@W)
                    namesV <- names (wt@V)

                    ## Align coefficients in time
                    wt <- align(wt)

                    ## Give original names which desapear after the align command
                    names (wt@W) <- namesW
                    names (wt@V) <- namesV
                    
                    ## length of the last level
                    remove <- 1
                    num.coefs <- length(wt@W[[.dwt.nb.levels]])
                    ## take coefficients from the borders

                    wt@V[[.dwt.nb.levels]] <- wt@V[[.dwt.nb.levels]][(remove+1):(num.coefs-remove), 1]
                    for (x in seq(.dwt.nb.levels,1)){

                        wt@W[[x]] <- wt@W[[x]][(remove+1):(num.coefs-remove), 1]
                        remove <- remove*2
                        num.coefs <- num.coefs*2

                    }
                    
                    wt <- dwt.to.coefs (wt)

                    ## ****** Cumulate coefficients
                    coef.mat <- rbind (coef.mat, wt$coefs)

                } # k

                ## ***** Store matrices
                response [[channel]] <- rbind (response [[channel]], coef.mat)
                signal.mean [[channel]] <- rbind (signal.mean [[channel]],
                                                  colMeans (signal.mat))
                signal.sd [[channel]] <- rbind (signal.sd [[channel]],
                                                apply (signal.mat, 2, sd))

                ## ***** Progress meter
                if (verbose) {
                    exp.label <- sprintf ("%25s, ",
                                          sprintf ("%s/%s",
                                                   exp.feature, exp.type))
                    cat (sprintf ("\r%ssubject %2d, stim %2d, channel %2d",
                                  exp.label, subj, stim, channel))
                    flush (stdout ())
                }

            } #channel

        } # stim

        ## *** Store results
        dwt.coefs <- list (stimulus = stimulus,
                           response = response,
                           mean = signal.mean,
                           sd = signal.sd,
                           nb.trials = nb.trials)

        ## *** Save the results
        save (dwt.coefs,
              file = dwt.coefs.filename (exp.feature, exp.type, subj))

    } # subj

    ## ** Clean the progress meter
    cat ("\n")

}
