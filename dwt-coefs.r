### Compute the discrete wavelet coefficients for the 17 subjects
##  ungrouped.

rm(list=ls())
### Load the needed libraries
library (R.matlab)
library (wavelets)

source ("parameters.r")

addresses=c("/VOT/Passivo/","/VOT/Ativo/","/Formantes/Passivo/","/Formantes/Ativo/")
groups <- 0

########################################################################
for (ad in 1:4){
    for (side in seq (1, 2)) {

        for (subj in subjects) {
### Initialize the output matrix to store the DWT coefficients and the
### position counter
### Compute the number of trials per subject, across all stimuli
        nb.trials <- 0
        for (stim in seq (1, nstim)) {
           eeg.file <- file.path (sprintf (paste("../Sujeito%d",addresses[ad],"Stim2_%d.mat",sep=""), subj, stim))
            file <- readMat(eeg.file)
            nb.trials <- nb.trials + ncol(file$samples2[side,,])
        }
            pos <- 1
            dwt.coefs <- matrix (NA, nrow = nb.trials, ncol = nb.points+1)
            for (stim in seq (1, nstim)) {

                ## Read individual stimulus
                eeg.file <- file.path (sprintf (paste("../Sujeito%d",addresses[ad],"Stim2_%d.mat",sep=""), subj, stim))
                file <- readMat(eeg.file)
                
                trials<- ncol(file$samples2[side,,])

                for (trial in 1:trials) {

                    ## Compute DWT
                    wt <- dwt (file$samples2[side,1:nb.points,trial], n.levels = nb.dwt.levels, filter=mw, boundary = "periodic")

                    ## Extract DWT coefficients
                    cf <- stim
                    for (i in 1 : nb.dwt.levels) {
                        band <- sprintf ("W%d", i)
                        cf <- c (cf, slot (wt, "W") [[band]])
                    }
                    band <- sprintf ("V%d", nb.dwt.levels)
                    cf <- c (cf, slot (wt, "V") [[band]])
                    
                    dwt.coefs [pos, ] <- cf
                    pos <- pos+1

                    cat (sprintf ("\rAddress: %d subject: %2d  side: %d   stimulus: %d    trial: %3d",ad, subj, side, stim, trial))
                    flush (stdout ())
                }#trials

            }#stim
###filename=sprintf("../Sujeito%d/VOT/Passivo/grouped_%s.mat",subj,c("left","right")[side])
###writeMat(filename,dwt_coefs=dwt.coefs)

            save(dwt.coefs, file=file.path(sprintf(paste("../Sujeito%d",addresses[ad],"grouped_%s.dat",sep="") ,subj,c("left","right")[side])))
        }#subject

    }#side
    
}#adresses

cat ("\n")
flush (stdout ())

############################################################################################################
