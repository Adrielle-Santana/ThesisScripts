### Compute the discrete wavelet coefficients for the 17 subjects
##  averaging groups of the same stimulus by subject to reduce noise.

rm(list=ls())
### Load the needed libraries
library (R.matlab)
library (wavelets)

source ("parameters_sig.r")

addresses=c("/VOT/Passivo/","/VOT/Ativo/","/Formantes/Passivo/","/Formantes/Ativo/")
groups <- 0
########################################################################
for (ad in 1:4){
    for (side in seq (1, 2)) {

        for (subj in subjects) {
### Initialize the output matrix to store the DWT coefficients and the
### position counter
            nb.trials <- 0
        for (stim in seq (1, nstim)) {
            eeg.file <- file.path (sprintf (paste("../Sujeito%d",addresses[ad],"Stim2_%d.mat",sep=""), subj, stim))
            file <- readMat(eeg.file)
            nb.trials <- nb.trials + ncol(file$samples2[side,,])
        }
            pos <- 1
            dwt.coefs <- matrix (NA, nrow = nb.trials, ncol = nb.points+1)
            for (stim in seq (1, nstim)) {

                ## Read individual stimulus (Stim2 is baseline corrected and begins at the time 0s of the trigger)
                eeg.file <- file.path (sprintf (paste("../Sujeito%d",addresses[ad],"Stim2_%d.mat",sep=""), subj, stim))
                file <- readMat(eeg.file)

                trials<- ncol(file$samples2[side,,])

                for (trial in 1:trials) {

                    cf <- stim
                   
                    cf <- c (cf, file$samples2[side,1:nb.points,trial])
                    
                    dwt.coefs [pos, ] <- cf
                    pos <- pos+1

                    cat (sprintf ("\rAddress: %d subject: %2d  side: %d   stimulus: %d    trial: %3d",ad, subj, side, stim, trial))
                    flush (stdout ())
                }#trial

            }#stim

            save(dwt.coefs, file=file.path(sprintf(paste("../Sujeito%d",addresses[ad],"sig_grouped_%s.dat",sep="") ,subj,c("left","right")[side])))
        }#subject

    }#side
    
}#adresses

cat ("\n")
flush (stdout ())

############################################################################################################
