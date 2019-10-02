### Compute the discrete wavelet coefficients for the 17 subjects
##  averaging groups of the same stimulus by subject to reduce noise.

rm(list=ls())
### Load the needed libraries
library (R.matlab)
library (wavelets)

source ("parameters_sig.r")

addresses=c("/VOT/Passivo/","/VOT/Ativo/","/Formantes/Passivo/","/Formantes/Ativo/")

########################################################################
for (ad in 1:4){
    for (side in seq (1, 2)) {

        for (subj in subjects) {
### Initialize the output matrix to store the DWT coefficients and the
### position counter
            pos <- 1
            dwt.coefs <- matrix (NA, nrow = groups * nstim, ncol = nb.points+1)
            for (stim in seq (1, nstim)) {

                ## Read individual stimulus (Stim2 is baseline corrected and begins at the time 0s of the trigger)
                eeg.file <- file.path (sprintf (paste("../Sujeito%d",addresses[ad],"Stim2_%d.mat",sep=""), subj, stim))
                file <- readMat(eeg.file)
                signals <- file$samples2[,1:nb.points,]

### Compute amount of trials by group
                g <- ncol(signals[side,,])%%groups
                tr.group<- (ncol(signals[side,,])-g)/groups

                begin <- 1
                end <- tr.group
                
                for (subgroup in 1:groups){

                    ## Compute mean value
                    mean.resp <- colMeans (t (signals[side,,begin:end]))

                    begin <- end+1
                    end <- end+tr.group

                    cf <- stim
                   
                    cf <- c (cf, mean.resp)
                    
                    dwt.coefs [pos, ] <- cf
                    pos <- pos+1

                    cat (sprintf ("\rAddress: %d subject: %2d  side: %d   stimulus: %d    subgroup: %2d",ad, subj, side, stim, subgroup))
                    flush (stdout ())
                }#subgroup

            }#stim

            save(dwt.coefs, file=file.path(sprintf(paste("../Sujeito%d",addresses[ad],"sig_grouped_%s.dat",sep="") ,subj,c("left","right")[side])))
        }#subject

    }#side
    
}#adresses

cat ("\n")
flush (stdout ())

############################################################################################################
