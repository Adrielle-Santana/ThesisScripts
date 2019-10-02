### Compute the discrete wavelet coefficients for subjects 1 to 11 IC's to denoise and save IC's denoised

rm(list=ls())
### Load the needed libraries
library (R.matlab)
library (wavelets)

########################################################################
### Number of levels of the DWT
nb.dwt.levels <- 8
### Mother wavelet to be used
mw <- "la8"
### Initial values
subjects=seq(1:11)
nstim=5
nb.points=4096
nchannels=17
### Until what band to filter (from the first one), including it?
cutBand=6
########################################################################

## Initialize the output matrix to store the DWT coefficients and the
pastas=c('VOT/Passivo/','VOT/Ativo/','Formantes/Passivo/','Formantes/Ativo/')
Xica_DWT <- matrix (NA, nrow = nchannels, ncol = nb.points)

for (subj in subjects) {
    for (p in 1:4){
        for (stim in 1:nstim) {
            Xica.dir=paste(sprintf("../Sujeito%d/",subj),pastas[p],sprintf("Xica%d.mat",stim),sep="")
            
            ## Read individual stimulus
            file <- readMat(Xica.dir)

            for (ch in 1:nchannels){

                ## Compute DWT
                wt <- dwt (file$Xica[ch,], n.levels = nb.dwt.levels, filter=mw, boundary = "periodic")

### Noise elimination procedure recommended in the paper below for comparison of mothe wavelets through cross-correlation
### Biomedical Signal Analysis by DWT Signal De-noising with Neural Networks (2014)
#########################################################################################################################
                cf <- c()
                for (i in 1 : nb.dwt.levels) {
                    band <- sprintf ("W%d", i)
                    cf <- c (cf, slot (wt, "W") [[band]])
                }
                band <- sprintf ("V%d", nb.dwt.levels)
                cf <- c (cf, slot (wt, "V") [[band]])

                cf.sort <- sort(abs(cf))
                sigma <- sqrt(mean(cf.sort[length(cf)/2], cf.sort[(length(cf)/2)+1])/0.6745)

### Essa aproximação do sigma calculada pela mediana é robusta a outliers
### T is the threshold corresponding to the standard deviation of noise.
### By deploying wavelet thresholding, noise in EEG signals can be eliminated. 
                T <- sigma*sqrt(2*log(length(cf)))
                
######################################################################################################
                threshold <-  function(W,T){
                    if (abs(W)>=T){
                        X <-(W/(abs(W)))*(abs(W)-T)
                    }
                    else{
                        X <- 0
                    }
                }
#########################################################################################
###Applies treshold in each level of the wt object
                for (i in 1 : nb.dwt.levels) {
                    level=c()
                    band <- sprintf ("W%d", i)
                    level <- slot (wt, "W") [[band]]
                    
                    for (l in 1:length(level)){
                        slot(wt, "W")[[band]][l,1] <- threshold(level[l],T)
                    }
                }
                
                level=c()
                band <- sprintf ("V%d", nb.dwt.levels)
                level <- slot (wt, "V") [[band]]
                
                for (l in 1:length(level)){
                    slot(wt, "V")[[band]][l,1] <- threshold(level[l],T)
                }

#################################################################
### Zeroing the high frequency bands
                for (i in 1 : cutBand) {
                    level=c()
                    band <- sprintf ("W%d", i)
                    level <- slot (wt, "W") [[band]]
                    
                    for (l in 1:length(level)){
                        slot(wt, "W")[[band]][l,1] <- 0
                    }
                }

                for (i in 1 : cutBand) {
                    level=c()
                    band <- sprintf ("V%d", i)
                    level <- slot (wt, "V") [[band]]
                    
                    for (l in 1:length(level)){
                        slot(wt, "V")[[band]][l,1] <- 0
                    }
                }
#########################################################################################################
### Computes inverse DWT to obtain Xica denoised by channel
                Xica_DWT[ch,] <- idwt(wt)
                
            }#nchannels

            filename <- paste(sprintf("../Sujeito%d/",subj),pastas[p],sprintf("Xica_DWT%d.mat",stim),sep="")
            writeMat (filename, Xica_DWT=Xica_DWT, W=file$W, w0=file$w0, X=file$X)

            cat (sprintf ("\rSubject: %2d   folder: %2d   stimulus: %d", subj, p, stim))
            flush (stdout ())
        }#stim
    }#pasta
}#subject




