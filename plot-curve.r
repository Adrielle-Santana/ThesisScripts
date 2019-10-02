plot.curve <- function (file.name) {

    library(robust)

###############################################################################
### If Expyriment version is below 0.9.0, change skip parameter value to 9 ###
###############################################################################
    dat <- read.csv (file.name, skip = 10, na.strings = "None")
    dat$tecla[which(dat$tecla==113)] <- NA
    
    ### Remove lines with NA's
    dat <- dat[complete.cases(dat),]
### If the sequency pa-pe is flipped (st=1), swap the answers LCTRL with LALT to equal to the response to st=0      
    for (t in which(dat$st==1)){
        if (dat$tecla[t]==308) dat$tecla[t] <- 306 #LCTRL
        else if (dat$tecla[t]==306) dat$tecla[t] <- 308 #RCTRL
    }
### Creates a new column replacing the numbers of the keys by 0's and 1's
### Assign 1 to 308 and 0 to 306
    ### Trocar 0 com 1 quando usar table keyboard.
    dat$resp <- c (0, 1) [factor (dat$tecla)]
    dat$stim <- abs(dat$stim)

    fm <- glmRob (resp ~ stim, family = binomial, data=dat, method="misclass") 

### Plot results

    dev.new(width=11, height=9)
    plot (dat$stim, dat$resp, xlab= "Stimulus number", ylab= "Response: ba-pa", main=sprintf("Psicometric curve of subject: %s", dat$subject_id[1]))
    lines (dat$stim[order(dat$stim)], predict (fm, type = "response")[order(dat$stim)], col = "red")

    k <- predict (fm, type = "response")[order(dat$stim)]
    
    # Position of both ambiguous stimuli in the 200 continuum 
    pa.amb <- dat$stim[order(dat$stim)][which(round(k,2)==0.05)[1]]
    amb <- dat$stim[order(dat$stim)][which(round(k,1)==0.5)[1]]
    pe.amb <- dat$stim[order(dat$stim)][which(round(k,2)==0.95)[1]]
    
    points(pa.amb, k[which(round(k,2)==0.05)[1]], pch=19)
    points(amb, k[which(round(k,1)==0.5)[1]], pch=19)
    points(pe.amb, k[which(round(k,2)==0.95)[1]], pch=19)

    dev.copy2pdf(file =  sprintf ("%s.pdf", file.name))

    #graphics.off()

    rt.mean <- mean(dat$rt)
    
    save (pa.amb, amb, pe.amb, fm, rt.mean, file = file.path (sprintf("%s.dat", file.name)))
}

