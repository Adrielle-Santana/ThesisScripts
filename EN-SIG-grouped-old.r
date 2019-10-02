### Computes Elastic Net coefficients for left and right hemispheres separately and for physical
### and psychophysical responses for several alpha values. The best set of coefficients are choosen
### and the predictions of the response are computed. Coefficients and predictions are plotted and
### a paired t-test is computed to validate the existence of a categorization effect.

rm(list=ls())

library(R.matlab)
library(wavelets)
library(glmnet)

source ("parameters_sig.r")
source("sig-grouped.r")

########################################################################
### Values of alpha(sparcity factor) for glmnet 
alphas <- seq (0, 1, by = 0.02)
### Apply sum squares prenormalization or shuffling 0-nothing 1-ssqrt 2-shuff
sqr.rand=0
#######################################################################
addresses=c("/VOT/Passivo/","/VOT/Ativo/","/Formantes/Passivo/","/Formantes/Ativo/")

for (ad in 1:4){

    X <- c ()
    Y <- c ()

    mps.l <- matrix(rep(0, 5*nb.subjects), nrow=nb.subjects)
    mps.r <- matrix(rep(0, 5*nb.subjects), nrow=nb.subjects)
    mbs.l <- matrix(rep(0, 5*nb.subjects), nrow=nb.subjects)
    mbs.r <- matrix(rep(0, 5*nb.subjects), nrow=nb.subjects)

    AIC.coef <- data.frame(psy.left=rep(0,nb.subjects), psy.right=rep(0,nb.subjects), phy.left=rep(0,nb.subjects), phy.right=rep(0,nb.subjects))

    ## ROI
    idx <- seq (2, (0.35*fs.eeg)+1)

    aux <- length(idx)

    coefs.phy <- matrix(rep(0,nb.subjects*aux*2),nrow=nb.subjects)
    coefs.psy <- matrix(rep(0,nb.subjects*aux*2),nrow=nb.subjects)
    count <- 1

    for (subj in subjects){

#########################################################################################################################
### Encontra termos para a  normalização dos dados de cada sujeito nos seus 10 estímulos e nas condições ativa e passiva
        juntoL=c()
        juntoR=c()
        if (ad<=2){
            load (file.path (sprintf(paste("../Sujeito%d",addresses[1],"sig_grouped_left.dat",sep=""), subj)))
            left <- dwt.coefs
            load (file.path (sprintf(paste("../Sujeito%d",addresses[1],"sig_grouped_right.dat",sep=""), subj)))
            right <- dwt.coefs
            for (i in 1:nrow(left)){
                juntoL <- c(juntoL, left[i,2:(nb.points+1)])
                juntoR <- c(juntoR, right[i,2:(nb.points+1)])
            }
            load (file.path (sprintf(paste("../Sujeito%d",addresses[2],"sig_grouped_left.dat",sep=""), subj)))
            left <- dwt.coefs
            load (file.path (sprintf(paste("../Sujeito%d",addresses[2],"sig_grouped_right.dat",sep=""), subj)))
            right <- dwt.coefs
            for (i in 1:nrow(left)){
                juntoL <- c(juntoL, left[i,2:(nb.points+1)])
                juntoR <- c(juntoR, right[i,2:(nb.points+1)])
            }
        } else{
            load (file.path (sprintf(paste("../Sujeito%d",addresses[3],"sig_grouped_left.dat",sep=""), subj)))
            left <- dwt.coefs
            load (file.path (sprintf(paste("../Sujeito%d",addresses[3],"sig_grouped_right.dat",sep=""), subj)))
            right <- dwt.coefs
            for (i in 1:nrow(left)){
                juntoL <- c(juntoL, left[i,2:(nb.points+1)])
                juntoR <- c(juntoR, right[i,2:(nb.points+1)])
            }
            load (file.path (sprintf(paste("../Sujeito%d",addresses[4],"sig_grouped_left.dat",sep=""), subj)))
            left <- dwt.coefs
            load (file.path (sprintf(paste("../Sujeito%d",addresses[4],"sig_grouped_right.dat",sep=""), subj)))
            right <- dwt.coefs
            for (i in 1:nrow(left)){
                juntoL <- c(juntoL, left[i,2:(nb.points+1)])
                juntoR <- c(juntoR, right[i,2:(nb.points+1)])
            }
        }
        medL=mean(juntoL)
        medR=mean(juntoR)
        sdL=sd(juntoL)
        sdR=sd(juntoR)
#########################################################################################################################
        
        subj.dir <- file.path (sprintf(paste("../Sujeito%d",addresses[ad],"EN",sep=""), subj))
        dir.create (subj.dir, showWarnings = FALSE)

        load (file.path (sprintf(paste("../Sujeito%d",addresses[ad],"sig_grouped_left.dat",sep=""), subj)))
        left <- dwt.coefs

        load (file.path (sprintf(paste("../Sujeito%d",addresses[ad],"sig_grouped_right.dat",sep=""), subj)))
        right <- dwt.coefs

        ## ROI selection

        coefs.left <- left [, c (1, idx)]
        coefs.right <- right [, c (1,idx)]

        stimidx <- coefs.left [, 1]
        
        coefs.left <-(coefs.left-medL)/sdL
        coefs.right <- (coefs.right-medR)/sdR

###normalize by the root sum squares (sd) or randomize
        if (sqr.rand==1){
            fname="-sqrt"
            coefs.left <- (coefs.left[, 2 : ncol (coefs.left)]-mean(coefs.left[, 2 : ncol (coefs.left)]))/(sqrt(sum(coefs.left[, 2 : ncol (coefs.left)]^2)))
            coefs.right <- (coefs.right[, 2 : ncol (coefs.right)]-mean(coefs.right[, 2 : ncol (coefs.right)]))/(sqrt(sum(coefs.right[, 2 : ncol (coefs.right)]^2)))
        } else if (sqr.rand==2){
            fname="-shuff"
            coefs.left <- coefs.left[sample(nrow(coefs.left)),]
            coefs.right <- coefs.right[sample(nrow(coefs.right)),]
        }else{
            fname=""
        }
                
###############################################################################################
response1 <- rep (c (0, 0.05, 0.5, 0.95, 1), rep (2, 5)) [stimidx]
###############################################################################################        
        psy.EN.left <- list ()
        for (i in seq (1, length (alphas))) {
            psy.EN.left [[i]] <- cv.glmnet (coefs.left [, 2 : ncol (coefs.left)],
                                            response1,
                                            nfolds = nstim,
                                            foldid = rep (seq (1, groups), nstim),
                                            alpha = alphas [i])
            cat (sprintf ("\raddress: %d   subject: %02d   alpha: %.2f  psy.left", ad, subj, alphas [i]))
            flush (stdout ())
        }

        cat ("\n")
        flush (stdout ())

        ###Model with less cross-validation mean error (cvm)
        psy.EN.left.best <- psy.EN.left [[which.min (sapply (psy.EN.left, function (x) min (x$cvm)))]]

        fit <- glmnet(coefs.left [, 2 : ncol (coefs.left)], response1, lambda=psy.EN.left.best$lambda.min)

        tLL <- fit$nulldev - deviance(fit)
        k <- fit$df
        n <- fit$nobs

        AIC.coef$psy.left[count] <- (-tLL+2*k+2*k*(k+1))/(n-k-1)

        ### lambda.min is the value of lambda that provides minimum cvm 
        coef.psy.left <- coef (psy.EN.left.best, lambda = "lambda.min")

##################################################################################################
        psy.EN.right <- list ()
        for (i in seq (1, length (alphas))) {
            psy.EN.right [[i]] <- cv.glmnet (coefs.right [, 2 : ncol (coefs.right)],
                                             response1,
                                             nfolds = nstim,
                                             foldid = rep (seq (1, groups), nstim),
                                             alpha = alphas [i])
            cat (sprintf ("\raddress: %d   subject: %02d   alpha: %.2f  psy.right", ad, subj, alphas [i]))
            flush (stdout ())
        }

        cat ("\n")
        flush (stdout ())

        psy.EN.right.best <- psy.EN.right [[which.min (sapply (psy.EN.right, function (x) min (x$cvm)))]]

        fit <- glmnet(coefs.right [, 2 : ncol (coefs.right)], response1, lambda=psy.EN.right.best$lambda.min)

        tLL <- fit$nulldev - deviance(fit)
        k <- fit$df
        n <- fit$nobs

        AIC.coef$psy.right[count] <- (-tLL+2*k+2*k*(k+1))/(n-k-1)

        coef.psy.right <- coef (psy.EN.right.best, lambda = "lambda.min")
##################################################################################################

        if (ad<3){
            IDnome="idVOT.csv"
        }
        else{
            IDnome="idFormantes.csv"
        }
        
        idphon <- read.csv (IDnome,
                            sep = ",", skip = 1, header = FALSE,
                            col.names = c ("subj", "stim2", "stim3","stim4"))

        aux <- idphon[subj,c(2,3,4)]
        response2 <- rep (c (1,aux[[1]],aux[[2]],aux[[3]], 200)/200, rep (2, 5)) [stimidx]

################################################################################################################
        phy.EN.left <- list ()
        for (i in seq (1, length (alphas))) {
            phy.EN.left [[i]] <- cv.glmnet (coefs.left [, 2 : ncol (coefs.left)],
                                            response2,
                                            nfolds = nstim,
                                            foldid = rep (seq (1, groups), nstim),
                                            alpha = alphas [i])
            cat (sprintf ("\raddress: %d   subject: %02d   alpha: %.2f  phy.left", ad, subj, alphas [i]))
            flush (stdout ())
        }

        cat ("\n")
        flush (stdout ())

        phy.EN.left.best <- phy.EN.left [[which.min (sapply (phy.EN.left, function (x) min (x$cvm)))]]

        fit <- glmnet(coefs.left [, 2 : ncol (coefs.left)], response2, lambda=phy.EN.left.best$lambda.min)

        tLL <- fit$nulldev - deviance(fit)
        k <- fit$df
        n <- fit$nobs

        AIC.coef$phy.left[count] <- (-tLL+2*k+2*k*(k+1))/(n-k-1)

        coef.phy.left <- coef (phy.EN.left.best, lambda = "lambda.min")
####################################################################################################
        phy.EN.right <- list ()
        for (i in seq (1, length (alphas))) {
            phy.EN.right [[i]] <- cv.glmnet (coefs.right [, 2 : ncol (coefs.right)],
                                             response2,
                                             nfolds = nstim,
                                             foldid = rep (seq (1, groups), nstim),
                                             alpha = alphas [i])
            cat (sprintf ("\raddress: %d   subject: %02d   alpha: %.2f  phy.right", ad, subj, alphas [i]))
            flush (stdout ())
        }

        cat ("\n")
        flush (stdout ())

        phy.EN.right.best <- phy.EN.right [[which.min (sapply (phy.EN.right, function (x) min (x$cvm)))]]

        fit <- glmnet(coefs.right [, 2 : ncol (coefs.right)], response2, lambda=phy.EN.right.best$lambda.min)

        tLL <- fit$nulldev - deviance(fit)
        k <- fit$df
        n <- fit$nobs

        AIC.coef$phy.right[count] <- (-tLL+2*k+2*k*(k+1))/(n-k-1)

        coef.phy.right <- coef (phy.EN.right.best, lambda = "lambda.min")
#######################################################################################################
### Predictions

###Predict responses using the model found and the input values, i.e., project input values in model curve and take the y values.
###Input is a matrix where each line is a trial with a scalar response.
        proj.phy.left <- predict (phy.EN.left.best, coefs.left[, 2 : ncol (coefs.left)], s="lambda.min")
        proj.phy.right <- predict (phy.EN.right.best, coefs.right[, 2 : ncol (coefs.right)], s="lambda.min")

        proj.psy.left <- predict (psy.EN.left.best, coefs.left[, 2 : ncol (coefs.left)], s="lambda.min")
        proj.psy.right <- predict (psy.EN.right.best, coefs.right[, 2 : ncol (coefs.right)], s="lambda.min")

######################################################################################################
### Original psychometric curve
        ry <- matrix(response1, nrow = 2*groups)
        rx <- matrix(response2, nrow = 2*groups)
### Predictions plots left

        mb.l <- matrix(proj.psy.left, nrow = 2*groups)
        mp.l <- matrix(proj.phy.left, nrow = 2*groups)

        mi <- min(mp.l, mb.l)
        ma <- max(mp.l, mb.l)
        adjust <- round(abs(mi),2)+0.05

        if  (mi>0){
            adjust <- 0.04
            mi <- 0
        }
        
        if (ma<1) ma <- 1

        d.phy2 <- density(mp.l[,2])
        d.phy4 <- density(mp.l[,4])
        d.psy2 <- density(mb.l[,2])
        d.psy4 <- density(mb.l[,4])

        dev.new()
        plot (0, 0, type = "n", xlim = c(mi,ma), ylim = c (mi,ma), las = 1,
              xlab = "physical predictions (left)", ylab = "psychophycal predictions (left)")
        for (s in seq (1, nrow(mb.l)))
            points (mp.l[s,], mb.l[s,], pch = 19, col = seq(1,5))

        lines (colMeans (mp.l), colMeans (mb.l), lwd=2)
        points (colMeans (mp.l), colMeans (mb.l), cex = 2)

        lines (colMeans (rx), colMeans (ry), lwd=1, col="blue")
        points (colMeans (rx), colMeans (ry), pch = 8, cex=1)
        
        abline(0,1,col="gray")

        lines((d.phy2$x),(d.phy2$y/(3*max(d.phy2$y)))-adjust, col="red", lwd=2)
        lines((d.phy4$x),(d.phy4$y/(3*max(d.phy4$y)))-adjust, col="blue",lwd=2)

        lines((d.psy2$y/(3*max(d.psy2$y)))-adjust,(d.psy2$x), col="red", lwd=2)
        lines((d.psy4$y/(3*max(d.psy4$y)))-adjust,(d.psy4$x), col="blue", lwd=2)

        ### Statistic analysis

        test <- t.test ((mb.l[,4]-mb.l[,2]), (mp.l[,4]-mp.l[,2]), paired = TRUE, alternative="greater")

        text(0.3,1,labels=sprintf("p-value of distances psy x phy: %g", test$p.value), cex=0.8)

### Fisher analysis
        m2=mean(mp.l[,2])
        m4=mean(mp.l[,4])
        S2=sum((mp.l[,2]-m2)*(mp.l[,2]-m2))/(length(mp.l[,2])-1)
        S4=sum((mp.l[,4]-m4)*(mp.l[,4]-m4))/(length(mp.l[,4])-1)
        fphy=((m4-m2)^2)/sum(S2,S4)

        m2=mean(mb.l[,2])
        m4=mean(mb.l[,4])
        S2=sum((mb.l[,2]-m2)*(mb.l[,2]-m2))/(length(mb.l[,2])-1)
        S4=sum((mb.l[,4]-m4)*(mb.l[,4]-m4))/(length(mb.l[,4])-1)
        fpsy=((m4-m2)^2)/sum(S2,S4)
        
        text(0.3,0.95,labels=sprintf("Fisher psy: %g   Fisher phy:%g", fpsy, fphy), cex=0.8)        

        dev.copy2pdf(file = file.path (subj.dir, sprintf("SIG-Proj-left%s.pdf",fname)))

        graphics.off()
########################################################################################################
### Predictions plots right

        mb.r <- matrix(proj.psy.right, nrow =2* groups)
        mp.r <- matrix(proj.phy.right, nrow =2* groups)

        mi <- min(mp.r, mb.r)
        ma <- max(mp.r, mb.r)
        adjust <- round(abs(mi),2)+0.05

        if  (mi>0){
            adjust <- 0.04
            mi <- 0
        }
        
        if (ma<1) ma <- 1

        d.phy2 <- density(mp.r[,2])
        d.phy4 <- density(mp.r[,4])
        d.psy2 <- density(mb.r[,2])
        d.psy4 <- density(mb.r[,4])

        dev.new()
        plot (0, 0, type = "n", xlim = c(mi,ma), ylim = c (mi,ma), las = 1,
              xlab = "physical predictions (right)", ylab = "psychophycal predictions (right)")
        grid(col = "lightgray", lty = "dotted")
        for (s in seq (1, nrow(mb.l)))
            points (mp.r[s,], mb.r[s,], pch = 19, col = seq(1,5))

        lines (colMeans (mp.r), colMeans (mb.r), lwd=2)
        points (colMeans (mp.r), colMeans (mb.r), cex = 2)

        lines (colMeans (rx), colMeans (ry), lwd=1, col="blue")
        points (colMeans (rx), colMeans (ry), pch = 8, cex=1)
        
        abline(0,1,col="gray")

        lines((d.phy2$x),(d.phy2$y/(3*max(d.phy4$y)))-adjust, col="red", lwd=2)
        lines((d.phy4$x),(d.phy4$y/(3*max(d.phy4$y)))-adjust, col="blue",lwd=2)

        lines((d.psy2$y/(3*max(d.psy2$y)))-adjust,(d.psy2$x), col="red", lwd=2)
        lines((d.psy4$y/(3*max(d.psy4$y)))-adjust,(d.psy4$x), col="blue", lwd=2)

### Statistic analysis

        test <- t.test ((mb.r[,4]-mb.r[,2]), (mp.r[,4]-mp.r[,2]), paired = TRUE, alternative="greater")
        
        text(0.3,1,labels=sprintf("p-value of distances psy x phy: %g", test$p.value),cex=0.8)

### Fisher analysis
        m2=mean(mp.r[,2])
        m4=mean(mp.r[,4])
        S2=sum((mp.r[,2]-m2)*(mp.r[,2]-m2))/(length(mp.r[,2])-1)
        S4=sum((mp.r[,4]-m4)*(mp.r[,4]-m4))/(length(mp.r[,4])-1)
        fphy=((m4-m2)^2)/sum(S2,S4)

        m2=mean(mb.r[,2])
        m4=mean(mb.r[,4])
        S2=sum((mb.r[,2]-m2)*(mb.r[,2]-m2))/(length(mb.r[,2])-1)
        S4=sum((mb.r[,4]-m4)*(mb.r[,4]-m4))/(length(mb.r[,4])-1)
        fpsy=((m4-m2)^2)/sum(S2,S4)
        
        text(0.3,0.95,labels=sprintf("Fisher psy: %g   Fisher phy:%g", fpsy, fphy), cex=0.8)

        dev.copy2pdf(file = file.path (subj.dir, sprintf("SIG-Proj-right%s.pdf",fname)))

        graphics.off()
##########################################################################################################
        coef.phy <- c(coef.phy.left[2:length(coef.phy.left)], coef.phy.right[2:length(coef.phy.right)])
        coef.psy <- c(coef.psy.left[2:length(coef.psy.left)], coef.psy.right[2:length(coef.psy.right)])

        save(phy.EN.right.best, psy.EN.right.best, phy.EN.left.best, psy.EN.left.best, coef.phy.left, coef.phy.right, coef.psy.left, coef.psy.right, mb.l, mb.r, mp.l, mp.r, file=file.path(subj.dir, sprintf("SIG_subject%02d%s.dat",subj,fname)))
        
##########################################################################################################    

        coefs.phy[count,] <- coef.phy
        coefs.psy[count,] <- coef.psy

        mps.l[count,] <- colMeans(mp.l)
        mps.r[count,] <- colMeans(mp.r)
        mbs.l[count,] <- colMeans(mb.l)
        mbs.r[count,] <- colMeans(mb.r)

        count <- count+1
        
    }#subjects

    save(AIC.coef, coefs.phy, coefs.psy, mps.l, mps.r, mbs.l, mbs.r, file = file.path(sprintf("Todos/SIG-matrices-all-%s%s-%d-groups.dat",c("VOTpass","VOTact","FormPass","FormAct")[ad],fname,groups)))

filename=sprintf("Todos/SIG-matrices-all-%s%s-%d-groups.mat",c("VOTpass","VOTact","FormPass","FormAct")[ad],fname,groups)
writeMat(filename,coefs.phy=coefs.phy, coefs.psy=coefs.psy, AIC.coef=AIC.coef)

###########################################################################################################################
### Plot predictions of everyone

    d.phy2 <- density(mps.l[,2])
    d.phy4 <- density(mps.l[,4])
    d.psy2 <- density(mbs.l[,2])
    d.psy4 <- density(mbs.l[,4])

    N <- length(mps.l[,2])

    plot(mps.l, mbs.l, col=c("black", "red", "yellow","green3","blue")[rep(c(1,2,3,4,5),rep(nb.subjects,5))], ylim=c(0,1), xlim=c(0,1), pch=19, xlab="Mean predictions phisical response left", ylab="Mean predictions psychophysical response left", main="All Subjects Mean Predictions")
    lines(colMeans(mps.l),colMeans(mbs.l))
    points(colMeans(mps.l),colMeans(mbs.l), col="magenta3", pch=17)

    abline(0,1,col="gray")

    lines((d.phy2$x),(d.phy2$y/(2*max(d.phy4$y)))-adjust, col="red", lwd=2)
    lines((d.phy4$x),(d.phy4$y/(2*max(d.phy4$y)))-adjust, col="green3",lwd=2)

    lines((d.psy2$y/(2*max(d.psy2$y)))-adjust,(d.psy2$x), col="red", lwd=2)
    lines((d.psy4$y/(2*max(d.psy4$y)))-adjust,(d.psy4$x), col="green3", lwd=2)

    ### Statistic analysis

    test <- t.test ((mbs.l[,4]-mbs.l[,2]), (mps.l[,4]-mps.l[,2]), paired = TRUE, alternative="greater")

    text(0.3,1,labels=sprintf("p-value of distances psy x phy: %g", test$p.value), cex=0.8)

    ### Pooled variance

    s_p = sqrt(((N - 1) * var(mbs.l[,4]-mbs.l[,2]) +(N - 1) * var(mps.l[,4]-mps.l[,2]) / (N + N - 2)))

### Power

    power.t.test(n = N, delta = abs(mean(c((mbs.l[,4]-mbs.l[,2]), (mps.l[,4]-mps.l[,2])))), sd = s_p,
                 sig.level = 0.05,type = "paired", alternative = "one.sided")

### Fisher analysis
        m2=mean(mps.l[,2])
        m4=mean(mps.l[,4])
        S2=sum((mps.l[,2]-m2)*(mps.l[,2]-m2))/(length(mps.l[,2])-1)
        S4=sum((mps.l[,4]-m4)*(mps.l[,4]-m4))/(length(mps.l[,4])-1)
        fphy=((m4-m2)^2)/sum(S2,S4)

        m2=mean(mbs.l[,2])
        m4=mean(mbs.l[,4])
        S2=sum((mbs.l[,2]-m2)*(mbs.l[,2]-m2))/(length(mbs.l[,2])-1)
        S4=sum((mbs.l[,4]-m4)*(mbs.l[,4]-m4))/(length(mbs.l[,4])-1)
        fpsy=((m4-m2)^2)/sum(S2,S4)
        
        text(0.3,0.95,labels=sprintf("Fisher psy: %g   Fisher phy:%g", fpsy, fphy), cex=0.8)    

    dev.copy2pdf(file = file.path (sprintf("Todos/SIG-Predictions-all-%s-left%s.pdf", c("VOTpass","VOTact","FormPass","FormAct")[ad],fname)))

    graphics.off()
###############################################################################################################
    d.phy2 <- density(mps.r[,2])
    d.phy4 <- density(mps.r[,4])
    d.psy2 <- density(mbs.r[,2])
    d.psy4 <- density(mbs.r[,4])

    plot(mps.r, mbs.r, col=c("black", "red", "yellow", "green3","blue")[rep(c(1,2,3,4,5),rep(nb.subjects,5))], ylim=c(0,1), xlim=c(0,1), pch=19, xlab="Mean predictions phisical response right", ylab="Mean predictions psychophysical response right", main="All Subjects Mean Predictions")
    lines(colMeans(mps.r),colMeans(mbs.r))
    points(colMeans(mps.r),colMeans(mbs.r), col="magenta3", pch=17)

    abline(0,1,col="gray")

    lines((d.phy2$x),(d.phy2$y/(2*max(d.phy4$y)))-adjust, col="red", lwd=2)
    lines((d.phy4$x),(d.phy4$y/(2*max(d.phy4$y)))-adjust, col="green3",lwd=2)

    lines((d.psy2$y/(2*max(d.psy2$y)))-adjust,(d.psy2$x), col="red", lwd=2)
    lines((d.psy4$y/(2*max(d.psy4$y)))-adjust,(d.psy4$x), col="green3", lwd=2)

    ### Statistic analysis

    test <- t.test ((mbs.r[,4]-mbs.r[,2]), (mps.r[,4]-mps.r[,2]), paired = TRUE, alternative="greater")

    text(0.3,1,labels=sprintf("p-value of distances psy x phy: %g", test$p.value), cex=0.8)

    ### Pooled variance

    s_p = sqrt(((N - 1) * var(mbs.r[,4]-mbs.r[,2]) +(N - 1) * var(mps.r[,4]-mps.r[,2]) / (N + N - 2)))

### Power

    power.t.test(n = N, delta = abs(mean(c((mbs.r[,4]-mbs.r[,2]), (mps.r[,4]-mps.r[,2])))), sd = s_p,
                 sig.level = 0.05,type = "paired", alternative = "one.sided")


### Fisher analysis
        m2=mean(mps.r[,2])
        m4=mean(mps.r[,4])
        S2=sum((mps.r[,2]-m2)*(mps.r[,2]-m2))/(length(mps.r[,2])-1)
        S4=sum((mps.r[,4]-m4)*(mps.r[,4]-m4))/(length(mps.r[,4])-1)
        fphy=((m4-m2)^2)/sum(S2,S4)

        m2=mean(mbs.r[,2])
        m4=mean(mbs.r[,4])
        S2=sum((mbs.r[,2]-m2)*(mbs.r[,2]-m2))/(length(mbs.r[,2])-1)
        S4=sum((mbs.r[,4]-m4)*(mbs.r[,4]-m4))/(length(mbs.r[,4])-1)
        fpsy=((m4-m2)^2)/sum(S2,S4)
        
        text(0.3,0.95,labels=sprintf("Fisher psy: %g   Fisher phy:%g", fpsy, fphy), cex=0.8)    

    dev.copy2pdf(file = file.path (sprintf("Todos/SIG-Predictions-all-%s-right%s.pdf", c("VOTpass","VOTact","FormPass","FormAct")[ad],fname)))

}#addresses

