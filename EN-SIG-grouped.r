### Computes Elastic Net coefficients for left and right hemispheres separately and for physical
### and psychophysical responses for several alpha values. The best set of coefficients are choosen
### and the predictions of the response are computed. Coefficients and predictions are plotted and
### a paired t-test is computed to validate the existence of a categorization effect.

rm(list=ls())

library(R.matlab)
library(wavelets)
library(glmnet)
library(effsize)
library(equivalence)
library(ANOVAreplication)

source ("parameters_sig.r")
###source("sig-grouped.r")

########################################################################
### Values of alpha(sparcity factor) for glmnet 
alphas <- seq (0, 1, by = 0.02)
### Apply shuffling 0-nothing 1-shuff
sqr.rand=0
#######################################################################
addresses=c("/VOT/Passivo/","/VOT/Ativo/","/Formantes/Passivo/","/Formantes/Ativo/")
nome=c("VOTpass", "VOTact", "Formpass", "Formact")

###beta <- readMat("beta.mat")
###btVOT <- sort(beta$bt[1,subjects],index.return=TRUE)
###btForm <- sort(beta$bt[2,subjects],index.return=TRUE)
phyVOT <- c(1, 34, 74, 125, 200)
###phyVOT <- c(-52, -40, -26, 8, 17)
phyForm <- c(1, 55, 97, 146, 200)
###phyForm <- c(533, 961, 1103, 1297, 1387)#F2-F1
###phyForm <- c(1.703, 2.2779, 2.538, 3.0785, 3.5356)#F2/F1

for (ad in 2:4){

    X <- c ()
    Y <- c ()

    AIC.coef <- data.frame(psy.left=0, psy.right=0, phy.left=0, phy.right=0)

    ## ROI
    idx <- seq (2, (0.35*fs.eeg)+1)

    aux <- length(idx)

    matrix.phy.left <- c()
    matrix.psy.left <- c()
    matrix.phy.right <- c()
    matrix.psy.right <- c()

    aux.left1=c()
    aux.left2=c()
    aux.left3=c()
    aux.left4=c()
    aux.left5=c()

    aux.right1=c()
    aux.right2=c()
    aux.right3=c()
    aux.right4=c()
    aux.right5=c()

###Organiza beta na ordem crescente
    
    if (ad<3){
        ix=btVOT$ix
    }else{
        ix=btForm$ix
    }
    
    for (subj in subjects[ix]){

#########################################################################################################################
### Encontra termos para a  normalização dos dados de cada sujeito nos seus 10 estímulos e nas condições ativa e passiva
### A normalização é feita por aquisição, ou seja, VOT ou Formantes        
        juntoL=c()
        juntoR=c()
        if (ad<=2){
            load (file.path (sprintf(paste("../../Sujeito%d",addresses[1],"sig_grouped_left.dat",sep=""), subj)))
            left <- dwt.coefs
            load (file.path (sprintf(paste("../../Sujeito%d",addresses[1],"sig_grouped_right.dat",sep=""), subj)))
            right <- dwt.coefs
            for (i in 1:nrow(left)){
                juntoL <- c(juntoL, left[i,2:(nb.points+1)])
                juntoR <- c(juntoR, right[i,2:(nb.points+1)])
            }
            load (file.path (sprintf(paste("../../Sujeito%d",addresses[2],"sig_grouped_left.dat",sep=""), subj)))
            left <- dwt.coefs
            load (file.path (sprintf(paste("../../Sujeito%d",addresses[2],"sig_grouped_right.dat",sep=""), subj)))
            right <- dwt.coefs
            for (i in 1:nrow(left)){
                juntoL <- c(juntoL, left[i,2:(nb.points+1)])
                juntoR <- c(juntoR, right[i,2:(nb.points+1)])
            }
        } else{
            load (file.path (sprintf(paste("../../Sujeito%d",addresses[3],"sig_grouped_left.dat",sep=""), subj)))
            left <- dwt.coefs
            load (file.path (sprintf(paste("../../Sujeito%d",addresses[3],"sig_grouped_right.dat",sep=""), subj)))
            right <- dwt.coefs
            for (i in 1:nrow(left)){
                juntoL <- c(juntoL, left[i,2:(nb.points+1)])
                juntoR <- c(juntoR, right[i,2:(nb.points+1)])
            }
            load (file.path (sprintf(paste("../../Sujeito%d",addresses[4],"sig_grouped_left.dat",sep=""), subj)))
            left <- dwt.coefs
            load (file.path (sprintf(paste("../../Sujeito%d",addresses[4],"sig_grouped_right.dat",sep=""), subj)))
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
        
        load (file.path (sprintf(paste("../../Sujeito%d",addresses[ad],"sig_grouped_left.dat",sep=""), subj)))
        left <- dwt.coefs

        load (file.path (sprintf(paste("../../Sujeito%d",addresses[ad],"sig_grouped_right.dat",sep=""), subj)))
        right <- dwt.coefs

        ## ROI selection

        coefs.left <- left [, c (1, idx)]
        coefs.right <- right [, c (1,idx)]

        stimidx <- coefs.left [, 1]
        
        coefs.left <-(coefs.left-medL)/sdL
        coefs.right <- (coefs.right-medR)/sdR

        aux=groups*2
        
        matrix.psy.left <-rbind(matrix.psy.left, coefs.left)
        matrix.psy.right <-rbind(matrix.psy.right, coefs.right)

        aux.left1 <-rbind(aux.left1, coefs.left[1:aux,])
        aux.right1 <-rbind(aux.right1, coefs.right[1:aux,])

        aux.left2 <-rbind(aux.left2, coefs.left[(aux+1):(2*aux),])
        aux.right2 <-rbind(aux.right2, coefs.right[(aux+1):(2*aux),])

        aux.left3 <-rbind(aux.left3, coefs.left[((2*aux)+1):(3*aux),])
        aux.right3 <-rbind(aux.right3, coefs.right[((2*aux)+1):(3*aux),])

        aux.left4 <-rbind(aux.left4, coefs.left[((3*aux)+1):(4*aux),])
        aux.right4 <-rbind(aux.right4, coefs.right[((3*aux)+1):(4*aux),])

        aux.left5 <-rbind(aux.left5, coefs.left[((4*aux)+1):(5*aux),])
        aux.right5 <-rbind(aux.right5, coefs.right[((4*aux)+1):(5*aux),])

    }#subjects

    matrix.phy.left <- rbind(aux.left1, aux.left2, aux.left3, aux.left4, aux.left5)
    matrix.phy.right <- rbind(aux.right1, aux.right2, aux.right3, aux.right4, aux.right5)
###############################################
    matrix.psy.left <-  matrix.phy.left
    matrix.psy.right <-  matrix.phy.right
###############################################

###normal or randomized
    if (sqr.rand==1){
        fname="-Shuff"
        matrix.phy.left <- matrix.phy.left[sample(nrow(matrix.phy.left)),]
        matrix.phy.right <- matrix.phy.right[sample(nrow(matrix.phy.right)),]
        matrix.psy.left <- matrix.psy.left[sample(nrow(matrix.psy.left)),]
        matrix.psy.right <- matrix.psy.right[sample(nrow(matrix.psy.right)),]
    }else{
        fname=""
    }

###############################################################################################
### Repete o beta de cada sujeito pela quantidade de estímulos dele(a)
    if (ad<=2){
###response1 <- rep (btVOT$x, rep (length(stimidx),nb.subjects))
        response1 <- rep(c (0, 0.05, 0.5, 0.95, 1), rep(groups*2*nb.subjects,5))
    }else{
###response1 <- rep (btForm$x, rep (length(stimidx), nb.subjects))
        response1 <- rep(c (0, 0.05, 0.5, 0.95, 1), rep(groups*2*nb.subjects,5))
    }
###############################################################################################        
    psy.EN.left <- list ()
    for (i in seq (1, length (alphas))) {
        psy.EN.left [[i]] <- cv.glmnet (matrix.psy.left [, 2 : ncol (matrix.psy.left)],
                                        response1,
                                        nfolds = nstim,
                                        foldid = rep (seq (1, groups), nstim*nb.subjects),
                                        alpha = alphas [i])
        cat (sprintf ("\raddress: %d   alpha: %.2f  psy.left", ad, alphas [i]))
        flush (stdout ())
    }

    cat ("\n")
    flush (stdout ())

###Model with less cross-validation mean error (cvm)
    psy.EN.left.best <- psy.EN.left [[which.min (sapply (psy.EN.left, function (x) min (x$cvm)))]]

    fit <- glmnet(matrix.psy.left [, 2 : ncol (matrix.psy.left)], response1, lambda=psy.EN.left.best$lambda.min)

    tLL <- fit$nulldev - deviance(fit)
    k <- fit$df
    n <- fit$nobs

    AIC.coef$psy.left <- (-tLL+2*k+2*k*(k+1))/(n-k-1)

### lambda.min is the value of lambda that provides minimum cvm 
    coef.psy.left <- coef (psy.EN.left.best, lambda = "lambda.min")

##################################################################################################
    psy.EN.right <- list ()
    for (i in seq (1, length (alphas))) {
        psy.EN.right [[i]] <- cv.glmnet (matrix.psy.right [, 2 : ncol (matrix.psy.right)],
                                         response1,
                                         nfolds = nstim,
                                         foldid = rep (seq (1, groups), nstim*nb.subjects),
                                         alpha = alphas [i])
        cat (sprintf ("\raddress: %d   alpha: %.2f   psy.right", ad, alphas [i]))
        flush (stdout ())
    }

    cat ("\n")
    flush (stdout ())

    psy.EN.right.best <- psy.EN.right [[which.min (sapply (psy.EN.right, function (x) min (x$cvm)))]]

    fit <- glmnet(matrix.psy.right [, 2 : ncol (matrix.psy.right)], response1, lambda=psy.EN.right.best$lambda.min)

    tLL <- fit$nulldev - deviance(fit)
    k <- fit$df
    n <- fit$nobs

    AIC.coef$psy.right <- (-tLL+2*k+2*k*(k+1))/(n-k-1)

    coef.psy.right <- coef (psy.EN.right.best, lambda = "lambda.min")
##################################################################################################

    if (ad<3){
        response2 <- rep (phyVOT, rep (groups*2*nb.subjects,length(phyVOT)))/200
    }else{
        response2 <- rep (phyForm, rep (groups*2*nb.subjects,length(phyForm)))/200
    }
################################################################################################################
    phy.EN.left <- list ()
    for (i in seq (1, length (alphas))) {
        phy.EN.left [[i]] <- cv.glmnet (matrix.phy.left [, 2 : ncol (matrix.phy.left)],
                                        response2,
                                        nfolds = nstim,
                                        foldid = rep (seq (1, groups), nstim*nb.subjects),
                                        alpha = alphas [i])
        cat (sprintf ("\raddress: %d   alpha: %.2f  phy.left", ad, alphas [i]))
        flush (stdout ())
    }

    cat ("\n")
    flush (stdout ())

    phy.EN.left.best <- phy.EN.left [[which.min (sapply (phy.EN.left, function (x) min (x$cvm)))]]

    fit <- glmnet(matrix.phy.left [, 2 : ncol (matrix.phy.left)], response2, lambda=phy.EN.left.best$lambda.min)

    tLL <- fit$nulldev - deviance(fit)
    k <- fit$df
    n <- fit$nobs

    AIC.coef$phy.left <- (-tLL+2*k+2*k*(k+1))/(n-k-1)

    coef.phy.left <- coef (phy.EN.left.best, lambda = "lambda.min")
####################################################################################################
    phy.EN.right <- list ()
    for (i in seq (1, length (alphas))) {
        phy.EN.right [[i]] <- cv.glmnet (matrix.phy.right [, 2 : ncol (matrix.phy.right)],
                                         response2,
                                         nfolds = nstim,
                                         foldid = rep (seq (1, groups), nstim*nb.subjects),
                                         alpha = alphas [i])
        cat (sprintf ("\raddress: %d   alpha: %.2f    phy.right", ad, alphas [i]))
        flush (stdout ())
    }

    cat ("\n")
    flush (stdout ())

    phy.EN.right.best <- phy.EN.right [[which.min (sapply (phy.EN.right, function (x) min (x$cvm)))]]

    fit <- glmnet(matrix.phy.right [, 2 : ncol (matrix.phy.right)], response2, lambda=phy.EN.right.best$lambda.min)

    tLL <- fit$nulldev - deviance(fit)
    k <- fit$df
    n <- fit$nobs

    AIC.coef$phy.right <- (-tLL+2*k+2*k*(k+1))/(n-k-1)

    coef.phy.right <- coef (phy.EN.right.best, lambda = "lambda.min")
#######################################################################################################
### Predictions

###Predict responses using the model found and the input values, i.e., project input values in model curve and take the y values.
###Input is a matrix where each line is a trial with a scalar response.
    proj.phy.left <- predict (phy.EN.left.best, matrix.phy.left[, 2 : ncol (matrix.phy.left)], s="lambda.min")
    proj.phy.right <- predict (phy.EN.right.best, matrix.phy.right[, 2 : ncol (matrix.phy.right)], s="lambda.min")

    proj.psy.left <- predict (psy.EN.left.best, matrix.psy.left[, 2 : ncol (matrix.psy.left)], s="lambda.min")
    proj.psy.right <- predict (psy.EN.right.best, matrix.psy.right[, 2 : ncol (matrix.psy.right)], s="lambda.min")

#################################################################################################################
    R1=cor(proj.psy.left,response1)
    R2=cor(proj.psy.right,response1)
    R3=cor(proj.phy.left,response2)
    R4=cor(proj.phy.right,response2)
#################################################################################################################
    par(mfrow=c(2,2))

    R=cor(proj.psy.left,response1)
    test=t.test(proj.psy.left, response1, alternartive="two.sided", paired = TRUE, conf.level=0.95)
    plot(response1,proj.psy.left,xlab="Psy values", ylab="Predictions of psy values", main=sprintf("Psy-left, R2=%g",R))
    abline(lm(response1~proj.psy.left))
    
    
###################################################################################################################
    R=cor(proj.psy.right,response1)
    test=t.test(proj.psy.right, response1, alternartive="two.sided", paired = TRUE, conf.level=0.95)
    plot(response1,proj.psy.right,xlab="Psy values", ylab="Predictions of psy values",  main=sprintf("Psy-right, R2=%g",R))
    abline(lm(response1~proj.psy.right))
    
##########################################################################################################
    R=cor(proj.phy.left,response2)
    d=cohen.d(as.numeric(proj.psy.left), response2,pooled=TRUE,paired=TRUE)$estimate
    tost(as.numeric(proj.psy.left), response2,epsilon=abs(d),paired=TRUE, var.equal=TRUE)
    test=t.test(proj.psy.left, response2, alternartive="two.sided", paired = TRUE, conf.level=0.95)
    plot(response2,proj.phy.left,xlab="Stimulus characteristic", ylab="Predictions of stimulus values", main=sprintf("Phy-left, R2=%g",R))
    abline(lm(response2~proj.phy.left))
    
###################################################################################################################
    R=cor(proj.phy.right,response2)
    test=t.test(proj.phy.right, response1, alternartive="two.sided", paired = TRUE, conf.level=0.95)
    plot(response2,proj.phy.right,xlab="Stimulus characteristic", ylab="Predictions of stimulus values", main=sprintf("Phy-right, R2=%g",R))
    abline(lm(response2~proj.phy.right))
    
    dev.copy2pdf(file = file.path (sprintf("Predictions-%s%s-psyold.pdf", c("VOTpass","VOTact","FormPass","FormAct")[ad],fname)))
###############################################################################################################################
    graphics.off()

    ry <- matrix(response1, nrow = 2*groups*nb.subjects)
    rx <- matrix(response2, nrow = 2*groups*nb.subjects)
### Predictions plots right

    mb.r <- matrix(proj.psy.right, nrow =2* groups*nb.subjects)
    mp.r <- matrix(proj.phy.right, nrow =2* groups*nb.subjects)

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

### Statistic analysis

    test <- t.test ((mp.r[,2]-mp.r[,1]), (mb.r[,2]-mb.r[,1]), paired = TRUE, alternative="greater")

### Fisher analysis
    m1=mean(mp.r[,1])
    m2=mean(mp.r[,2])
    S1=sum((mp.r[,1]-m1)*(mp.r[,1]-m1))/(length(mp.r[,1])-1)
    S2=sum((mp.r[,2]-m2)*(mp.r[,2]-m2))/(length(mp.r[,2])-1)
    fphy=((m2-m1)^2)/sum(S1,S2)

    m1=mean(mb.r[,1])
    m2=mean(mb.r[,2])
    S1=sum((mb.r[,1]-m1)*(mb.r[,1]-m1))/(length(mb.r[,1])-1)
    S2=sum((mb.r[,2]-m2)*(mb.r[,2]-m2))/(length(mb.r[,2])-1)
    fpsy=((m2-m1)^2)/sum(S1,S2)
    

    dev.new()
    plot (0, 0, type = "n",xlim=c(mi,ma),ylim=c(mi,ma), las = 1,
          xlab = "physical predictions (right)", ylab = "psychophycal predictions (right)", main=sprintf("Proj-right, Fisher psy=%g, Fisher phy=%g, t.test(dist2-1. psy>phy)=%g",fpsy, fphy, test$p.value),cex.main=0.9)
    grid(col = "lightgray", lty = "dotted")
    for (s in seq (1, nrow(mb.r)))
        points (mp.r[s,], mb.r[s,], pch = 19, col = seq(1,5))

    lines (colMeans (mp.r), colMeans (mb.r), lwd=2)
    points (colMeans (mp.r), colMeans (mb.r), cex = 2)

    lines (colMeans (rx), colMeans (ry), lwd=1, col="blue")
    points (colMeans (rx), colMeans (ry), pch = 8, cex=1)
    
    lines((d.phy2$x),(d.phy2$y/(4*max(d.phy4$y)))-adjust, col="red", lwd=2)
    lines((d.phy4$x),(d.phy4$y/(4*max(d.phy4$y)))-adjust, col="blue",lwd=2)

    lines((d.psy2$y/(4*max(d.psy2$y)))-adjust,(d.psy2$x), col="red", lwd=2)
    lines((d.psy4$y/(4*max(d.psy4$y)))-adjust,(d.psy4$x), col="blue", lwd=2)

    dev.copy2pdf(file = file.path (sprintf("Proj-all-right%s.pdf",nome[ad])))
################################################################################################################
###############################################################################################################################
    graphics.off()
### Predictions plots left

    mb.l <- matrix(proj.psy.left, nrow =2* groups*nb.subjects)
    mp.l <- matrix(proj.phy.left, nrow =2* groups*nb.subjects)

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

### Statistic analysis

    test <- t.test ((mp.l[,2]-mp.l[,1]), (mb.l[,2]-mb.l[,1]), paired = TRUE, alternative="greater")

### Fisher analysis
    m1=mean(mp.l[,1])
    m2=mean(mp.l[,2])
    S1=sum((mp.l[,1]-m1)*(mp.l[,1]-m1))/(length(mp.l[,1])-1)
    S2=sum((mp.l[,2]-m2)*(mp.l[,2]-m2))/(length(mp.l[,2])-1)
    fphy=((m2-m1)^2)/sum(S1,S2)

    m1=mean(mb.l[,1])
    m2=mean(mb.l[,2])
    S1=sum((mb.l[,1]-m1)*(mb.l[,1]-m1))/(length(mb.l[,1])-1)
    S2=sum((mb.l[,2]-m2)*(mb.l[,2]-m2))/(length(mb.l[,2])-1)
    fpsy=((m2-m1)^2)/sum(S1,S2)
    

    dev.new()
    plot (0, 0, type = "n",xlim = c(mi,ma), ylim = c (mi,ma), las = 1,
          xlab = "physical predictions (left)", ylab = "psychophycal predictions (left)", main=sprintf("Proj-left, Fisher psy=%g, Fisher phy=%g, t.test(dist2-1. psy>phy)=%g",fpsy, fphy, test$p.value),cex.main=0.9)
    grid(col = "lightgray", lty = "dotted")
    for (s in seq (1, nrow(mb.l)))
        points (mp.l[s,], mb.l[s,], pch = 19, col = seq(1,5))

    lines (colMeans (mp.l), colMeans (mb.l), lwd=2)
    points (colMeans (mp.l), colMeans (mb.l), cex = 2)

    lines (colMeans (rx), colMeans (ry), lwd=2, col="blue")
    points (colMeans (rx), colMeans (ry), pch = 8, cex=1)
    
    lines((d.phy2$x),(d.phy2$y/(4*max(d.phy4$y)))-adjust, col="red", lwd=2)
    lines((d.phy4$x),(d.phy4$y/(4*max(d.phy4$y)))-adjust, col="blue",lwd=2)

    lines((d.psy2$y/(4*max(d.psy2$y)))-adjust,(d.psy2$x), col="red", lwd=2)
    lines((d.psy4$y/(4*max(d.psy4$y)))-adjust,(d.psy4$x), col="blue", lwd=2)
    
    dev.copy2pdf(file = file.path (sprintf("Proj-all-left%s.pdf",nome[ad])))
#######################################################################################################################    
    
    coef.phy <- c(coef.phy.left[2:length(coef.phy.left)], coef.phy.right[2:length(coef.phy.right)])
    coef.psy <- c(coef.psy.left[2:length(coef.psy.left)], coef.psy.right[2:length(coef.psy.right)])

    save(AIC.coef,phy.EN.right.best, psy.EN.right.best, phy.EN.left.best, psy.EN.left.best, coef.phy.left, coef.phy.right, coef.psy.left, coef.psy.right, file=file.path(sprintf("SIG_all_%s%s-psyold.dat",c("VOTpass","VOTact","FormPass","FormAct")[ad],fname)))
    
##########################################################################################################    
    filename=sprintf("SIG-matrices-all-%s%s-psyold.mat",c("VOTpass","VOTact","FormPass","FormAct")[ad],fname)
    writeMat(filename,R1=R1, R2=R2, R3=R3, R4=R4,coef.phy=coef.phy, coef.psy=coef.psy,proj.psy.left=proj.psy.left, proj.psy.right=proj.psy.right, proj.phy.left=proj.phy.left, proj.phy.right=proj.phy.right, AIC.coef=AIC.coef)

}#addresses

