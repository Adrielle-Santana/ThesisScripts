### Computes Elastic Net coefficients for left and right hemispheres separately and for physical
### and psychophysical responses for several alpha values. The best set of coefficients are choosen
### and the predictions of the response are computed. Coefficients and predictions are plotted and
### a paired t-test is computed to validate the existence of a categorization effect.

rm(list=ls())

library(R.matlab)
library(wavelets)
library(glmnet)

source ("parameters_sig.r")
###source("sig-grouped.r")

########################################################################
### Values of alpha(sparcity factor) for glmnet 
alphas <- seq (0, 1, by = 0.02)

#######################################################################
addresses=c("/VOT/Passivo/","/VOT/Ativo/","/Formantes/Passivo/","/Formantes/Ativo/")

beta <- readMat("beta.mat")
btVOT <- sort(beta$bt[1,subjects],index.return=TRUE)
btVOT$x <- btVOT$x/max(btVOT$x) #representa em escala entre 0 e 1
btForm <- sort(beta$bt[2,subjects],index.return=TRUE)
btForm$x <- btForm$x/max(btForm$x)
###phyVOT <- c(1, 34, 74, 125, 200)
###phyVOT <- c(-52, -40, -26, 8, 17)
st1 <- rep(-52,11)
st2 <- c(-29,-26,-51,-30,-50,-36,-20,-51,-36,-51,-23)
st3 <- c(-19,-20,-37,-18,-24,-20,-9,-37,-19,-37,-9)
st4 <- c(-9,-12,-12,8,8,8,10,-11,8,-14,10)
st5 <- rep(17,11)
ix=btVOT$ix
phyVOT <- c(st1[subjects[ix]],st2[subjects[ix]],st3[subjects[ix]],st4[subjects[ix]],st5[subjects[ix]])
###phyForm <- c(1, 55, 97, 146, 200)
ix=btForm$ix
###phyForm <- c(533, 961, 1103, 1297, 1387)#F2-F1
st1 <- rep(533,11)
st2 <- c(693,973,1061,1185,805,852,1036,865,988,533,1136)
st3 <- c(911,1197,1228,1147,1071,1067,1230,1067,1078,1186,1263)
st4 <- c(1228,1302,1222,1247,1176,1335,1209,1312,1246,1322,1315)
st5 <- rep(1387,11)
phyForm <- c(st1[subjects[ix]],st2[subjects[ix]],st3[subjects[ix]],st4[subjects[ix]],st5[subjects[ix]])
###phyForm <- c(1.703, 2.2779, 2.538, 3.0785, 3.5356)#F2/F1
for (stim in 1:5){

    for (ad in 1:4){

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
        matrix.psy.left <-eval(parse(text=sprintf("aux.left%d",stim)))
        matrix.psy.right <-eval(parse(text=sprintf("aux.right%d",stim)))


###############################################################################################
### Repete o beta de cada sujeito pela quantidade de estímulos dele(a)
        if (ad<=2){
            response1 <- rep (btVOT$x, rep (length(stimidx)/5,nb.subjects))
        }else{
            response1 <- rep (btForm$x, rep (length(stimidx)/5, nb.subjects))
        }
###############################################################################################        
        psy.EN.left <- list ()
        for (i in seq (1, length (alphas))) {
            psy.EN.left [[i]] <- cv.glmnet (matrix.psy.left [, 2 : ncol (matrix.psy.left)],
                                            response1,
                                            nfolds = 2,
                                            foldid = rep (seq (1, groups), 2*nb.subjects),
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
                                             nfolds = 2,
                                             foldid = rep (seq (1, groups), 2*nb.subjects),
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
            response2 <- rep (phyVOT, rep (groups*2,length(phyVOT)))
        }else{
            response2 <- rep (phyForm, rep (groups*2,length(phyForm)))
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
        par(mfrow=c(2,2))

        R1=round(cor(proj.psy.left,response1),3)
        test=t.test(proj.psy.left, response1, alternartive="two.sided", paired = TRUE, conf.level=0.95)
        plot(response1,proj.psy.left,main=sprintf("Psy-left R: %g t-test: %g", R1, test$p.value),xlab="Beta values", ylab="Predictions of beta values")
        abline(lm(response1~proj.psy.left))

###################################################################################################################
        R2=round(cor(proj.psy.right,response1),3)
        test=t.test(proj.psy.right, response1, alternartive="two.sided", paired = TRUE, conf.level=0.95)
        plot(response1,proj.psy.right,main=sprintf("Psy-right R: %g t-test: %g", R2, test$p.value),xlab="Beta values", ylab="Predictions of beta values")
        abline(lm(response1~proj.psy.right))

##########################################################################################################
        R3=round(cor(proj.phy.left,response2),3)
        test=t.test(proj.phy.left, response2, alternartive="two.sided", paired = TRUE, conf.level=0.95)
        plot(response2,proj.phy.left,main=sprintf("Phy-left R: %g t-test: %g", R3, test$p.value),xlab="Stimulus characteristic", ylab="Predictions of stimulus values")
        abline(lm(response2~proj.phy.left))

###################################################################################################################
        R4=round(cor(proj.phy.right,response2),3)
        test=t.test(proj.phy.right, response2, alternartive="two.sided", paired = TRUE, conf.level=0.95)
        plot(response2,proj.phy.right,main=sprintf("Phy-right R: %g t-test: %g", R4, test$p.value),xlab="Stimulus characteristic", ylab="Predictions of stimulus values")
        abline(lm(response2~proj.phy.right))

        dev.copy2pdf(file = file.path (sprintf("psy%d/Predictions-%s.pdf", stim, c("VOTpass","VOTact","FormPass","FormAct")[ad])))

###############################################################################################################################    
        coef.phy <- c(coef.phy.left[2:length(coef.phy.left)], coef.phy.right[2:length(coef.phy.right)])
        coef.psy <- c(coef.psy.left[2:length(coef.psy.left)], coef.psy.right[2:length(coef.psy.right)])

        save(AIC.coef,phy.EN.right.best, psy.EN.right.best, phy.EN.left.best, psy.EN.left.best, coef.phy.left, coef.phy.right, coef.psy.left, coef.psy.right, file=file.path(sprintf("psy%d/SIG_all_%s.dat",stim,c("VOTpass","VOTact","FormPass","FormAct")[ad])))
        
##########################################################################################################    
        filename=sprintf("psy%d/SIG-matrices-all-%s.mat",stim,c("VOTpass","VOTact","FormPass","FormAct")[ad])
        writeMat(filename,R1=R1, R2=R2, R3=R3, R4=R4,coef.phy=coef.phy, coef.psy=coef.psy, AIC.coef=AIC.coef)

    }#addresses
}#stim
