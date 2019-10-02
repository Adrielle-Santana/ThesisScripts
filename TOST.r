### Esse scritp testa a igualdade dos coeficientes EN para os dados psicofísicos obtidos individualmente para cada estímulo. Espera-se que sejam semelhantes para os stim1=stim2 e stim4=stim5. Stim3 é uma incógnita aqui e talvez pode ser explicado pela percepção física desse estímulo. Melhor ainda se não for próximo de nenhum das extremidades.
rm(list=ls())

library(effsize)
library(equivalence)
library(ANOVAreplication)

nomes=c("SIG_all_VOTact.dat", "SIG_all_VOTPass.dat", "SIG_all_FormAct.dat", "SIG_all_Formpass.dat")

for (arq in 1:4){
    
    load(sprintf("psy1/%s",nomes[arq]))
    left=as.matrix(coef.psy.left)
    left=left[2:length(left)]
    right=as.matrix(coef.psy.right)
    right=right[2:length(right)]

    load(sprintf("psy2/%s",nomes[arq]))
    left2=as.matrix(coef.psy.left)
    left2=left2[2:length(left2)]
    right2=as.matrix(coef.psy.right)
    right2=right2[2:length(right2)]

    load(sprintf("psy3/%s",nomes[arq]))
    left3=as.matrix(coef.psy.left)
    left3=left3[2:length(left3)]
    right3=as.matrix(coef.psy.right)
    right3=right3[2:length(right3)]

    load(sprintf("psy4/%s",nomes[arq]))
    left4=as.matrix(coef.psy.left)
    left4=left4[2:length(left4)]
    right4=as.matrix(coef.psy.right)
    right4=right4[2:length(right4)]

    load(sprintf("psy5/%s",nomes[arq]))
    left5=as.matrix(coef.psy.left)
    left5=left5[2:length(left5)]
    right5=as.matrix(coef.psy.right)
    right5=right5[2:length(right5)]

    t=seq(0:(length(left)-1))/5000

###stim1 com stim2
    sp = sqrt(((length(left) - 1) * var(left) +(length(left2) - 1) * var(left2) / (length(left)+length(left2) - 2)))
    d=(mean(left)-mean(left2))/sp
    ###d=cohen.d(left,left2,pooled=TRUE,paired=TRUE)$estimate
    tost(left,left2,epsilon=abs(d),paired=TRUE, var.equal=TRUE, alpha=0.025)

    sp =  sqrt(((length(right) - 1) * var(right) +(length(right2) - 1) * var(right2) / (length(right)+length(right2) - 2)))
    d=(mean(right)-mean(right2))/sp
    ###d=abs(cohen.d(right,right2,pooled=TRUE,paired=TRUE)$estimate)
    tost(right,right2,epsilon=abs(d),paired=TRUE, var.equal=TRUE)

###stim4 com stim5
    sp = sqrt(((length(left4) - 1) * var(left4) +(length(left5) - 1) * var(left5) / (length(left4)+length(left5) - 2)))
    d=(mean(left4)-mean(left5))/sp
    d=cohen.d(left4,left5,pooled=TRUE,paired=TRUE)$estimate
    tost(left4,left5,epsilon=abs(d),paired=TRUE, var.equal=TRUE,alpha=0.025)

    sp =  sqrt(((length(right4) - 1) * var(right4) +(length(right5) - 1) * var(right5) / (length(right4)+length(right5) - 2)))
    d=(mean(right4)-mean(right5))/sp
    d=abs(cohen.d(right4,right5,pooled=TRUE,paired=TRUE)$estimate)
    tost(right4,right5,epsilon=d,paired=TRUE, var.equal=TRUE)

### stim4 com stim2
    sp = sqrt(((length(left4) - 1) * var(left4) +(length(left2) - 1) * var(left2) / (length(left4)+length(left2) - 2)))
    d=(mean(left4)-mean(left2))/sp
    d=cohen.d(left2,left4,pooled=TRUE,paired=TRUE)$estimate
    tost(left2,left4,epsilon=abs(d),paired=TRUE, var.equal=TRUE,alpha=0.025)

    sp =  sqrt(((length(right4) - 1) * var(right4) +(length(right2) - 1) * var(right2) / (length(right4)+length(right2) - 2)))
    d=(mean(right4)-mean(right2))/sp
    d=abs(cohen.d(right2,right4,pooled=TRUE,paired=TRUE)$estimate)
    tost(right2,right4,epsilon=d,paired=TRUE, var.equal=TRUE)

### stim1 com stim4
    sp = sqrt(((length(left) - 1) * var(left) +(length(left4) - 1) * var(left4) / (length(left)+length(left4) - 2)))
    d=(mean(left)-mean(left4))/sp
    ###d=cohen.d(left,left4,pooled=TRUE,paired=TRUE)$estimate
    tost(left,left4,epsilon=abs(d),paired=TRUE, var.equal=TRUE,,alpha=0.025)

    sp = sqrt(((length(right4) - 1) * var(right4) +(length(right) - 1) * var(right) / (length(right4)+length(right) - 2)))
    d=(mean(right)-mean(right4))/sp
    #d=abs(cohen.d(right,right4,pooled=TRUE,paired=TRUE)$estimate)
    tost(right,right4,epsilon=d,paired=TRUE, var.equal=TRUE)

### stim1 com stim3
    sp = sqrt((((length(left) - 1) * var(left)) +((length(left3) - 1) * var(left3))) / (length(left)+length(left3) - 2))
    d=(mean(left)-mean(left3))/sp
    d=cohen.d(left,left3,pooled=TRUE,paired=TRUE)$estimate
    tost(left,left3,epsilon=abs(d),paired=TRUE, var.equal=TRUE, alpha=0.025)

    sp = sqrt(((length(right3) - 1) * var(right3) +(length(right) - 1) * var(right) / (length(right3)+length(right) - 2)))
    d=(mean(right)-mean(right3))/sp
    #d=abs(cohen.d(right,right4,pooled=TRUE,paired=TRUE)$estimate)
    tost(right,right3,epsilon=abs(d),paired=TRUE, var.equal=TRUE)

###stim1 com stim5
    sp = sqrt(((length(left) - 1) * var(left) +(length(left5) - 1) * var(left5) / (length(left)+length(left5) - 2)))
    d=(mean(left)-mean(left5))/sp
    d=cohen.d(left,left5,pooled=TRUE,paired=TRUE)$estimate
    tost(left,left5,epsilon=abs(d),paired=TRUE, var.equal=TRUE,0.025)

    sp = sqrt(((length(right5) - 1) * var(right5) +(length(right) - 1) * var(right) / (length(right5)+length(right) - 2)))
    d=(mean(right)-mean(right5))/sp
    d=abs(cohen.d(right,right5,pooled=TRUE,paired=TRUE)$estimate)
    tost(right,right5,epsilon=d,paired=TRUE, var.equal=TRUE)

###stim3 com stim5
    sp = sqrt(((length(left3) - 1) * var(left3) +(length(left5) - 1) * var(left5)) / (length(left3)+length(left5) - 2))
    d=(mean(left3)-mean(left5))/sp
    #d=cov(left3,left5)/(sd(left3)*sd(left5)) #Pearson
    #d=cohen.d(left3,left5,pooled=TRUE,paired=TRUE)$estimate
    tost(left3,left5,epsilon=abs(d),paired=TRUE, var.equal=TRUE,0.025)

    t.test(left3,left5,alternative="two.sided",paired=TRUE,var.equal=TRUE,mu=abs(d))

    sp = sqrt(((length(right3) - 1) * var(right3) +(length(right5) - 1) * var(right5) / (length(right3)+length(right5) - 2)))
    d=(mean(right5)-mean(right3))/sp
    d=cohen.d(right3,right5,pooled=TRUE,paired=TRUE)$estimate
    tost(right3,right5,epsilon=d,paired=TRUE, var.equal=TRUE)

}
