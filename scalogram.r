### Function to plot scalogram
### Inputs:
###         coef.phy and coef.psy: DWT coefficiet vectors without intercept value

###         from.level: Provide the DWT level from which coefficients will be represented.

###         mw: string with the mother wavelet used. Ex.: "la8"

### Adrielle C. Santana

scalogram <-  function(coef.phy, coef.psy, from.level, mw, sound){

    library(R.matlab)
    library(wavelets)
    library(glmnet)
    library(graphics)
    library (MASS)
    library (RColorBrewer)
    library (binhf)
    library (SDMTools)
    library(tuneR)
    
    start <- 0
    points <- nb.points
    
    if (from.level>1){
        start <- sum ((2 ^ seq ((log2(nb.points)-1), 1, by = -1)) [1 : (from.level - 1)])
        points <- nb.points-start
    }
    
###Creating the object to compute the idwt for the physical response
    obj.phy.left <- dwt(rep(0,nb.points), n.levels=nb.dwt.levels, filter=mw)
    obj.phy.right <- dwt(rep(0,nb.points), n.levels=nb.dwt.levels, filter=mw)

    initial <-1
    final <-nb.points/2
    leng <- nb.points/2

###Filling in the slots W and V for the physical response

    coef.p <- coef.phy

    wp <- c(rep(0,start),coef.p[1:points],rep(0,start),coef.p[(points+1):(points*2)])

    for (k in seq (1, nb.dwt.levels)) {
        if (k>=from.level){
            obj.phy.left@W[[k]] <- matrix(wp [initial : final], ncol=1)
            obj.phy.left@V[[k]] <- matrix(rep(0,leng), ncol=1)
            obj.phy.right@W[[k]] <- matrix(wp [(initial+nb.points) : (final+nb.points)], ncol=1)
            obj.phy.right@V[[k]] <- matrix(rep(0,leng), ncol=1)
        }
        initial <- final + 1
        final <- final + (leng/2)
        leng <- leng / 2 
    }
    obj.phy.left@V[[nb.dwt.levels]] <- matrix(wp [initial : nb.points], ncol=1)
    obj.phy.right@V[[nb.dwt.levels]] <- matrix(wp [(initial+nb.points) : (nb.points*2)], ncol=1)

###Creating the object to compute the idwt for the behavioral response

    obj.psy.left <- dwt(rep(0,nb.points), n.levels=nb.dwt.levels, filter=mw)
    obj.psy.right <- dwt(rep(0,nb.points), n.levels=nb.dwt.levels, filter=mw)

    initial <-1
    final <-nb.points/2
    leng <- nb.points/2

###Filling in the slots W and V for the behavioral response

    coef.b <- coef.psy

    wb <- c(rep(0,start),coef.b[1:points],rep(0,start),coef.b[(points+1):(points*2)])

    for (u in seq (1, nb.dwt.levels)) {
        if (u>=from.level){
            obj.psy.left@W [[u]] <- matrix(wb [initial : final], ncol=1)
            obj.psy.left@V[[u]] <- matrix(rep(0,leng), ncol=1)
            obj.psy.right@W [[u]] <- matrix(wb [(initial+nb.points) : (final+nb.points)], ncol=1)
            obj.psy.right@V[[u]] <- matrix(rep(0,leng), ncol=1)
        }
        initial <- final + 1
        final <- final + (leng/2)
        leng <- leng / 2 
    }
    obj.psy.left@V[[nb.dwt.levels]] <- matrix(wb [initial : nb.points], ncol=1)
    obj.psy.right@V[[nb.dwt.levels]] <- matrix(wb [(initial+nb.points) : (nb.points*2)], ncol=1)

###Calculates the inverse wavelets and fill the field "series" in the dwt object 
    inv.obj.phy.left <- idwt(obj.phy.left)
    obj.phy.left@series <- matrix(inv.obj.phy.left, ncol=1)
    inv.obj.phy.right <- idwt(obj.phy.right)
    obj.phy.right@series <- matrix(inv.obj.phy.right, ncol=1)


    inv.obj.psy.left <- idwt(obj.psy.left)
    obj.psy.left@series <- matrix(inv.obj.psy.left, ncol=1)
    inv.obj.psy.right <- idwt(obj.psy.right)
    obj.psy.right@series <- matrix(inv.obj.psy.right, ncol=1)

##########################################################################################
### Values to adjust width and height of the plot and the color bar
hgt <-nb.dwt.levels-(2*(from.level-1))
xneg <- -(16/(2^(from.level-1)))*4
wdt <-c(60, 30, 15, 10, 7, 6, 4, 4)*1.4
wdt2 <-c(8,7,6,5,3,2,1.5,1.2)
xneg2 <- xneg*0.9
divs <- c(9, 8, 7, 6, 5, 4, 3, 2)

    if (sound==1){
        st1 <- readWave("da.wav")
        nst1 <- "da"
        st2 <- readWave("ta.wav")
        nst2 <- "ta"
        div <- 30000
    }
    else{
        st1 <- readWave("pa.wav")
        nst1 <- "pa"
        st2 <- readWave("pe.wav")
        nst2 <- "pe"
        div <- 10000
    }
    
start <- sum ((2 ^ seq (10, 1, by = -1)) [1 : (from.level - 1)])
points <- nb.points-start

### Shift in wavelet coeficients to be considered for the plot
### It will be the same fot the other objects
obj.phy.left <- align(obj.phy.left)
obj.phy.right <- align(obj.phy.right)
obj.psy.left <- align(obj.psy.left)
obj.psy.right <- align(obj.psy.right)

### It's important to keep this sequency as it is below so the colors will be right
obj.list <- list(phy.left=obj.phy.left, phy.right=obj.phy.right, psy.left=obj.psy.left, psy.right=obj.psy.right)

graphics.off()
dev.new(width=11, height=9)
layout (t (matrix (seq(1,4), ncol = 2)))

### Construction of the wavelets graphs Time X Frequency 

for (c in 1:length(obj.list)){

### Ordering the coefficients in crescent order to tune color scale
    ma1 <- ma2 <- ma3 <- ma4<- 0
    for (i in from.level:nb.dwt.levels ){
        aux <- max(abs(obj.phy.left@W[[i]]))
        if (aux>ma1) ma1 <- aux

        aux <- max(abs(obj.phy.right@W[[i]]))
        if (aux>ma2) ma2 <- aux

        aux <- max(abs(obj.psy.left@W[[i]]))
        if (aux>ma3) ma3 <- aux

        aux <- max(abs(obj.psy.right@W[[i]]))
        if (aux>ma4) ma4 <- aux
    }
    
###ma.list <- list(ma1, ma2, ma3, ma4)
    ma.list <- rep(list(max(ma1, ma2, ma3, ma4)), 4)
    
### Color vector
    
    colors <- colorRampPalette (c ("blue", "white", "red")) (50)

    serie.temp <- seq(0,(length(obj.list[[c]]@series)-1))/fs.eeg 

### Maximum amount of squares necessay in the one line of the graph
    last.level.length <- length(obj.list[[c]]@W[[from.level]]) 

### Creating the plot window
###            dev.new(width=11, height=9)

    plot (c (xneg, (last.level.length+wdt[from.level])), c (0, (from.level+hgt)), type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "time (s)", ylab = "frequency (Hz)", main=names(obj.list[c]))
    
    pnts = cbind(x =c((last.level.length+wdt2[from.level]),(last.level.length+1),(last.level.length+1),(last.level.length+wdt2[from.level])), y =c((from.level+hgt),(from.level+hgt),0,0))

    legend.gradient(pnts,colors,limits=c(round(-(ma.list[[c]]),2),round(ma.list[[c]],2)), title="")

    count <- 0

    for (h in nb.dwt.levels:from.level) {

###Apply the coresponding shift to each W level
        
        amount.squares <- length(obj.list[[c]]@W[[h]])
        width.square <- last.level.length/amount.squares
        
        x <- c (0, 0, width.square, width.square)

        if (h == nb.dwt.levels){
            
            y <- c (0, 1, 1, 0)
        }
###Apply the coresponding shift to the last V level only

        else{
            y <- c (0, 1, 1, 0) + count
        }
        
        for (w in 1:amount.squares){
            
            if (h == nb.dwt.levels){

                if(obj.list[[c]]@V[[h]][w] != 0){

                    idx <- round(((obj.list[[c]]@V[[h]][w]+(ma.list[[c]]))/(2*ma.list[[c]]))*50)
                    
                    
                    if (idx<1){
                        polygon (x, y, col=colors[1], border="grey")
                    }
                    else if (idx>50){
                        polygon (x, y, col=colors[50], border="grey")
                    }
                    else{
                        polygon (x, y, col=colors[idx], border="grey")
                    }
                    
                }
                else{
                    polygon (x, y, col="white", border="grey")
                }
            }
            if (obj.list[[c]]@W[[h]][w] != 0){

                idx <- round(((obj.list[[c]]@W[[h]][w]+(ma.list[[c]]))/(2*ma.list[[c]]))*50)

                if (idx<1){
                    polygon (x, y+1, col=colors[1], border="grey")
                }
                else if (idx>50){
                    polygon (x, y+1, col=colors[50], border="grey")
                }
                else{
                    polygon (x, y+1, col=colors[idx], border="grey")
                }
            }
            else{
                polygon (x, y+1, col="white", border="grey")
            }

            x <- x + c (width.square, width.square, width.square, width.square)
        }
        if (h == nb.dwt.levels){
            
            text(x=xneg,y=(y[2])/2, labels=sprintf("V%d",h))
        }
        
        text(x=xneg,y=(y[2]+0.5), labels=sprintf("W%d",h))
        
        count <- count+1
    }
    
    lines((serie.temp*last.level.length)/(max(serie.temp)), (obj.list[[c]]@series/(ma.list[[c]]*0.2))+3, type='l')

    aux2 <- fs.eeg
    
    for (f in (nb.dwt.levels+1):0){

        aux2 <- round(aux2/2,2)

        if (f <= divs[from.level]){
            
            if (f == 0){
                text(x=xneg2,y=f, labels= sprintf("%g", 0), col="red")
            }
            else{
                text(x=xneg2,y=f, labels= sprintf("%g", aux2), col="red")
            }
        }
    }
    

    temp <- seq(1:length(st1@left))/st1@samp.rate

    leng <- (last.level.length*(length(temp)/st1@samp.rate))/(length(serie.temp)/fs.eeg)

    lines((temp*leng)/max(temp), (st1@left/div)+1.5, type='l', col="#00000090")
    text(x=leng-1,y=1.6, labels=nst1)
    lines((temp*leng)/max(temp), (st2@left/div)+0.5, type='l', col="#00000090")
    text(x=leng-1,y=0.6, labels=nst2)

    nb.div <- length(obj.list[[1]]@W[[nb.dwt.levels]])+1
    lab <- round(seq(0,serie.temp[length(serie.temp)],l=nb.div),2)
    
    for (s in seq(1:nb.div)){
        lab[s] <- toString(lab[s])
    }
    
    axis(1,at=seq(0, last.level.length, l=nb.div),labels=lab,lwd = 0.2)
    
    axis(4,at=3,pos=c(last.level.length+0.5,3),labels="Coef.Values", tick = FALSE)
}
    
}#function
