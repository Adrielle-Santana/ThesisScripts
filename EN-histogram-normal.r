### Generates the frequency (histogram) representation of the wavelets for all subjects
### in a time frequency wavelets representation
### plot histogram with vertical axis near normal sizes expected

histogram <- function(coefs.phy, coefs.psy, mw, sound){

    library(R.matlab)
    library(wavelets)
    library(glmnet)
    library(graphics)
    library (MASS)
    library (RColorBrewer)
    library (binhf)
    library (SDMTools)
    library(tuneR)
    
### Values to adjust width and height of the plot and the color bar
    hgt <-nb.dwt.levels-(2*(from.level-1))
    xneg <- -(16/(2^(from.level-1)))*4
    wdt <-c(60, 30, 15, 10, 7, 6, 4, 4)*1.4
    wdt2 <-c(8,7,6,5,3,2,1.5,1.2)
    xneg2 <- xneg*0.9
    divs <- c(9, 8, 7, 6, 5, 4, 3, 2)

### Generated by EN-matrix-coefs.r

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
            
    serie.temp <- seq(0,nb.points-1)/fs.eeg

### ROI
    start <- 0
    if (from.level > 1) start <- sum ((2 ^ seq ((log2(nb.points)-1), 1, by = -1)) [1 : (from.level - 1)])
    aux <- nb.points-start

### Places 1 where the lasso coefficients are different of zero
    psy.left <- cbind(matrix(rep(0, nb.subjects*start), nrow = nb.subjects), coefs.psy[,1:aux])
    phy.left <- cbind(matrix(rep(0, nb.subjects*start), nrow = nb.subjects), coefs.phy[,1:aux])
    psy.right <- cbind(matrix(rep(0, nb.subjects*start), nrow = nb.subjects), coefs.psy[,(aux+1):(aux*2)])
    phy.right <- cbind(matrix(rep(0, nb.subjects*start), nrow = nb.subjects), coefs.phy[,(aux+1):(aux*2)])

    psy.left[which (psy.left != 0)] <- 1
    phy.left[which (phy.left != 0)] <- 1
    psy.right[which (psy.right != 0)] <- 1
    phy.right[which (phy.right != 0)] <- 1

    t.psy.left <- table (colSums(psy.left))

    t.phy.left <- table (colSums(phy.left))

    t.psy.right <- table (colSums(psy.right))

    t.phy.right <- table (colSums(phy.right))

### Amount of colors necessary
    ma <-max(c(max(as.numeric(attributes(t.psy.left)$dimnames[[1]])), max(as.numeric(attributes(t.psy.right)$dimnames[[1]])),
               max(as.numeric(attributes(t.phy.left)$dimnames[[1]])),max(as.numeric(attributes(t.phy.right)$dimnames[[1]]))))

    colors <- rev(topo.colors(ma))
### Put color just from 2 repetitions of the coefficient. 0 and 1 repetition are represented by "white".
    colors[1] <- "white"
#######################################################################################################
### Creating new dwt objects with the frequency of each coefficient in all subjects
    obj.phy.left <- dwt(rep(0,nb.points), n.levels=nb.dwt.levels)
    obj.phy.right <- dwt(rep(0,nb.points), n.levels=nb.dwt.levels)

    initial <-1
    final <-nb.points/2
    leng <- nb.points/2

###Filling in the slots W and V for the physical response after unormalize data (not adding the average here)
    wp <- c(colSums(phy.left),colSums(phy.right))

    for (k in seq (1, nb.dwt.levels)) {
        if (k>=from.level){
            obj.phy.left@W[[k]] <- matrix(wp [initial : final], ncol=1)
            obj.phy.left@V[[k]] <- matrix(rep(0,leng), ncol=1)
            obj.phy.right@W[[k]] <- matrix(wp [(initial+(nb.points)) : (final+(nb.points))], ncol=1)
            obj.phy.right@V[[k]] <- matrix(rep(0,leng), ncol=1)
        }
        initial <- final + 1
        final <- final + (leng/2)
        leng <- leng / 2 
    }
    obj.phy.left@V[[nb.dwt.levels]] <- matrix(wp [initial : (nb.points)], ncol=1)
    obj.phy.right@V[[nb.dwt.levels]] <- matrix(wp [(initial+(nb.points)) : (nb.points*2)], ncol=1)
########################################################################################################
    obj.psy.left <- dwt(rep(0,nb.points), n.levels=nb.dwt.levels)
    obj.psy.right <- dwt(rep(0,nb.points), n.levels=nb.dwt.levels)

    initial <-1
    final <-nb.points/2
    leng <- nb.points/2

###Filling in the slots W and V for the behavioral response after unormalize data (not adding the average here)
    wb <- c(colSums(psy.left),colSums(psy.right))

    for (u in seq (1, nb.dwt.levels)) {
        if (u>=from.level){
            obj.psy.left@W [[u]] <- matrix(wb [initial : final], ncol=1)
            obj.psy.left@V[[u]] <- matrix(rep(0,leng), ncol=1)
            obj.psy.right@W [[u]] <- matrix(wb [(initial+(nb.points)) : (final+(nb.points))], ncol=1)
            obj.psy.right@V[[u]] <- matrix(rep(0,leng), ncol=1)
        }
        initial <- final + 1
        final <- final + (leng/2)
        leng <- leng / 2 
    }
    obj.psy.left@V[[nb.dwt.levels]] <- matrix(wb [initial : (nb.points)], ncol=1)
    obj.psy.right@V[[nb.dwt.levels]] <- matrix(wb [(initial+(nb.points)) : (nb.points*2)], ncol=1)
#############################################################################################################
### Shift in wavelet coeficients to be considered for the plot

    obj.phy.left <- align(obj.phy.left)
    obj.phy.right <- align(obj.phy.right)
    obj.psy.left <- align(obj.psy.left)
    obj.psy.right <- align(obj.psy.right)

    obj.list <- list(phy.left=obj.phy.left, phy.right=obj.phy.right, psy.left=obj.psy.left, psy.right=obj.psy.right)

    graphics.off()
    dev.new(width=11, height=9)
    layout (t (matrix (seq(1,4), ncol = 2)))

### Construction of the wavelets graphs Time X Frequency 

    for (c in 1:length(obj.list)){
        
### Maximum amount of squares necessay in the one line of the graph
        last.level.length <- length(obj.list[[c]]@W[[from.level]]) 

### Creating the plot window

        plot (c (xneg, (last.level.length+wdt[from.level])), c (0, (from.level+hgt)), type = "n", xaxt = "n", yaxt = "n", bty = "n", xlab = "time (s)", ylab = "frequency (Hz)", main=names(obj.list[c]))
        
        pnts = cbind(x =c((last.level.length+wdt2[from.level]),(last.level.length+1),(last.level.length+1),(last.level.length+wdt2[from.level])), y =c((from.level+hgt),(from.level+hgt),0,0))
        
        legend.gradient(pnts,c("white",colors),limits=" ", title="")

        count <- 1
        df <- 0.2

        for (h in nb.dwt.levels:from.level) {

###Apply the coresponding shift to each W level
            
            amount.squares <- length(obj.list[[c]]@W[[h]])
            width.square <- last.level.length/amount.squares
            
            x <- c (0, 0, width.square, width.square)

            if (h == nb.dwt.levels){
                
                y <- c (0, 1, 1, 0)*df

###Apply the coresponding shift to the last V level only
            }
            else{
                y <- c (0, 1, 1, 0) * (count*df)
            }
            
            for (w in 1:amount.squares){
                
                if (h == nb.dwt.levels){

                    if(obj.list[[c]]@V[[h]][w] != 0){    

                        idx <- obj.list[[c]]@V[[h]][w]

                        polygon (x, y, col=colors[idx], border="grey")
                    }
                    else{
                        polygon (x, y, col="white", border="grey")
                    }
                }
                if (obj.list[[c]]@W[[h]][w] != 0){
                    
                    idx <- obj.list[[c]]@W[[h]][w]

                    polygon (x, y+(count*df), col=colors[idx], border="grey")
                }
                else{
                    polygon (x, y+(count*df), col="white", border="grey")
                }

                x <- x + c (width.square, width.square, width.square, width.square)
            }
            if (h == nb.dwt.levels){
                
                text(x=xneg,y=(y[2])/2, labels=sprintf("V%d",h))
            }
            
            text(x=xneg,y=(y[2]+0.5), labels=sprintf("W%d",h))
            
            count <- 2*count
        }

                                       
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

        axis(4,at=seq(0.2, (from.level+hgt-0.2), l=ma+1),pos=c(last.level.length+2,0),labels=seq(0,ma),lwd = 0,las=1)
        
    }
}#function

