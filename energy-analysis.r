### Analyses of regression directions obtained from RoLDSIS

### * Load local libraries
source ("paths.r")
source ("dwt-lib.r")
source ("experiments.r")
source ("sink.r")
source ("contrasts.r")

### * Load system libraries
load.pkgs (c ("lme4", "lmerTest", "emmeans"))
load.pkgs ("wavelets")

### Directory to save figures
roldsis.dir <- file.path (figures.dir, "Roldsis-analysis")
force.dir.create (roldsis.dir)

### * Initialize data frames
unique.df <- data.frame ()
roi.df <- data.frame ()

### * Create factors: electrodes and wavelet bands
electrodes <- c ("F7", "F8", "Fz", "TP9", "TP10")
wav.band <- c ("W9","W8","W7","W6W5")

### * Load previously computed direction vectors (by RoLDSIS)
load (file = file.path (results.dir, sprintf("dirVec-%d.dat",eeg.sampfreq)))

### * Fill data frame
for (ef in experiment.features){

    for (et in experiment.types){

        for (subj in cohort){

            S <- sprintf("subj%d",subj)

            for (reg in c("psy","phy")){

                for (elt in electrodes) {

                    ## ** Select direction vector and convert to wavelet object
                    aux <- dirVec [[ef]] [[et]] [[S]] [[reg]] [[elt]]
                    obj <- vec.to.dwt(aux, n=dwt.length)
                    
                    ## ** Define ROIs

                    #e1 <-  (mean ((obj@W [[9]] [1:2])))
                    #e2 <-  (mean ((obj@W [[9]] [3])))
                    #e3 <-   (mean ((obj@W [[8]] [1:3])))
                    #e4 <-  (mean ((obj@W [[8]] [4:6])))
                    #e5 <-   (mean ((obj@W [[7]] [1:6])))
                    #e6 <-  (mean ((obj@W [[7]] [7:12])))
                    #e7 <- mean (c(mean ((obj@W [[6]] [1:12])),  mean ((obj@W [[5]] [1:24]))))
                    #e8 <- mean (c(mean ((obj@W [[6]] [13:24])), mean ((obj@W [[5]] [25:46]))))

                    e1 <-  (mean ((obj@W [[9]] [1])))
                    e2 <-  (mean ((obj@W [[9]] [2])))
                    e3 <-  (mean ((obj@W [[8]] [1:2])))
                    e4 <-  (mean ((obj@W [[8]] [3:4])))
                    e5 <-  (mean ((obj@W [[7]] [1:4])))
                    e6 <-  (mean ((obj@W [[7]] [5:8])))
                    e7 <-  mean(c(mean ((obj@W [[6]] [1:8])), mean((obj@W [[5]] [1:16]))))
                    e8 <-  mean(c(mean ((obj@W [[6]] [9:16])), mean ((obj@W [[5]] [17:32]))))


                    ## ** Store data in the data frame 1
                    roi.df <- rbind (roi.df, data.frame(subject = subj,
                                                        electrode = elt,
                                                        regression = reg,
                                                        feature = ef,
                                                        type = et,
                                                        roi1 = e1,
                                                        roi2 = e2,
                                                        roi3 = e3,
                                                        roi4 = e4,
                                                        roi5 = e5,
                                                        roi6 = e6,
                                                        roi7 = e7,
                                                        roi8 = e8))

                    aux2 <- 1
                    ## ** Create auxiliary variables to store size and coefficients
                    ## from each band of the wavelet object (for data frame unique.df)
                    for(bd in seq(1,length(wav.band))){

                        for (tm in c ("early", "late")) {

                            
                            e <- eval(parse(text = sprintf("e%d", aux2)))
                            
                            aux2 <- aux2+1

                            ## *** Store data in the data frame 2
                            unique.df <- rbind (unique.df,
                                                data.frame(subject = subj,
                                                           electrode = elt,
                                                           regression = reg,
                                                           feature = ef,
                                                           band = wav.band [bd],
                                                           type = et,
                                                           time = tm,
                                                           var = e))

                        } ## tm
                    } ## bd
                } ## elt
            } ## reg
        } ## subj
    } ## et
} ## ef

### One-way contrasts
contr <- list ("feature" = build.contr (roi.df$feature),
               "type" = build.contr (roi.df$type),
               "regression" = build.contr(roi.df$regression),
               "electrode" = list ("medial - lateral" = c (-1/4, -1/4, 1, -1/4, -1/4),
                                   "frontal - temporal" = c (1/2, 1/2, 0, -1/2, -1/2),
                                   "left - right" = c (1/2, -1/2, 0, -1/2, 1/2)))

### Two-way contrasts

n <- names (contr)
for (i in seq (1, length (n) - 1)) {
    ci <- contr [[n [i]]]
    for (j in seq (i + 1, length (n))) {
        cj <- contr [[n [j]]]
        m <- sprintf ("%s:%s", n [[i]], n [[j]])
        l <- list ()
        for (p in names (ci))
            for (q in names (cj))
                l [[sprintf ("%s : %s", p, q)]] <- prod.contr (ci [[p]], cj [[q]])
        contr [[m]] <- l
    }
}

### Check contrasts
banner ("Check contrasts", "*")
for (n in names (contr)) {
    banner (n, "=")
    check.contrast (contr [[n]])
}

### * Compute mixed effects models for both data frames

results.filestem <- "energy-analysis"
sink.open (file.path (results.dir, sprintf ("%s.txt", results.filestem)))

## Complete data frame unique.df
#banner ("Full model", "*")

#fm <- lmer (var ~ feature * type * electrode * regression * time * band
#            + (1 | subject), unique.df)
#step.val <- step (fm)
#fm.full <- get_model (step.val)
#show (anova (fm.full))
#show (ranova (fm.full))

## Data frame separating ROIs
for (var in c ("roi1", "roi2", "roi3", "roi4", "roi5", "roi6", "roi7", "roi8")) {

    banner (var, "*")

    frm <- as.formula (sprintf ("%s ~ feature * type * electrode * regression + (1 | subject)",
                                var))

    fm <- lmer (frm, roi.df)

     ## *** Global ANOVA table
    banner ("ANOVA table for full model", "=")
    show (anova (fm))

    ## *** Plot the original data
    pdf (file.path (roldsis.dir, sprintf ("%s-0full.pdf", var)))
    print (emmip (fm, regression ~ electrode | feature * type,
                  CIs = TRUE, ylab = var))
    dummy <- dev.off ()

    ## Stepwise selection
    step.val <- step (fm)
    fm <- get_model (step.val)

    ## Random effect
    banner ("Random effect ANOVA table for step-wise selected model", "=")
    show (ranova (fm))

    ## ANOVA table
    banner ("ANOVA table for step-wise selected model", "=")
    aov <- anova (fm)
    show (aov)

    ## Check effects

    ##  Get names of the variables in the model
    n <- names (attributes (aov) $ hypotheses)

    ##  Loop over effects
    for (i in seq (1, length (n))){

        ##  Only consider factors with significant effects
        if (aov [["Pr(>F)"]] [i] < 0.05) {

            ##  Get factor name
            eff <- n [[i]]

            ##  Only consider factors for which effects have been defined
            if (eff %in% names (contr)) {

                ##  Banner for effect
                banner (sprintf ("Effect: %s", eff), "=")

                ##  Get emmeans formula
                frm.em <- as.formula (sprintf ("~ %s", gsub (":", " * ", eff)))

                ##  Test the effect
                em <- emmeans (fm, frm.em)
                show (summary (em))
                show (contrast (em, method = contr [[eff]]))

                ## Plot the effect
                if (grepl (":", eff)){
                    frm.emmip <- as.formula (gsub (":", " ~ ", eff))
                } else{
                    frm.emmip <- as.formula (sprintf ("~ %s", eff))
                }

                pdf (file.path (roldsis.dir, sprintf ("%s-%s.pdf", var, eff)))
                print (emmip (fm, frm.emmip, ylab = var, CIs = TRUE))
                dummy <- dev.off ()

            } ## if eff
        } ## if aov
    } ## for i
} ## for var

sink.close ()

##  Compose the PDF file with the results for all ROIs
system (sprintf ("pdftk %s/*.pdf cat output %s/roldsis-analysis.pdf",
                 roldsis.dir, figures.dir))
