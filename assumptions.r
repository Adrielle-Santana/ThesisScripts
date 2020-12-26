### This script was developed to test the assumptions for the ANOVA computation
### Those assumptions include: homocedasticity and  normality of residues from fixed
### and random effects. Independency of data is assumed to be true because
### the data was acquired from different subjects in a randomized experiment.

### ** Load local library
source ("paths.r")
source ("experiments.r")
source ("dwt-lib.r")
source ("sink.r")

### ** Load system libraries
load.pkgs (c ("lme4", "lmerTest", "emmeans", "lattice", "sjPlot", "effects", "glmmTMB"))

### ** Load data
load (file.path (results.dir, "N1-P2.dat"))

### ** Transform stimulus into discrete factor
n1.p2.df$stimulus <- factor (n1.p2.df$stimulus)

### ** Eliminate bad cases
bad.cases <- read.csv (file.path (data.dir, "bad-cases.csv"))
bad.cases$electrode <- as.character (bad.cases$electrode)
for (i in seq (1, nrow (bad.cases))) {
    idx <- rep (TRUE, nrow (n1.p2.df))
    for (n in names (bad.cases))
        idx <- idx & (n1.p2.df [[n]] == bad.cases [[n]] [i])
    n1.p2.df <- n1.p2.df [which (! idx), ]
}

levels(n1.p2.df$feature)[levels(n1.p2.df$feature)=="Formantes"] <- "Formants"
levels(n1.p2.df$type)[levels(n1.p2.df$type)=="Passivo"] <- "Passive"
levels(n1.p2.df$type)[levels(n1.p2.df$type)=="Ativo"] <- "Active"

### * Redirect output to text file
results.filestem <- "Assumptions"
sink.open (file.path (results.dir, sprintf ("%s.txt", results.filestem)))

### * Directory for saving plots
assumptions.fig.dir <- file.path (figures.dir, results.filestem)
unlink (assumptions.fig.dir, recursive = TRUE)
force.dir.create (assumptions.fig.dir)

dep.vars <- names (n1.p2.df) [seq (6, 11)]

### ** Loop over dependent variables for N1-P2
for (var in dep.vars) {

    banner (sprintf("Variable: %s", var),"=")
    
    ## *** Model formula for the current dependent variable
    frm <- as.formula (sprintf ("%s  ~ feature * type * electrode * stimulus + (1 | subject)", var))
    ## *** Fit general model
    fm <- lmer (frm, n1.p2.df)

    if (var == "P2.N1") {

        ## **** Compensate for different ISI between types of experiment
        ## (Ativo × Passivo)

        ## **** Compute the multiplicative factor between Active and Passive
        ## Do it for each feature/subject/electrode/stimulus.
        ## The data frame below contains the ratio Active/Passive for each case.
        frm.scale <- as.formula (
            sprintf ("%s  ~ subject * feature * stimulus * electrode", var))
        df <- aggregate (frm.scale, n1.p2.df, function (x) x [2] / x [1])

        ## **** Get the median value
        act.pass.factor <- median (df [[var]], na.rm = TRUE)

        ## **** Correction factor for the Ativo cases
        idx <- which (n1.p2.df$type == "Active")
        n1.p2.df [[var]] [idx] <- n1.p2.df [[var]] [idx] / act.pass.factor

        ## **** Fit the model with “ISI correction”
        ## There should be no effect now for factor type
        fm <- lmer (frm, n1.p2.df)
    }

    banner("Normality test for residues of full model","=")
    show(shapiro.test(resid(fm)))
    pdf (file.path (assumptions.fig.dir, sprintf ("%s-QQfull.pdf", var)))
    qqmath(fm)
    dummy <- dev.off ()

    pdf (file.path (assumptions.fig.dir, sprintf ("%s-normal.pdf", var)))
    plot_model(fm, type='diag')[[3]]
    dummy <- dev.off ()    

    banner("Normality test for random effect","=")
    r=ranef(fm)$subject$`(Intercept)`
    show(shapiro.test(r))
    pdf (file.path (assumptions.fig.dir, sprintf ("%s-QQrandom.pdf", var)))
    qqnorm(r)
    qqline(r)
    dummy <- dev.off ()

    pdf (file.path (assumptions.fig.dir, sprintf ("%s-variance.pdf", var)))
    plot_model(fm, type='diag')[[4]]
    dummy <- dev.off () 
    
}
####################################################################################
### Analysis of assumptions for ROIs
### * Initialize data frame
roi.df <- data.frame ()

### * Create factors: electrodes and wavelet bands
electrodes <- c ("F7", "F8", "Fz", "TP9", "TP10")
wav.band <- c ("W9","W8","W7","W6W5")

### * Load previously computed direction vectors (by RoLDSIS)
load (file = file.path (results.dir, sprintf("dirVec-%d.dat",eeg.sampfreq)))

importance <- function (x, y)
    sum (x ^ 2) - sum (y ^ 2)

### * Fill data frame
for (ef in experiment.features){
    if (ef == "Formantes"){
        efAux = "Formants"
    }else{
        efAux = "VOT"
    }
    
    for (et in experiment.types){

        if (et == "Passivo"){
            etAux = "Passive"
        }else{
            etAux = "Active"
        }
        
        for (subj in cohort){

            S <- sprintf("subj%d",subj)

            for (reg in c("psy","phy")){

                for (elt in electrodes) {

                    ## ** Select direction vector and convert to wavelet object
                    aux.phy <- dirVec [[ef]] [[et]] [[S]] $phy [[elt]]
                    obj.phy <- vec.to.dwt (aux.phy, n=dwt.length)

                    aux.psy <- dirVec [[ef]] [[et]] [[S]] $psy [[elt]]
                    obj.psy <- vec.to.dwt (aux.psy, n=dwt.length)
                    
                    e1 <-  importance (obj.phy@W [[9]] [1], obj.psy@W [[9]] [1])
                    e2 <-  importance (obj.phy@W [[9]] [2], obj.psy@W [[9]] [2])
                    e3 <-  importance (obj.phy@W [[8]] [1:2], obj.psy@W [[8]] [1:2])
                    e4 <-  importance (obj.phy@W [[8]] [3:4], obj.psy@W [[8]] [3:4])
                    e5 <-  importance (obj.phy@W [[7]] [1:4], obj.psy@W [[7]] [1:4])
                    e6 <-  importance (obj.phy@W [[7]] [5:8], obj.psy@W [[7]] [5:8])
                    e7 <-  (importance (obj.phy@W [[6]] [1:8], obj.psy@W [[6]] [1:8])
                            + importance (obj.phy@W [[5]] [1:16], obj.psy@W [[5]] [1:16]))
                    e8 <-  (importance (obj.phy@W [[6]] [9:16], obj.psy@W [[6]] [9:16])
                            + importance (obj.phy@W [[5]] [17:32], obj.psy@W [[5]] [17:32]))


                    ## ** Store data in the data frame 1
                    roi.df <- rbind (roi.df, data.frame(subject = subj,
                                                        electrode = elt,
                                                        regression = reg,
                                                        feature = efAux,
                                                        type = etAux,
                                                        roi1 = e1,
                                                        roi2 = e2,
                                                        roi3 = e3,
                                                        roi4 = e4,
                                                        roi5 = e5,
                                                        roi6 = e6,
                                                        roi7 = e7,
                                                        roi8 = e8))
                } ## elt
            } ## reg
        } ## subj
    } ## et
} ## ef

### Loop over dependent variable ROI

## Data frame separating ROIs
for (var in c ("roi1", "roi2", "roi3", "roi4", "roi5", "roi6", "roi7", "roi8")) {

    banner (sprintf("Variable: %s", var),"=")

    frm <- as.formula (sprintf ("%s ~ feature * type * electrode + (1 | subject)", var))

    fm <- lmer (frm, roi.df)

    banner("Normality test for residues of full model","=")
    show(shapiro.test(resid(fm)))
    pdf (file.path (assumptions.fig.dir, sprintf ("%s-QQfull.pdf", var)))
    qqmath(fm)
    dummy <- dev.off ()

    pdf (file.path (assumptions.fig.dir, sprintf ("%s-normal.pdf", var)))
    plot_model(fm, type='diag')[[3]]
    dummy <- dev.off ()    

    banner("Normality test for random effect","=")
    r=ranef(fm)$subject$`(Intercept)`
    show(shapiro.test(r))
    pdf (file.path (assumptions.fig.dir, sprintf ("%s-QQrandom.pdf", var)))
    qqnorm(r)
    qqline(r)
    dummy <- dev.off ()

    pdf (file.path (assumptions.fig.dir, sprintf ("%s-variance.pdf", var)))
    plot_model(fm, type='diag')[[4]]
    dummy <- dev.off ()
}

### * Close the output text file
sink.close ()
