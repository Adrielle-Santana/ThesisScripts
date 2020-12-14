### Analyses of regression directions obtained from RoLDSIS

### * Load local libraries
source ("paths.r")
source ("dwt-lib.r")
source ("experiments.r")
source ("sink.r")
source ("contrasts.r")

### * Load system libraries
load.pkgs (c ("lme4", "lmerTest", "emmeans", "Cairo", "ggplot2", "wavelets"))

### Directory to save figures
roldsis.dir <- file.path (figures.dir, "Roldsis-analysis")
force.dir.create (roldsis.dir, clear = TRUE)

### * Initialize data frames
roi.df <- data.frame ()

### * Electrode names
electrodes <- c ("F7", "F8", "Fz", "TP9", "TP10")

### * Load previously computed direction vectors (by RoLDSIS)
load (file = file.path (results.dir, sprintf ("dirVec-%d.dat", eeg.sampfreq)))

discrepancy.angle <- function (x, y) {
    nx <- x / sqrt (sum (x ^ 2))
    ny <- y / sqrt (sum (y ^ 2))
    acos (sum (nx * ny))
}

discrepancy.distance <- function (x, y)
    sqrt (sum ((x - y) ^ 2))

discrepancy.distance.ang <- function (x,y){
    nx <- x / sqrt (sum (x ^ 2))
    ny <- y / sqrt (sum (y ^ 2))
    (acos (sum (nx * ny)) * 180 / pi) * discrepancy.distance (x, y)
}

discrepancy.diff.norm <- function (x, y)
    sqrt (sum (x ^ 2)) - sqrt (sum (y ^ 2))

discrepancy.diff.abs <- function (x, y)
    sum (abs (x)) - sum (abs (y))

discrepancy.diff.norm.times.distance <- function (x, y)
    discrepancy.diff.norm (x, y) * discrepancy.distance (x, y)

discrepancy <- discrepancy.distance

### * Names of the ROI
roi.names <- c ("r1.early.theta", "r2.late.theta",
                "r3.early.alpha", "r4.late.alpha",
                "r5.early.beta", "r6.late.beta",
                "r7.early.gamma", "r8.late.gamma")

### * Fill data frame
for (ef in experiment.features) {

    efAux <- ifelse (ef == "Formantes", "Formants", "VOT")

    for (et in experiment.types) {

        etAux <- ifelse (et == "Passivo", "Passive", "Active")

        for (subj in cohort) {

            S <- sprintf ("subj%d", subj)

            for (elt in electrodes) {

                ## ** Select direction vector and convert to wavelet object
                ## aux <- dirVec [[ef]] [[et]] [[S]] [[reg]] [[elt]]
                ## obj <- vec.to.dwt(aux, n=dwt.length)

                aux.phy <- dirVec [[ef]] [[et]] [[S]] $phy [[elt]]
                obj.phy <- vec.to.dwt (aux.phy, n = dwt.length)

                aux.psy <- dirVec [[ef]] [[et]] [[S]] $psy [[elt]]
                obj.psy <- vec.to.dwt (aux.psy, n = dwt.length)

                ## ** Store data in the data frame 1
                tmp.df <- data.frame (subject = subj,
                                      electrode = elt,
                                      feature = efAux,
                                      type = etAux)

                ## ** Define ROIs

                roi.val <- c (discrepancy (obj.phy@W [[9]] [1],
                                           obj.psy@W [[9]] [1]),
                              discrepancy (obj.phy@W [[9]] [2],
                                           obj.psy@W [[9]] [2]),
                              discrepancy (obj.phy@W [[8]] [1:2],
                                           obj.psy@W [[8]] [1:2]),
                              discrepancy (obj.phy@W [[8]] [3:4],
                                           obj.psy@W [[8]] [3:4]),
                              discrepancy (obj.phy@W [[7]] [1:4],
                                           obj.psy@W [[7]] [1:4]),
                              discrepancy (obj.phy@W [[7]] [5:8],
                                           obj.psy@W [[7]] [5:8]),
                              discrepancy (obj.phy@W [[6]] [1:8],
                                           obj.psy@W [[6]] [1:8])
                                  + discrepancy (obj.phy@W [[5]] [1:16],
                                                 obj.psy@W [[5]] [1:16]),
                              discrepancy (obj.phy@W [[6]] [9:16],
                                           obj.psy@W [[6]] [9:16])
                                  + discrepancy (obj.phy@W [[5]] [17:32],
                                                 obj.psy@W [[5]] [17:32]))

                for (i in seq (1, length (roi.names)))
                    tmp.df [[roi.names [i]]] <- roi.val [i]

                roi.df <- rbind (roi.df, tmp.df)

                cat (sprintf ("\rFeature: %s  Type: %s  Subject: %s  Electrode: %s",
                              ef, et, subj, elt))
                flush (stdout ())

                aux2 <- 1

            } ## elt
        } ## subj
    } ## et
} ## ef

cat ("\n")
flush (stdout ())

roi.df$electrode <- as.factor (roi.df$electrode)
roi.df$feature <- as.factor (roi.df$feature)
roi.df$type <- as.factor (roi.df$type)

### One-way contrasts
contr <- list ("feature" = build.contr (roi.df$feature),
               "type" = build.contr (roi.df$type),
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

results.filestem <- "energy-analysis"
sink.open (file.path (results.dir, sprintf ("%s.txt", results.filestem)))
sink (type = "message")

### Check contrasts
banner ("Check contrasts", "*")
for (n in names (contr)) {
    banner (n, "=")
    check.contrast (contr [[n]])
}

### * Maximum discrepancy (for y axis of emmip plots)
min.discrepancy <- 0
max.discrepancy <- 0.11

### * Aux function for prettifying VD name
pretty.name <- function (name) {
    name <- sub ("^r.\\.", "", name)
    name <- sub ("\\.theta", " θ", name)
    name <- sub ("\\.alpha", " α", name)
    name <- sub ("\\.beta", " β", name)
    name <- sub ("\\.gamma", " γ", name)
    return (name)
}

### * Compute mixed effects models for both data frames

for (var in roi.names) {

    banner (var, "*")

    frm <- as.formula (sprintf ("%s ~ feature * type * electrode + (1 | subject)",
                                 var))

    fm <- lmer (frm, roi.df)

     ## *** Global ANOVA table
    banner ("ANOVA table for full model", "=")
    show (anova (fm))

    ## *** Plot the original data
    CairoPDF (file.path (roldsis.dir, sprintf ("%s-0full.pdf", var)))
    print (emmip (fm, ~ electrode | feature * type,
                  CIs = TRUE, ylab = pretty.name (var)))
    dummy <- dev.off ()

    ## Treatment in case there is no significant effect
    if (all (anova (fm)$Pr > 0.05))
        next

    ## Stepwise selection
    step.val <- step (fm)
    fm.step <- get_model (step.val)

    ## Check return of step function:
    ##    if the random effect is significant: returns class 'lmerModLmerTest'
    ##    if the random effect is not significant: returns a list with 13 elements
    if (length (fm.step) > 1) {
        terms <- paste (attr (fm.step$terms, "term.labels"),
                        sep = " ", collapse = "+")
        frm <- as.formula (sprintf ("%s ~ %s + (1 | subject)", var, terms))
        fm.step <- lmer (frm, roi.df)
    }

    ## Random effect
    banner ("Random effect ANOVA table for step-wise selected model", "=")
    show (ranova (fm.step))

    ## ANOVA table
    banner ("ANOVA table for step-wise selected model", "=")
    aov <- anova (fm.step)
    show (aov)

    ## Check effects

    ##  Get names of the variables in the model
    n <- names (attributes (aov) $ hypotheses)

    ## Skip case with no hypotheses
    if (length (n) == 0)
        next

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
                em <- emmeans (fm.step, frm.em)
                show (summary (em))
                show (contrast (em, method = contr [[eff]]))

                ## Plot the effect
                if (grepl (":", eff)){
                    frm.emmip <- as.formula (gsub (":", " ~ ", eff))
                } else{
                    frm.emmip <- as.formula (sprintf ("~ %s", eff))
                }

                CairoPDF (file.path (roldsis.dir,
                                     sprintf ("%s-%s.pdf", var, eff)),
                          width = 5, height = 4)
                print (emmip (fm.step, frm.emmip, ylab = pretty.name (var),
                              CIs = TRUE)
                       + ylim (min.discrepancy, max.discrepancy))
                dummy <- dev.off ()

                message (sprintf ("ROI: %s  factor: %s", var, eff))

            } ## if eff
        } ## if aov
    } ## for i
} ## for var

sink.close ()

##  Compose the PDF file with the results for all ROIs
system (sprintf ("pdftk %s/*.pdf cat output %s/roldsis-analysis.pdf",
                 roldsis.dir, figures.dir))
