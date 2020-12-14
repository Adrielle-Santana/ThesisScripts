### Analysis of N1/P2 complex

### * Load libraries

### ** Load local library
source ("paths.r")
source ("sink.r")
source ("contrasts.r")

### ** Load system libraries
load.pkgs (c ("lme4", "lmerTest", "emmeans"))

### * N1/P2 data

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
results.filestem <- "N1-P2-analysis"
sink.open (file.path (results.dir, sprintf ("%s.txt", results.filestem)))

### * Build contrasts

### ** Generate polynomial contrasts for stimulus factor
##poly <- contr.poly (5)

### ** One-way contrasts
contr <- list ("feature" = build.contr (n1.p2.df$feature),
               "type" = build.contr (n1.p2.df$type),
               "electrode" = list ("medial - lateral" = c (-1/4, -1/4, 1, -1/4, -1/4),
                                   "frontal - temporal" = c (1/2, 1/2, 0, -1/2, -1/2),
                                   "left - right" = c (1/2, -1/2, 0, -1/2, 1/2)),
               "stimulus" = list ("linear" = c(-1, -1/2, 0, 1/2, 1),
                                  "phy - psy" = c (-1/4, 1/2, 0, -1/2, 1/4),
                                  "ambiguity" = c (1/4, 1/4, -1, 1/4, 1/4)))

               ##"stimulus" = list ("linear" = poly [, 1],
               ##                   "quadratic" = poly [, 2],
               ##                   "cubic" = poly [, 3],
               ##                   "fourth" = poly [, 4]))

### ** Two-way contrasts
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


### ** Check contrasts
banner ("Check contrasts", "*")
for (n in names (contr)) {
    banner (n, "=")
    check.contrast (contr [[n]])
}

### * Directory for saving plots
n1.p2.fig.dir <- file.path (figures.dir, results.filestem)
unlink (n1.p2.fig.dir, recursive = TRUE)
force.dir.create (n1.p2.fig.dir)

### * Dependent variables analysis

### ** Get variable names
### We should get here:
### T1, N1, T2, P2, T2.T1, P2.N1, T1.boot, N1.boot, T2.boot, P2.boot,
### T2.T1.boot, and P2.N1.boot
dep.vars <- names (n1.p2.df) [seq (6, 17)]

### ** Loop over dependent variables
for (var in dep.vars) {

    ## *** Banner for dependent variable
    banner (sprintf ("Dependent variable: %s", var), "*")

    ## *** Model formula for the current dependent variable
    frm <- as.formula (sprintf ("%s  ~ feature * type * electrode * stimulus + (1 | subject)", var))
    
    ## *** Fit general model
    fm <- lmer (frm, n1.p2.df)

    ## *** Global ANOVA table
    banner ("ANOVA table for full model", "=")
    show (anova (fm))

    ## *** Plot the original data
    pdf (file.path (n1.p2.fig.dir, sprintf ("%s-0full.pdf", var)))
    print (emmip (fm, electrode ~ stimulus | feature * type,
                  CIs = TRUE, ylab = var))
    dummy <- dev.off ()

    ## *** ISI correction
    if (var == "P2.N1" | var == "P2.N1.boot") {

        ## **** Compensate for different ISI between types of experiment
        ## (Ativo × Passivo)

        ## **** Compute the multiplicative factor between Active and Passive
        ## Do it for each feature/subject/electrode/stimulus.
        ## The data frame below contains the ratio Active/Passive for each case.
        frm.scale <- as.formula (sprintf ("%s  ~ subject * feature * stimulus * electrode", var))
        df <- aggregate (frm.scale, n1.p2.df, function (x) x [2] / x [1])

        ## **** Get the median value and plot the values
        act.pass.factor <- median (df [[var]], na.rm = TRUE)
        plot (df [[var]], pch = 19, col = "#00000080")
        abline (h = act.pass.factor, col = "#ff000080", lwd = 3)
        plot (density (df [[var]], na.rm = TRUE),
              main = "",
              xlab = "active-passive ratio for feature/subject/electrode/stimulus")
        axis (1, at = act.pass.factor, labels = sprintf ("%.2f", act.pass.factor))
        abline (v = act.pass.factor, col = "red")

        ## **** Correction factor for the Ativo cases
        idx <- which (n1.p2.df$type == "Active")
        n1.p2.df [[var]] [idx] <- n1.p2.df [[var]] [idx] / act.pass.factor

        ## **** Fit the model with “ISI correction”
        ## There should be no effect now for factor type
        fm <- lmer (frm, n1.p2.df)

         ## *** Global ANOVA table
        banner ("ANOVA table for full corrected  model", "=")
        show (anova (fm))

        ## *** Plot the original data
        pdf (file.path (n1.p2.fig.dir, sprintf ("%s-0full-Corrected.pdf", var)))
        print (emmip (fm, electrode ~ stimulus | feature * type,
                      CIs = TRUE, ylab = var))
        dummy <- dev.off ()

    }

    ## *** Stepwise selection
    step.val <- step (fm)
    fm <- get_model (step.val)

    ## *** Random effect
    banner ("Random effect ANOVA table for step-wise selected model", "=")
    show (ranova (fm))

    ## *** ANOVA table
    banner ("ANOVA table for step-wise selected model", "=")
    aov <- anova (fm)
    show (aov)

    ## *** Check effects

    ## *** Get names of the variables in the model
    n <- names (attributes (aov) $ hypotheses)

    ## *** Loop over effects
    for (i in seq (1, length (n)))

        ## **** Only consider factors with significant effects
        if (aov [["Pr(>F)"]] [i] < 0.05) {

            ## ***** Get factor name
            eff <- n [[i]]

            ## ***** Only consider factors for which effects have been defined
            if (eff %in% names (contr)) {

                ## ****** Banner for effect
                banner (sprintf ("Effect: %s", eff), "=")

                ## ****** Get emmeans formula
                frm.em <- as.formula (sprintf ("~ %s", gsub (":", " * ", eff)))

                ## ****** Test the effect
                em <- emmeans (fm, frm.em)
                show (summary (em))
                show (contrast (em, method = contr [[eff]]))

                ## ****** Plot the effect
                if (grepl (":", eff))
                    frm.emmip <- as.formula (gsub (":", " ~ ", eff))
                else
                    frm.emmip <- as.formula (sprintf ("~ %s", eff))
                pdf (file.path (n1.p2.fig.dir, sprintf ("%s-%s.pdf", var, eff)), width = 5, height = 4)
                print (emmip (fm, frm.emmip, ylab = var, CIs = TRUE))
                dummy <- dev.off ()

            }

        }

}

### * Close the output text file
sink.close ()

### * Compose the PDF file with the results for all subjects
system (sprintf ("pdftk %s/*.pdf cat output %s/%s.pdf",
                 n1.p2.fig.dir, figures.dir, results.filestem))

