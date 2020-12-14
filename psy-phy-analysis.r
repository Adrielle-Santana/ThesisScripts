source ("paths.r")
source ("parameters.r")
source ("dwt-lib.r")
source ("scalogram.r")
source ("experiments.r")

### * Load the system packages
load.pkgs (c ("shape", "Cairo"))

load (file = file.path (results.dir, sprintf("dirVec-%d.dat",eeg.sampfreq)))

load (file = file.path (results.dir, "id-slope.dat"))

angle <- list ()

fig.roldsis.dir <- file.path (figures.dir, "roldsis")
force.dir.create (fig.roldsis.dir)

nb.wavelets <- 2 * length (dwt (rep (0, dwt.length))@W [[dwt.start.level]])

plot.roldsis.dir <- function (dir, palette, limits, main) {
    layout (matrix (c (2, 1, 3, 4, 6, 5), nrow = 2, byrow = TRUE))
    for (i in seq (1, 5)) {
        obj <- vec.to.dwt (dir [((i - 1) * nb.wavelets + 1) : (i * nb.wavelets)],
                           n = dwt.length)
        plot.scalogram (obj, palette = palette, limits = limits)
    }
    par (mar = c (0, 0, 5, 0))
    plot (0, 0, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
          bty = "n", main = main)
}

cont <- 0
for (ef in experiment.features) {

    angle [[ef]] <- list ()
    cont <- cont+1

    for (et in experiment.types) {

        psy <- c ()
        for (subj in cohort) {
            x <- c ()
            for (i in seq (1, 5))
                x <- c (x, dirVec [[ef]] [[et]] [[subj]] $ psy [[i]])
            psy <- rbind (psy, x)
            pdf (file = file.path (fig.roldsis.dir,
                                   sprintf ("%s-%s-s%02d-phy.pdf",
                                            ef, et, subj)))
            plot.roldsis.dir (x, palette.bwr, NA,
                              sprintf ("%s %s S%d Phy", ef, et, subj))
            dummy <- dev.off ()
        }

        phy <- c ()
        for (subj in cohort) {
            x <- c ()
            for (i in seq (1, 5))
                x <- c (x, dirVec [[ef]] [[et]] [[subj]] $ phy [[i]])
            phy <- rbind (phy, x)
            pdf (file = file.path (fig.roldsis.dir,
                                   sprintf ("%s-%s-s%02d-psy.pdf",
                                            ef, et, subj)))
            plot.roldsis.dir (x, palette.bwr, NA,
                              sprintf ("%s %s S%d Psy", ef, et, subj))
            dummy <- dev.off ()
        }

        psy.hist <- sqrt (colSums (psy ^ 2))
        phy.hist <- sqrt (colSums (phy ^ 2))
        diff.psy.phy <- psy.hist - phy.hist

        pdf (file = file.path (fig.roldsis.dir,
                               sprintf ("%s-%s-phy.pdf", ef, et)))
        plot.roldsis.dir (phy.hist, palette.wb, c (0, max (phy.hist)),
                          sprintf ("%s %s Phy", ef, et))
        dummy <- dev.off ()

        pdf (file = file.path (fig.roldsis.dir,
                               sprintf ("%s-%s-psy.pdf", ef, et)))
        plot.roldsis.dir (psy.hist, palette.wb, c (0, max (psy.hist)),
                          sprintf ("%s %s Psy", ef, et))
        dummy <- dev.off ()

        pdf (file = file.path (fig.roldsis.dir,
                               sprintf ("%s-%s-diff.pdf", ef, et)))
        plot.roldsis.dir (diff.psy.phy, palette.bwr,
                          c (min (diff.psy.phy), max (diff.psy.phy)),
                          sprintf ("%s %s (Psy - Phy)", ef, et))
        dummy <- dev.off ()

        ang <- rep (NA, length (cohort))
        for (i in cohort)
            ang  [i] <- acos (sum (phy [i, ] * psy [i, ])) * 180 / pi

        beta <- id.slope[cont,]


        pdf (file = file.path (fig.roldsis.dir,
                               sprintf ("%s-%s-slope-angle.pdf", ef, et)))
        par(mar=c(5,5,4,1))
        plot (beta, ang, pch = 19, las = 1, cex=2, cex.main = 2, cex.axis=2, cex.lab=2, main = sprintf("\nr: %g   p.value: %g", round(cor.test(beta,ang)$estimate,3),round(cor.test(beta,ang)$p.value,3)),
               ylab = "phy/psy angle (degrees)", xlab = "slope (%/ms)")
        abline (lm (ang ~ beta))
        dummy <- dev.off ()

        angle [[ef]] [[et]] <- ang

        cat (sprintf ("\r%s %s", ef, et))
        flush (stdout ())

    }

}

cat ("\n")

system (sprintf ("pdftk %s/*-*-*.pdf cat output %s/all.pdf",
                 fig.roldsis.dir, fig.roldsis.dir))
