source ("paths.r")
source ("dwt-lib.r")
source ("experiments.r")

load.pkgs ("boot")

electrodes <- read.table (file.path (data.dir, "electrodes.tab")) [, 1]
selected.electrodes <- sapply (c ("Fz", "F7", "F8", "TP9", "TP10"),
                               function (x) which (electrodes == x))

n1.p2.dir <- file.path (figures.dir, "N1-P2")
dir.create (n1.p2.dir, showWarnings = FALSE)

time.basis <- seq (0, by = 1 / eeg.sampfreq, length.out = dwt.length)

idx.lpf <- seq (1, 2016)
n1.region <- seq (round (0.05 * eeg.sampfreq), round (0.17 * eeg.sampfreq))
p2.region <- seq (round (0.10 * eeg.sampfreq), round (0.26 * eeg.sampfreq))

boot.statistic <- function (x, i) {
    m <- colMeans (x [i, ])
    m [idx.lpf] <- 0
    sig <- vec.to.signal (m)
    i1 <- n1.region [1] - 1 + which.min (sig [n1.region])
    n1 <- sig [i1]
    t1 <- time.basis [i1]
    i2 <- p2.region [1] - 1 + which.max (sig [p2.region])
    p2 <- sig [i2]
    t2 <- time.basis [i2]
    n1.p2 <- p2 - n1
    t1.t2 <- t2 - t1
    return (c (t1, n1, t2, p2, t1.t2, n1.p2))
}

boot.R <- 200

stim.cols <- c ("red", "orange", "gray", "green3", "blue")

stat.names <- c ("T1", "N1", "T2", "P2", "T2-T1", "P2-N1")
unit.names <- rep (c ("s", "uV"), 3)

n1.p2.df <- data.frame ()

for (ef in experiment.features)

    for (et in experiment.types)

        for (subj in cohort) {

            load (dwt.coefs.filename (ef, et, subj))

            for (elt in selected.electrodes) {

                n1 <- n1.0 <- p2 <- p2.0 <- n1.p2 <- n1.p2.0 <- c ()

                elt.stem <- sprintf ("%s-%s-S%02d-%s",
                                     ef, et, subj, electrodes [elt])

                elt.t <- elt.t0 <- list ()

                for (stim in seq (1, 5)) {

                    idx <- which (dwt.coefs$stimulus == stim)
                    cf <- dwt.coefs$response [[elt]] [idx, ]

                    b <- boot (cf, boot.statistic, boot.R)

                    file.stem <- sprintf ("%s-stim%d", elt.stem, stim)
                    file.name <- sprintf ("pics-boot-%s.pdf", file.stem)
                    pdf (file = file.path (n1.p2.dir, file.name),
                         width = 10, height = 6)

                    m <- colMeans (cf)
                    plot (time.basis, vec.to.signal (m), type = "l", las = 1,
                          col = "#00000040", xlab = "time (s)",
                          ylab = "amplitude (uV)", main = file.stem,
                          ylim = c (min (b$t [, 2]), max (b$t [, 4])))
                    m [idx.lpf] <- 0
                    lines (time.basis, vec.to.signal (m))

                    points (b$t [, 1], b$t [, 2], col = "#ff808040", pch = 19,
                            cex = 0.5)
                    points (b$t [, 3], b$t [, 4], col = "#8080ff40", pch = 19,
                            cex = 0.5)
                    points (b$t0 [1], b$t0 [2], pch = 21, cex = 2, bg = "red")
                    points (b$t0 [3], b$t0 [4], pch = 21, cex = 2, bg = "blue")

                    dummy <- dev.off ()

                    cat (sprintf ("\r%60s", file.name))
                    flush (stdout ())

                    elt.t0 [[stim]] <- b$t0
                    elt.t [[stim]] <- b$t

                    n <- length (stat.names)
                    est.boot <- rep (NA, )
                    for (i in seq (1, n))
                        if (i %% 2 == 1)
                            ## Temporal measurement
                            est.boot [i] <- mean (b$t [, i])
                        else {
                            ## Pic amplitude measurement
                            d <- density (b$t [, i])
                            est.boot [i] <- d$x [which.max (d$y)]
                        }

                    n1.p2.df <- rbind (n1.p2.df,
                                       data.frame (feature = ef,
                                                   type = et,
                                                   subject = subj,
                                                   electrode = electrodes [elt],
                                                   stimulus = stim,
                                                   T1 = b$t0 [1],
                                                   N1 = b$t0 [2],
                                                   T2 = b$t0 [3],
                                                   P2 = b$t0 [4],
                                                   T2.T1 = b$t0 [5],
                                                   P2.N1 = b$t0 [6],
                                                   T1.boot = est.boot [1],
                                                   N1.boot = est.boot [2],
                                                   T2.boot = est.boot [3],
                                                   P2.boot = est.boot [4],
                                                   T2.T1.boot = est.boot [5],
                                                   P2.N1.boot = est.boot [6]))

                }

                for (i in seq (1, length (stat.names))) {

                    file.name <- file.path (n1.p2.dir, sprintf ("%s-dst-%s.pdf",
                                                                stat.names [i],
                                                                elt.stem))
                    pdf (file = file.name, width = 10, height = 6)

                    x.min <- Inf
                    x.max <- y.max <- -Inf
                    for (j in seq (1, 5)) {
                        v <- elt.t [[j]] [, i]
                        d <- density (v)$y
                        x.min <- min (x.min, min (v))
                        x.max <- max (x.max, max (v))
                        y.max <- max (y.max, max (d))
                    }

                    for (j in seq (1, 5)) {
                        v <- elt.t [[j]] [, i]
                        v0 <- elt.t0 [[j]] [i]
                        if (j == 1)
                            plot (density (v), col = stim.cols [1], lwd = 2, las = 1,
                                  xlim = c (x.min, x.max), ylim = c (0, y.max),
                                  xlab = sprintf ("%s (%s)", stat.names [i],
                                                  unit.names [i]),
                                  ylab = "density", main = elt.stem)
                        else
                            lines (density (v), col = stim.cols [j], lwd = 3)
                        abline (v = v0, col = stim.cols [j], lwd = 2, lty = 2)
                    }

                    legend ("topright", ins = 0.05, col = stim.cols, lwd = 2,
                            legend = sapply (seq (1, 5),
                                             function (x)
                                                 sprintf ("stim #%d", x)))

                    dummy <- dev.off ()

                    cat (sprintf ("\r%60s", file.name))
                    flush (stdout ())

                }

            }

        }

cat ("\n")
flush (stdout ())

n1.p2.df$electrode <- droplevels (n1.p2.df$electrode)

save (n1.p2.df, file = file.path (results.dir, "N1-P2.dat"))

system (sprintf ("pdftk %s/pics-boot-*S*.pdf cat output %s/pics-boot-all.pdf",
                 n1.p2.dir, figures.dir))

for (i in seq (1, length (stat.names)))
    system (sprintf ("pdftk %s/%s-dst-*S*.pdf cat output %s/%s-dst-all.pdf",
                     n1.p2.dir, stat.names [i], figures.dir, stat.names [i]))
