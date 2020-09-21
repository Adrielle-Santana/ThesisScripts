### * Utility functions for generating contrasts for multiple comparison analysis

### * Function for checking contrast orthogonality
check.contrast <- function (contr) {
    k <- length (contr)
    n <- names (contr)
    for (i in seq (1, k)) {
        cat (sprintf ("sum (#%d) = %f\n", i, sum (contr [[i]])))
        if (i < k)
            for (j in seq (i + 1, k))
                cat (sprintf ("sum (#%d * #%d) = %f\n",
                              i, j, sum (contr [[i]] * contr [[j]])))
    }
}

### * Function for building contrast for one-way effetcs with two levels
build.contr <- function (var) {
    lev <- levels (var)
    name <- sprintf ("%s - %s", lev [1], lev [2])
    retval <- list ()
    retval [[name]] = c (-1, 1)
    return (retval)
}

### * Function for building contrast for two-way effetcs
prod.contr <- function (c1, c2) {
    retval <- c ()
    for (i in seq (1, length (c2)))
        retval <- c (retval, c1 * c2 [i])
    return (retval)
}
