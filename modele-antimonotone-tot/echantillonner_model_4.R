###
### Travail 1
### Fonctions pour échantillonner un modèle avec antimonotones
### entre les composantes du modèle.
###
##

echantillonner_modele_4 <- function(n, qDistN, qDistX, rseed = 201)
{
    realisations <- list()
    realisation_i <- vector(mode = "numeric")
    U_i <- 0

    max_N_i <- 0

    for (i in seq(n))
    {
        U_i <- runif(1)

        realisation_i <- vector(mode = "numeric")
        realisation_i[1] <- qDistN(U_i)

        if (realisation_i[1] != 0)
        {
            realisation_i <- c(realisation_i[1], rep(qDistX(1 - U_i), realisation_i[1]))
        }

        if (realisation_i[1] > max_N_i)
        {
            max_N_i <- realisation_i[1]
        }

        realisations[[i]] <- realisation_i
    }

    list(realisations = realisations,
         max_N_i = max_N_i)
}

structurer_echantillon <- function(ech_list, max_N_i)
{
    n <- length(ech_list)

    realisations <- matrix(rep(0, n * (max_N_i + 1)), nrow = n, ncol = max_N_i + 1)

    for (i in seq(n))
    {
        realisations[i, 1:(ech_list[[i]][1] + 1)] <- ech_list[[i]]
    }

    realisations
}
