###
### Travail 1
### Fonctions pour échantillonner un modèle avec comonotonicité
### entre les composantes du modèle.
### Présentement, des distributions sont fixées pour la
### fréquence et la sévérité, mais il sera possible de
### rendre la fonction plus dynamique selon plusieurs choix
### de distributions de fréquence et de sévérité.
###
##

echantillonner_modele_3 <- function(n, rseed = 201)
{
    realisations <- list()
    realisation_i <- vector(mode = "numeric")
    U_i <- 0

    max_N_i <- 0

    for (i in seq(n))
    {
        U_i <- runif(1)
        realisation_i[1] <- qpois(U_i, 10)

        if (realisation_i[1] != 0)
        {
            realisation_i <- c(realisation_i[1], rep(qexp(U_i, 0.001), realisation_i[1]))
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
