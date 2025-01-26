###
### Travail 1, ACT-7119
### Algorithme 13 de [Cossette et al., 2019]
###
##


# À tester.
echantillonner_crm_cop_archi_hiera <- function(n, qDistTheta0, tlsTheta0, qDistB, tlsTheta1, qDistN, qDistX)
{
    realisations <- list()
    realisation_i <- vector(mode = "numeric")

    # Simuler les réalisation de N
    realisations_N <- qDistN(tlsTheta0(rexp(n) / qDistTheta0(runif(n))))

    R_X <- 0
    U_X <- 0
    Theta01 <- 0

    max_N_i <- 0

    for (i in seq(n))
    {
        # Stocker la réalisation i de N
        realisation_i[1] <- realisations_N[i]

        if (realisation_i[1] != 0)
        {
            for (j in seq(realisations_N[i]))
            {
                Theta01 <- sum(qDistB(runif(realisations_N[i])))
                R_X <- rexp(1)
                U_X <- tlsTheta1(R_X /Theta01)

                realisation_i <- c(realisation_i, qDistX(U_X))
            }
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
