###
### Travail 1, ACT-7119
### Fonction pour approximer la valeur de kmax.
###
### Exemple d'utilisation :
### approximer_kmax(function(x) ppois(x, 4), 0.00001)
###
##

approximer_kmax <- function(pDistN, threshold)
{
    kmax <- 0 # initialisation à 0

    while (pDistN(kmax) <= 1 - threshold)
    {
        kmax <- kmax + 1
    }

    kmax
}

approximer_kmax_version_2 <- function(dDistN, borne_inf_support, threshold)
{
    kmax <- borne_inf_support # initialisation à 0

    while (sum(dDistN(borne_inf_support:kmax)) <= 1 - threshold)
    {
        kmax <- kmax + 1
    }

    kmax
}
