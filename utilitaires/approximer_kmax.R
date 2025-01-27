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
