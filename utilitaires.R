###
### Travail 1, ACT-7119
### Utilitaires
###
##

# Fonction pour des lois de probabilité
qpareto <- function(x, al, la)
{
    # ...
}

dlogarithmique <- function(x, ga)
{
    -((ga^x) / (x * log(1 - ga)))
}

qlogarithmique <- function(u, ga)
{
    somme_prob <- 0
    i <- 0
    while (somme_prob < u)
    {
        i <- i + 1
        somme_prob <- somme_prob + dlogarithmique(i, ga)
    }

    i
}
