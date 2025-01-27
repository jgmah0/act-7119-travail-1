###
### Travail 1, ACT-7119
### Utilitaires
###
##

# Fonction indicatrice "geq"
calculer_indicatrice_geq <- function(xinf, xsup)
{
    xinf <= xsup
}

# Fonction pour des lois de probabilité
qpareto <- function(x, al, la)
{
    # ...
}

# u = 1 - (la / (la + x))^al
# (1 - u)^(1 / al) = la / (la + x)
# x = la (1 - (1 - u)^(1 / al)) / ((1 - u)^(1 / al))

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
