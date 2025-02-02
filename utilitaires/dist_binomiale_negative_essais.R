###
### Loi binomiale négative essais
### Source des formules : Notes de cours du chapitre 4 du cours ACT-1002
###                       donné à l'Université Laval.
###
##


dnbinom_essais <- function(k, n, q)
{
    choose(k - 1, n - 1) * (q^n) * ((1 - q)^(k - n))
}

