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


# TLS inverse pour le modèle avec copule Archimédienne hiérarchique
# geom-geom
fgp_inv_geom_essais <- function(u, q)
{
    u / (q + (1 - q) * u)
}

tls_inv_geom_essais <- function(u, q)
{
    log(( q + u * (1 - q) ) / u)
}

tls_inv_geom_essais_comp_geom_essais <- function(u, q0, q1)
{
    tls_inv_geom_essais(fgp_inv_geom_essais(u, q0), q1)
}

# TLS pour le modèle avec copule Archimédienne hiérarchique
# geom-geom
tls_geom_essais <- function(t, q)
{
    (q * exp(-t)) / (1 - (1 - q) * exp(-t))
}

fgp_geom_essais <- function(t, q)
{
    (q * t) / (1 - (1 - q) * t)
}


### Discrétiser
discr <- function(h, pDistX, qDistX,
                    tol = 10^-7, method = "lower")
{
    max = qDistX(1 - tol)

    length <- max / h

    if(method == "lower")
    {
        length <- length + 1
        return(c(0,  pDistX((1:length) * h) -
              pDistX((0:(length - 1)) * h)))
    }

    c(pDistX(h),
      pDistX((2:length) * h) - pDistX((1:(length - 1)) * h))
}


VaR_kapp_discr <- function(kapp, h, fx)
{
  (min(which(cumsum(fx) >= kapp) - 1)) * h
}
# fxx <- discr(0.01, function(x) plnorm(x, log(10) - 0.32, 0.8),
#              function(x) qlnorm(x, log(10) - 0.32, 0.8),
#              method = "lower")
#
# F_xx <- cumsum(fxx)
#
#
#
# VaR_kapp_methode(0.9999, 0.01, fxx)

TVaR_kapp_discr <- function(kapp, h, fx)
{
  1 / (1 - kapp) *
    sum(pmax(seq(length(fx)) - 1 - VaR_kapp_discr(kapp, h, fx), 0) * fx) +
    VaR_kapp_discr(kapp, h, fx)
}

Mesure_entropique_discr <- function(rho, fx)
{
  1/rho * log(sum(exp(rho * (seq(length(fx)) - 1))*fx))
}

