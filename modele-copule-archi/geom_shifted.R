# Construction de la géométrique essais
dgeo_shifted = function(k, alpha)
{
  # (k > 0)*(alpha)^(k-1)*(1-alpha)
  dgeom(k - 1, 1 - alpha)
}

sum(sapply(seq(50) - 1, function(y) dgeo_shifted(y, 0.2))) # Somme à 1

pgeo_shifted = function(k, alpha)
{
  # (k > 0)*(1 - alpha^k)
  pgeom(k - 1, 1 - alpha)
}

qgeo_shifted = function(u, alpha)
{
  qgeom(u, 1 - alpha) + 1
}