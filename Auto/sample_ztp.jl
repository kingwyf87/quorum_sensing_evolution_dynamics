# zero-truncated Poisson distribution
function sample_ztp(lambda::Float64)
  k = 1
  t = exp(-lambda) / (1 - exp(-lambda)) * lambda
  s = t
  u = rand()
  while s < u
    k = k + 1
    t = t * lambda / k
    s = s + t
  end
  return k
end
