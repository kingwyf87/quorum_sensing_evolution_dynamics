# random number drawn from Poisson distribution (Knuth)
function randNum_Poisson(lambda::Float64)

	L = exp(-lambda)
	k = 0
	p = 1.

	k = k+1
	u = rand()
	p = p*u

	while p > L
	  k = k+1
	  u = rand()
	  p = p*u
	end

	return k-1

end
