# calculate signal concentration at equilibrium state
#   syms p N u S K r
#   eqn = p(j).*N(i)-u.*S + p(j).*N(i).*r.*S^1/(K^1+S^1) == 0;
#   sol(i,j) = mean(double(vpasolve(eqn,S,[0 Inf]))); /solx = solve(eqn,S)
function cal_Signal_Concentration(p::Matrix{Float64}, N::Matrix{Float64}, u::Float64, r::Matrix{Float64}, K::Float64)
	
	return ((N.^2.*p.^2.*r.^2 + 2.*N.^2.*p.^2.*r + N.^2.*p.^2 - 2.*N.*K.*p.*r.*u + 2.*N.*K.*p.*u + K^2.*u^2).^(1/2) + N.*p - K.*u + N.*p.*r)/(2.*u)

end

