# mutation operation for evolving traits 
function mut_parameter(mut_Vector::Vector{Float64},mut_P::Float64,mut_SD::Float64,mut_Min::Float64,mut_Max::Float64,index_Cheats::Vector{Int64},size_Pop::Int64)

	num_Mut = randNum_Poisson(mut_P*size_Pop)
	if num_Mut!=0
	    index_Mut = randperm(size_Pop)[1:num_Mut]
	    index_Diff = index_Mut[index_Cheats[index_Mut].==0]
	    if !isempty(index_Diff)
	        mut_Vector[index_Diff] = mut_Vector[index_Diff] + randn(length(index_Diff)).*mut_SD
	        mut_Vector[index_Diff[find(mut_Vector[index_Diff].<mut_Min)]] = mut_Min
	        mut_Vector[index_Diff[find(mut_Vector[index_Diff].>mut_Max)]] = mut_Max
	    end
	end

end
