function hw_algorithm(p, n ,f)
	d = total_degree(f)
	vectors = collect(compositions(d, n))
	g = f^(p-1)
	collect_coeff = []

#collecting the desired coefficients

	for i in eachindex(vectors)
		push!(collect_coeff, [])
		for j in eachindex(vectors)
			push!(collect_coeff[i], coeff(g,p*vectors[i] - vectors[j]))
		end
	end

	return reduce(hcat, collect_coeff)

end
