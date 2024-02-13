module PolyData

using Oscar

# this is duplicated in another module... fix this later
exponent_vectors(poly) = leading_exponent_vector.(terms(poly))

function projective_space_fan_string(n)
  
  projstring = "[[" * join(fill(-1,n-1)," ") * "]"

  for i in n-1:-1:1
    coordarray = zeros(Int,n-1)
    coordarray[i] = 1
    
    projstring = projstring * "[" * join(coordarray, " ") * "]"
  end

  projstring = projstring * "]"

  projstring
end#function

function tcr_string(p,poly,name)
  tcrstring = name * ":"

  n = length(gens(parent(poly)))

  # put the exponent vectors here
  tcrstring = tcrstring * "["

  for exp_vec in exponent_vectors(poly)
    tcrstring = tcrstring * "[" * join(exp_vec[1:end-1], " ") * "]"
  end

  tcrstring = tcrstring * "]:"

  # put the coefficients here
  
  tcrstring = tcrstring * "[" * join(coefficients(poly)," ") * "]:"

  # The fan for projective space

  tcrstring = tcrstring * projective_space_fan_string(n) * ":"

  # the degree here

  degarr = zeros(Int,n)
  deg = sum(exponent_vectors(poly)[1]) # we assume homogeneous
  degarr[1] = deg
  
  tcrstring = tcrstring * "[" * join(degarr," ") * "]"

  # the prime number

  tcrstring = tcrstring * ":" * string(p)

  # println(tcrstring)

  tcrstring
end#function

function output_tcr_file(ps,polys,names,outputfilename)
  open(outputfilename,"w") do outfile
    for i in 1:length(ps)
      println(outfile,tcr_string(ps[i],polys[i],names[i]))
    end
  end
end#function

end#module
