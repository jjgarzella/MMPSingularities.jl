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

# MARK - Inputting data from census of cubic fourfolds in char=2



function read_orbit_file_to_byte_array(filename)
  s = open(filename,"r")
  data = UInt8[]
  data = readbytes!(s,data,typemax(Int))
  close(s)
  data
end#function

function read_cubic_fourfold_orbit_reps(dirname)
  result = []
  R, vars = polynomial_ring(GF(2),6)
  monomials = [] # put the monomials in order here
  # read the files

  for i in 1:85

    println("Starting file $i")
    
    # convert each one to array of cubics
    bytes = read_orbit_file_to_byte_array(dirname * "/orbitreps-$i.data")

    l = length(bytes)
    l % 7 != 0 && println("Number of bytes $l is not a multiple of 7... " * string(l % 7) * " mod 7")

    println("Found file with $(div(l,7)) cubic fourfolds")

    for j in 1:7:l

      sevenbytes = bytes[j:j+6]
      bits = []
      for k in 0:6
        # extract the bits in a fancy way
        single_byte_of_bits = [sevenbytes[k+1] & (0x1<<n) != 0 for n in 0:7]

        # in case of emergency, uncomment the following line of code
        # reverse(single_byte_of_bits)
        # ....er, I mean endianness, not emergency

        append!(bits, single_byte_of_bits)
      end

      poly = GF(2).(bits) .* monomials 

      # append to result
      push!(result,poly)
    end
  end 
  result
end#function

end#module
