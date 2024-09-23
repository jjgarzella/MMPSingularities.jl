#module PolyData
#
#using Oscar
#
#include("griffiths-dwork-construction/Utils.jl")

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

"""
Reads file with file name 'filename'
into a byte array.

This is not actually specific to fourfolds
"""
function read_orbit_file_to_byte_array(filename)
  s = open(filename,"r")
  data = UInt8[]
  readbytes!(s,data,typemax(Int))
  close(s)
  data
end#function

"""
Reads the data from the orbit representatives into 
a large array.

See the function LoadCubicOrbitData in 
https://github.com/JonahWeinbaum/cubic-fourfolds/blob/main/CubicLib/CubicLib.m

if firsttwo == true, then only use the first two files.

looks for files with the Census of Cubic fourfolds paper's
naming convention in the directory provided.
"""
function read_cubic_fourfold_orbit_reps(dirname,firsttwo=false)
  result = []
  R, vars = polynomial_ring(GF(2),6)
  monomials = compute_monomials(6,3,R,vars)
  # read the files

  N = 85
  firsttwo && (N = 2)

  for i in 1:N

    println("Starting file $i")
    
    # convert each one to array of cubics
    bytes = read_orbit_file_to_byte_array(dirname * "/orbitreps-$i.data")
    println(typeof(bytes))

    l = length(bytes)
    l % 7 != 0 && println("Number of bytes $l is not a multiple of 7... " * string(l % 7) * " mod 7")

    println("Found file with $(div(l,7)+1) cubic fourfolds")

    for j in 1:7:l

      sevenbytes = bytes[j:j+6]
      bits = zeros(UInt8,56)
      for k in 0:6
        # extract the bits in a fancy way
        single_byte_of_bits = zeros(UInt8,8)
        for n in 0:7
          single_byte_of_bits[n+1] = sevenbytes[k+1] & (0x1<<n) != 0 
        end

        # in case of emergency, uncomment the following line of code
        # reverse(single_byte_of_bits)
        # ....er, I mean endianness, not emergency

        bits[k*7 + 1:k*7 + 8] = single_byte_of_bits
        #append!(bits, single_byte_of_bits)
      end

      poly = sum(GF(2).(bits) .* monomials)

      poly == zero(R) && println("Found a zero polynomial at (i,j) = ($i,$j)")
      # append to result
      push!(result,poly)
      #println("Finished fourfold $(div(j,7)+1)")
    end
  end 
  result
end#function

"""
Reads the data from the csv of zeta functions of cubic fourfolds
that Jack Petok sent me (JJ). 

See ReadZetaFunctions in
https://github.com/JonahWeinbaum/cubic-fourfolds/blob/main/CubicLib/CubicLib.m

"""
function read_cubic_fourfold_zeta_functions(filename)
  lines = readlines(filename)
  zetas = Dict{Int,Any}()
  counter = 0
  for entry in lines
    components = split(entry,",",limit=2)
    key = parse(Int,components[1])
    withslashes = replace(components[2],"/" => "//")
    coefs = eval(Meta.parse(withslashes))
    zetas[key] = process_cubic_fourfold_zeta(coefs)
    counter = counter + 1
    if counter % 1000 == 0
      println("Read $counter fourfolds")
    end
  end

  zetas
end#function

"""
Processes the zeta functions as the code in 
newton-polygons-survey.m does in 
https://github.com/JonahWeinbaum/cubic-fourfolds/blob/main/CubicLib/CubicLib.m

"""
function process_cubic_fourfold_zeta(coefs)
  newcoefs = zeros(Int128,length(coefs)) # won't work for very very very large numbers,
  # ...but cubic fourfolds have L-polynomials of degree about 20
                 

  # evaluate the given polynmial at 4t, i.e. gives f(4t)
  for i in eachindex(coefs)
    c = 4^(i-1)
    newcoefs[i] = convert(Int,c*coefs[i])
  end

  # then check if first term is -1, otherwise swap the signs
  #
  # To be honest, I'm not exactly sure why the census does this...
  # it can't hurt though, so I'll do it just in case.
  if newcoefs[1] == -1
    newcoefs .* -1
  end
  newcoefs
end

# MARK - newton polygons

"""
returns the p-adic valuation of the integer num
"""
function padic_val(p,num)
  num == zero(num) && (return typemax(Int128)) # this may cause errors with big numbers
  fact = factor(num)
  if p in fact
    return fact[p]
  end

  # we are a p-adic unit!
  0
end#function

"""
Given a list of coefficients of a polynomial,
starting with the constant term and ending with the
highest degree term, gives the newton polygon
using a very naive algorithm.

This should work for newton polygons of varieties,
but I don't know exactly what I'm assuming about the
points here.

coefs is an array of coeffiecients,
such that the constant term is in coefs[1]
and the coefficients are integers
"""
function newton_polygon(p,coefs)
  vertices = [(0,0)]

  points = [(i-1,padic_val(2,coefs[i])) for i in eachindex(coefs)]

  x = 0
  while x < length(coefs) - 1 # i.e. the end of the hodge polygon
    #println("x=$x")
    startpoint = points[x+1]
    x2 = x + 1
    smallest_slope = typemax(Int128)
    finalx = x
    while x2 ≤ length(coefs) - 1
      #println("x2=$x2, y1=$(points[x2 + 1][2])")
      thispoint = points[x2 + 1]
      thisslope = (thispoint[2] - startpoint[2]) // (thispoint[1] - startpoint[1])
      if thisslope ≤ smallest_slope # we do \leq even though the slope won't change so we can update finalpoint
        smallest_slope = thisslope
        finalx = x2
      end
      x2 = x2 + 1
    end # postcondition: x2 is the next vertex

    # add x2 to the vertices
    push!(vertices,points[finalx+1])

    x = finalx

  end

  vertices
end#function


#end#module
