module Experiments

include("QuasiFSplitness.jl")
include("RandomPolynomials.jl")

using .QuasiFSplitness
using .RandomPolynomials

using Oscar

function test_Δ₁(p,polys)
  
  Δ₁.(p,polys) .== Δ₁l(p,polys)

end#function

"""
If outfilename is not the empty string,
outputs a file containing the data of arr.

Outerwise, does nothing.

"""
function output_array_to_file(arr,outfilename)

  if outfilename != ""
    open(outfilename, "w") do outfile
      for r in arr 
        println(outfile,r)
      end
    end
  end

end#function

function find_highest_quasi_F_split_height_random(p,n,d,cutoff,howhigh,numtries;outfilename="")
  R, vars = polynomial_ring(GF(p),n)

  results = []
  for i in 1:numtries
    f = RandomPolynomials.random_homog_poly_mod(p,vars,d)
    h = QuasiFSplitness.quasiFSplitHeight_CY_lift(p,f,cutoff)
    if howhigh ≤ h 
      println("Found variety with quasi-F-split height $h: $f")
      push!(results,(h,f))
    end
  end
  
  output_array_to_file(results,outfilename)

  results
end#function


function quasi_F_split_height_random_bargraph(p,n,d,cutoff,howhigh,numtries;outfilename="")
  R, vars = polynomial_ring(GF(p),n)

  results = []
  bars = zeros(Int,cutoff+2)
  for i in 1:numtries
    f = RandomPolynomials.random_homog_poly_mod(p,vars,d)
    h = QuasiFSplitness.quasiFSplitHeight_CY_lift(p,f,cutoff)
    bars[h] = bars[h] + 1
    if howhigh ≤ h 
      println("Found variety with quasi-F-split height $h: $f")
      push!(results,(h,f))
    end
  end

  output_array_to_file(results,outfilename)

  (results,bars)
end#function

function perturb_quasi_f_split_height(p,d,startf,cutoff,numtries;saveheight=-1,outfilename="")
  vars = gens(parent(startf))
  mons = RandomPolynomials.random_monomials(vars,d,numtries)
  coefs = rand(1:p-1,numtries)
  terms = coefs .* mons

  results = []
  bars = zeros(Int,cutoff+2)
  for term in terms
    fm = startf + term
    h = QuasiFSplitness.quasiFSplitHeight_CY_lift(p,fm,cutoff)
    bars[h] = bars[h] + 1
    if saveheight != -1 && saveheight ≤ h
      println("Found variety wtih quasi-F-split height $h: $fm")
      push!(results,(h,fm))
    end
  end

  output_array_to_file(results,outfilename)

  (results,bars)
end#function

function fermat_hypersurface_heights(mind,maxd,minn,maxn,minp,maxp,cutoff=100)
  result = Dict()

  for n = minn:maxn
    for p = minp:maxp

      R, vars = polynomial_ring(GF(p),n)

      for d = mind:maxd

        fermat = sum(vars .^ d)

        h = QuasiFSplitness.quasiFSplitHeight(p,fermat,cutoff)

        println("Height of $fermat is $h")

        result[(n,d,p)] = h
      end
    end
  end

end#function

function height_cubicfourfold(newtonpoly)
  firstbreak = newtonpoly[2]
  firstslope = firstbreak[2] / firstbreak[1]
  for h in 1:10
    if firstslope ≈ 2 - (1/h)
      return h
    end
  end
  if firstslope ≈ 2
    return 11 # really infinity!!
  end
  println("This newton polygon doesn't seem to be a cubic fourfold")
  return -1
end#function


function testquasi1Fsplit_guess(orbdata,newtpolydict)
  counter = 0
  for i in eachindex(orbdata)
    if haskey(newtpolydict,i)
      # compute the quasi 1 f split height
      q1fsH = QuasiFSplitness.quasi1FSplitHeight_2CY_lift(2,orbdata[i],10)
      # compute the newton polygon height
      newtH = height_cubicfourfold(newtpolydict[i])

      # take care for the case of infinity
      if 11 < q1fsH 
        q1fsH = 11
      end

      # compare them, if one doesn't match, then print it
      if q1fsH != newtH
        println("Fourfold at index $i is a counterexample, qfs ht = $q1fsH, ht = $newtH")
      end
      counter += 1
    end

    if counter % 100 == 0
      println("$counter surfaces checked")
    end
  end
end#function

end#module
