module QuasiFSplitness

using Oscar

exponent_vectors(poly) = leading_exponent_vector.(terms(poly))

"""
Calculates if the hypersuface defined by the 
polynomial poly is F-split

note that p must be prime for this to have mathematical meaning
"""
function isFSplit(p,poly)
  #maybe TODO: check that p is prime


  fpower = poly^(p-1)

  for exponent_vector in exponent_vectors(fpower)
    if all(exponent_vector .< p)
      # We are F-split!!
      return true
    end
  end

  false
end#function

end#module
