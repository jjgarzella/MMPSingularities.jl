
function test_cases_k3surf_2()
  R, (w,x,y,z) = polynomial_ring(GF(2),4)

  f3_1 = w^3*x + w^3*y + w^3*z + w^2*x^2 + w^2*y^2 + w^2*z^2 + w*x^3 + w*x^2*y + w*x*y^2 + w*y^2*z + x^3*z + x*y^3 + x*y^2*z + x*z^3 + y^4 + y^3*z + y^2*z^2 + y*z^3

  cases = [f3_1]
  heights = [3]

  (cases, heights)
end

function test_cases_cy3_2()
  p = 2
  R, (u,w,x,y,z) = polynomial_ring(GF(p),["u","w","x","y","z"])
  
  f60 = x^5 + y^5 + z^5 + w^5 + u^5 + x*z^3*w + y*z*w^3 + x^2*z*u^2 + y^2*z^2*w + x*y^2*w*u + y*z*w*u^2

  cases = [f60]
  heights = [60]

  (cases, heights)
end

function test_cases_k3surf_3()
  R, (w,x,y,z) = polynomial_ring(GF(3),["w","x","y","z"])
 
  # (corrected) table 2 from arXiv:2204.10076
  f1 = x^4 + y^4 + z^4 + 2w^4 + x^2*y*w + y*z^2*w
  f2 = x^4 + 2y^4 + 2z^4 + 2w^4 + x*y*z^2
  f3 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2 + z^3*w
  f4 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2
  f5 = x^4 + y^4 + z^4 + w^4 + x^3*z + z^3*w + y*z^2*w + y*z*w^2
  f6 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x^2*y*z + x*z^3
  f7 = x^4 + y^4 + z^4 + w^4 + x*y^2*z + x*z^2*w + y*z*w^2 + y^2*z*w
  f8 = x^4 + x^2*y*z + x^2*y*w + 2x^2*z^2 + x*y*w^2 + 2*y^4 + y^3*w + z^4 + w^4
  f9 = x^4 + y^4 + z^4 + w^4 + x*y^3 + y^3*w + z^2*w^2 + 2*x*y*z^2 + y*z*w^2
  f10 = x^4 + 2*x^2*y*z + x^2*y*w + x*y^2*w + y^4 + y^3*w + y^2*z^2 + 2*y^2*z*w + y^2*w^2 + y*z^3 + y*z^2*w + y*z*w^2 + z^4 + z*w^3
  finfty = x^4 + y^4 + z^4 + w^4

  # this example came from a personal commutication with the authors of arXiv:2204.10076
  fm10 = 2w^4 + x*w^3 + w^2*x^2 + w^2*x*y + w^2*y*z + 2w*x^3 + 2w*x^2*y + w*x*y^2 + w*x*y*z + w*x*z^2 + w*y^3 + 2w*y^2*z + w*y*z^2 + x^4 + 2*x^3*y + x^2*y*z + x^2*z^2 + x*y^3 + x*y^2*z + 2x*y*z^2 + x*z^3 + 2*y^4 + 2*y^3*z + y^2*z^2 + y*z^3 + 2*z^4

  # I don't remember where the following example is from, but I think it
  # might also have been in a personal communication from the authors of arXiv:2204.10076
  (x1,x2,x3,x4) = (w,x,y,z)
  fr8 = x1^3*x2 + 2*x1^2*x2^2 + 2*x1^2*x2*x3 + x1^2*x3^2 + x1^2*x3*x4 + 2*x1*x2^2*x3 + 2*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + x2^4 + x2^3*x3 + 2*x2^3*x4 + x2^2*x3^2 + 2*x2*x3^3 + x2*x3^2*x4 + x2*x3*x4^2 + x2*x4^3 + x3^4 + x3^2*x4^2 + x3*x4^3

  cases = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, finfty, fm10, fr8]
  heights = [1,2,3,4,5,6,7,8,9,10,11,10,8]

  (cases, heights)
end

function test_cases_k3surf_5()
  p = 5
  R, (x1,x2,x3,x4) = polynomial_ring(GF(p),4)

  ftwo1 = 4*x1^4 + 2*x1^3*x2 + x1^3*x4 + 4*x1^2*x2^2 + 2*x1^2*x2*x3 + 2*x1^2*x3^2 + x1^2*x3*x4 + 3*x1*x2^3 + 4*x1*x2^2*x3 + 4*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 4*x2^4 + 2*x2^3*x3 + 4*x2^3*x4 + 4*x2^2*x3^2 + x2^2*x3*x4 + 2*x2^2*x4^2 + 3*x2*x3^3 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 2*x2*x4^3 + 2*x3^4 + 2*x3^3*x4 + 2*x3^2*x4^2 + x3*x4^3 + 4*x4^4

  ftwo2 = x1^4 + x1^3*x2 + 2*x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x3^2 + 4*x1^2*x3*x4 + x1^2*x4^2 + 4*x1*x2^3 + 4*x1*x2^2*x3 + 3*x1*x2^2*x4 + 2*x1*x2*x3^2 + 4*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 4*x1*x3^3 + 2*x1*x4^3 + x2^3*x3 + x2^3*x4 + 4*x2^2*x3^2 + 3*x2^2*x3*x4 + 4*x2*x3^3 + x2*x3^2*x4 + x2*x3*x4^2 + 2*x2*x4^3 + x3^2*x4^2 + x3*x4^3

  fthree1 = 2*x1^4 + x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x2*x4 + x1^2*x3^2 + 4*x1^2*x3*x4 + 3*x1^2*x4^2 + 4*x1*x2^3 + 3*x1*x2^2*x3 + x1*x2^2*x4 + 2*x1*x2*x3^2 + 3*x1*x2*x3*x4 + x1*x3^3 + 4*x1*x3*x4^2 + 2*x1*x4^3 + x2^3*x3 + 3*x2^3*x4 + 4*x2^2*x3^2 + 4*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^4 + 3*x3^3*x4 + 2*x3*x4^3 + 3*x4^4

  ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3

  ffive1 = 4*x1^2*x2^2 + 4*x1^2*x2*x3 + 2*x1^2*x2*x4 + 3*x1^2*x3*x4 + 4*x1*x2^2*x3 + 3*x1*x2*x3*x4 + x1*x3^3 + x1*x3^2*x4 + 3*x1*x3*x4^2 + x1*x4^3 + x2^4 + 2*x2*x3^2*x4 + 4*x2*x3*x4^2

  fsix1 = 4*x1^3*x2 + 4*x1^3*x4 + 4*x2^4 + 2*x2*x3^2*x4 + x2*x3*x4^2 + x3^3*x4

  fseven1 = 4*x1^3*x4 + x1*x2^3 + 2*x1*x3^3 + 2*x2^3*x4 + x2^2*x4^2 + x3^3*x4 + 3*x3^2*x4^2 + 2*x3*x4^3

  feight1 = 2*x1^3*x3 + x1^2*x2^2 + x1^2*x4^2 + x1*x2^2*x4 + 3*x1*x3^2*x4 + x2^4 + 2*x2^3*x3 + 3*x2^2*x4^2 + 4*x2*x3^2*x4 + 3*x3*x4^3

  fnine1 = 3*x1^4 + 3*x1^3*x2 + 3*x1^3*x3 + x1^2*x2^2 + 3*x1^2*x2*x3 + 3*x1^2*x2*x4 + 3*x1^2*x3^2 + 2*x1^2*x3*x4 + 2*x1^2*x4^2 + 4*x1*x2^3 + 2*x1*x2^2*x3 + 4*x1*x2*x3^2 + 2*x1*x2*x3*x4 + 4*x1*x2*x4^2 + x1*x3^3 + 3*x1*x3^2*x4 + 3*x1*x3*x4^2 + x1*x4^3 + 3*x2^3*x3 + 4*x2^3*x4 + 3*x2^2*x3*x4 + x2^2*x4^2 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 4*x2*x4^3 + 3*x3*x4^3 + 4*x4^4
  
  ften1 = 2*x1^4 + 4*x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + x1^2*x2^2 + 2*x1^2*x2*x3 + 2*x1^2*x2*x4 + 4*x1^2*x3^2 + 4*x1^2*x3*x4 + 2*x1^2*x4^2 + x1*x2^3 + 4*x1*x2^2*x4 + 3*x1*x2*x3^2 + 3*x1*x2*x4^2 + 2*x1*x3^3 + 3*x1*x3^2*x4 + 2*x1*x3*x4^2 + x1*x4^3 + 3*x2^4 + 2*x2^3*x3 + 2*x2^3*x4 + 4*x2^2*x3^2 + 3*x2^2*x3*x4 + 3*x2^2*x4^2 + x2*x3^3 + 2*x2*x3*x4^2 + 2*x2*x4^3 + 4*x3^4 + x3^3*x4 + 3*x3^2*x4^2 + 4*x3*x4^3 + 3*x4^4

  finfty = x1^4 + x2^4 + x3^4 + x4^4 + x1*x2*x3*x4

  cases = [ftwo1,ftwo2,fthree1,ffour1,ffive1,fsix1,fseven1,feight1,fnine1,ften1,finfty]
  heights = [2,2,3,4,5,6,7,8,9,10,11]
  (cases, heights)
end

function test_cases_k3surf_7()
  ([],[])   
end

function test_cases_k3surf_11()
  ([],[])   
end

function test_cases_k3surf_13()

end
