
"""
Evaluates the functino qfs_fn on
cases[i] for all i and compares the height to 
heights[i]

qfs_fn is expected to have the function signature of
quasiFSplitHeight, that is:

qfs_fn(p,f,cutoff) returns the quasi-F-split height
h if h ≤ cutoff, otherwise it returns cutoff + 1 if the algorithm
does not know if the height is infinity and cutoff + 2 if the
algorithm does know that the height is infinity.

"""
function evaluate_test_set(qfs_fn, ps, cases, cutoffs, heights)
  for i in 1:length(cases)
    println("Testing example $i")
    @test qfs_fn(ps[i],cases[i],cutoffs[i]) == heights[i]
  end
end

function test_heights_cy3(p,qfs_fn)
  cases, heights = test_cases_cy3(p)
  l = length(cases)
  ps = fill(p,l)
  cutoffs = fill(101,l)
  evaluate_test_set(qfs_fn,ps,cases,cutoffs,heights)
end

function test_heights_k3surf(p,qfs_fn)
  cases, heights = test_cases_k3surf(p)
  l = length(cases)
  ps = fill(p,l)
  cutoffs = fill(10,l)
  evaluate_test_set(qfs_fn,ps,cases,cutoffs,heights)
end

"""
Note: calling this method is expensive
on the first call, as it will cause Julia
to compile all of the methods which get returned
in the array.
"""
function algorithms_k3surf(p)
  
  # create pregen manually, this will speed things up a lot
  pregen_k3 = MMPSingularities.pregen_delta1(4,p)

  stateoftheart = [MMPSingularities.quasiFSplitHeight_CY_lift,
                   MMPSingularities.quasiFSplitHeight_CY_lift_sort,
                   (p,f,cutoff) -> MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p,f,cutoff,pregen_k3),
                   (p,f,cutoff) -> MMPSingularities.quasiFSplitHeight_CY_gpu(p,f,cutoff,pregen_k3)]

  others = [MMPSingularities.quasiFSplitHeight_CY_naive_expansion,
            MMPSingularities.quasiFSplitHeight_CY_lift_lazy,
            MMPSingularities.quasiFSplitHeight_CY_lift_matrix,
            MMPSingularities.quasiFSplitHeight_CY_lift_matrix_combined,
            MMPSingularities.quasiFSplitHeight_CY_formula]

  algs = [stateoftheart; others]


  if p == 2
    algs
  elseif p == 3
    algs
  elseif p == 5
    algs[1:4]
  elseif p == 7
    algs[2:4]
  elseif p == 11
    [algs[3]]
  elseif p == 13
    [algs[3]]
  else
    []
  end
    
end

function algorithms_cy3(p)
  # create pregen manually, this will speed things up a lot
  pregen_k3 = MMPSingularities.pregen_delta1(5,p)

  stateoftheart = [MMPSingularities.quasiFSplitHeight_CY_lift,
                   MMPSingularities.quasiFSplitHeight_CY_lift_sort,
                   (p,f,cutoff) -> MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu(p,f,cutoff,pregen_k3),
                   (p,f,cutoff) -> MMPSingularities.quasiFSplitHeight_CY_gpu(p,f,cutoff,pregen_k3)]

  if p == 2
    stateoftheart
  else
    []
  end
end


"""
Manually select the algorithms that make sense to test for various p so that
the tests take a reasonable amount of time.
"""
function test_heights_all_performant()


  ps_k3 = [2,3,5,7,11,13]
  ps_cy3 = 2

  # K3 surfaces
  for p in ps_k3
    for alg in algorithms_k3(p)
      test_heights_k3surf(p,alg)
    end
  end

  # CY threefolds
  for p in ps_cy3
    for alg in algorithms_cy3(p)
      test_heights_cy3(p,alg)
    end
  end
end

"""
Runs the test case `case` and times it.

useful for regression testing, though actually setting
up regression testing seems like overkill right now.
"""
function test_time_qfs(p,case,n,qfs_fn)
  qfs_height_fn = qfs_fn#MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu

  pregen = MMPSingularities.pregen_delta1(n,p)


#  R, (x1,x2,x3,x4) = polynomial_ring(GF(p),n)

 # ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3

  # prime the jitter

  if n == 4
    cutoff = 10
  elseif n == 5 
    cutoff = 101
  else
    cutoff = -1
    error("unsupported n: $n")
  end
      
  println("PRIMING THE JITTER:")
  CUDA.@time qfs_height_fn(p,case,cutoff,pregen)

  println()
  CUDA.@time qfs_height_fn(p,case,cutoff,pregen)

  println()
  #Profile.Allocs.clear()
  #Profile.Allocs.@profile qfs_height_fn(p,case,cutoff,pregen)

  #Profile.clear()
  #Profile.@profile qfs_height_fn(p,case,cutoff,pregen)
end

function test_heights_all()
  #algs_cy3_2 = algorithms_cy3(2)

  #test_heights_cy3(2,algs_cy3[1]) # naive lift method

  #algs_k3_2 = algorithms_k3surf(2)
  #algs_k3_3 = algorithms_k3surf(3)
  #algs_k3_5 = algorithms_k3surf(5)

  #test_heights_k3surf(3,algs_k3_3[7]) # matrix
  #test_heights_k3surf(2,algs_k3_2[3]) # sort_gpu
  #test_heights_k3surf(3,algs_k3_3[3]) # sort_gpu
  #test_heights_k3surf(5,algs_k3_5[3]) # sort_gpu

  cases_2, _ = test_cases_k3surf(2)
  ff = cases_2[9]
  
  cases_2_cy, _ = test_cases_cy3(2)
  ff3 = cases_2_cy[2]
  #test_matrix_K3(2,ff,MMPSingularities.matrix_of_multiply_then_split_sortmodp_kronecker)
  test_time_qfs(2,ff3,5,
                (p,case,cutoff,pregen) -> 
                MMPSingularities.quasiFSplitHeight_CY_lift_sort(p,case,cutoff))

  #test_heights_k3surf(5,algs_k3[4]) # gpu (gpu delta1, no matrix)

end

#function test_high_height_CY3_char_2()
#  p = 2
#  best_qfs_height = f -> MMPSingularities.quasiFSplitHeight_CY_lift(2,f,101)
#  
#  R, (u,w,x,y,z) = polynomial_ring(GF(p),["u","w","x","y","z"])
#  
#  f60 = x^5 + y^5 + z^5 + w^5 + u^5 + x*z^3*w + y*z*w^3 + x^2*z*u^2 + y^2*z^2*w + x*y^2*w*u + y*z*w*u^2
#
#  @test best_qfs_height(f60) == 60
#end
#
#function test_K3_3()
#  p = 3
#  best_qfs_height = f -> MMPSingularities.quasiFSplitHeight_CY_lift_matrix(3,f,10)
#  #best_qfs_height = f -> MMPSingularities.quasiFSplitHeight_CY_lift(3,f,10)
#
#  R, (w,x,y,z) = polynomial_ring(GF(3),["w","x","y","z"])
# 
#  f1 = x^4 + y^4 + z^4 + 2w^4 + x^2*y*w + y*z^2*w
## 2*w^4 + w*x^2*y + w*y*z^2 + x^4 + y^4 + z^4
#  f2 = x^4 + 2y^4 + 2z^4 + 2w^4 + x*y*z^2
## 2*w^4 + x^4 + x*y*z^2 + 2*y^4 + 2*z^4
#  f3 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2 + z^3*w
## w^4 + w*z^3 + x^4 + x^2*z^2 + x*y*z^2 + y^4 + z^4
#  f4 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2
## w^4 + x^4 + x^2*z^2 + x*y*z^2 + y^4 + z^4
#  f5 = x^4 + y^4 + z^4 + w^4 + x^3*z + z^3*w + y*z^2*w + y*z*w^2
## w^4 + w^2*y*z + w*y*z^2 + w*z^3 + x^4 + x^3*z + y^4 + z^4
#  f6 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x^2*y*z + x*z^3
## w^4 + x^4 + x^2*y*z + x^2*z^2 + y^4 + z^4
#  f7 = x^4 + y^4 + z^4 + w^4 + x*y^2*z + x*z^2*w + y*z*w^2 + y^2*z*w
## w^4 + w^2*y*z + w*x*z^2 + w*y*z^2 + x^4 + x*y^2*z + y^4 + z^4
#  f8 = x^4 + x^2*y*z + x^2*y*w + 2x^2*z^2 + x*y*w^2 + 2*y^4 + y^3*w + z^4 + w^4
## w^4 + w^2*x*y + w*x^2*y + w*y^3 + x^4 + x^2*y*z + 2*x^2*z^2 + 2*y^4 + z^4
#  f9 = x^4 + y^4 + z^4 + w^4 + x*y^3 + y^3*w + z^2*w^2 + 2*x*y*z^2 + y*z*w^2
## w^4 + w^2*y*z + w^2*z^2 + w*y^3 + x^4 + x*y^3 + 2*x*y*z^2 + y^4 + z^4
#  f10 = x^4 + 2*x^2*y*z + x^2*y*w + x*y^2*w + y^4 + y^3*w + y^2*z^2 + 2*y^2*z*w + y^2*w^2 + y*z^3 + y*z^2*w + y*z*w^2 + z^4 + z*w^3
## w^3*z + w^2*y^2 + w^2*y*z + w*x^2*y + w*x*y^2 + w*y^3 + 2*w*y^2*z + w*y*z^2 + x^4 + 2*x^2*y*z + y^4 + y^2*z^2 + y*z^3 + z^4
#  finfty = x^4 + y^4 + z^4 + w^4
## w^4 + x^4 + y^4 + z^4
#  fm10 = 2w^4 + x*w^3 + w^2*x^2 + w^2*x*y + w^2*y*z + 2w*x^3 + 2w*x^2*y + w*x*y^2 + w*x*y*z + w*x*z^2 + w*y^3 + 2w*y^2*z + w*y*z^2 + x^4 + 2*x^3*y + x^2*y*z + x^2*z^2 + x*y^3 + x*y^2*z + 2x*y*z^2 + x*z^3 + 2*y^4 + 2*y^3*z + y^2*z^2 + y*z^3 + 2*z^4
#
#  (x1,x2,x3,x4) = (w,x,y,z)
#  fr8 = x1^3*x2 + 2*x1^2*x2^2 + 2*x1^2*x2*x3 + x1^2*x3^2 + x1^2*x3*x4 + 2*x1*x2^2*x3 + 2*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + x2^4 + x2^3*x3 + 2*x2^3*x4 + x2^2*x3^2 + 2*x2*x3^3 + x2*x3^2*x4 + x2*x3*x4^2 + x2*x4^3 + x3^4 + x3^2*x4^2 + x3*x4^3
#
#
#  @test best_qfs_height(f1) == 1
#  @test best_qfs_height(f2) == 2
#  @test best_qfs_height(f3) == 3
#  @test best_qfs_height(f4) == 4
#  @test best_qfs_height(f5) == 5
#  @test best_qfs_height(f6) == 6
#  @test best_qfs_height(f7) == 7
#  @test best_qfs_height(f8) == 8
#  @test best_qfs_height(f9) == 9
#  @test best_qfs_height(f10) == 10
#  infty_substitute = best_qfs_height(finfty)
#  @test infty_substitute == 11 || infty_substitute == 12
#  @test best_qfs_height(fm10) == 10 
#  @test best_qfs_height(fr8) == 8 
#end
#
#function test_K3_5()
#
#  p = 5
#  pregen = MMPSingularities.pregen_delta1(4, p)
#  qfs_height_fn(f) = MMPSingularities.quasiFSplitHeight_CY_gpu(p, f, 10, pregen)
#
#  R, (x1,x2,x3,x4) = polynomial_ring(GF(p),4)
#
#  ftwo1 = 4*x1^4 + 2*x1^3*x2 + x1^3*x4 + 4*x1^2*x2^2 + 2*x1^2*x2*x3 + 2*x1^2*x3^2 + x1^2*x3*x4 + 3*x1*x2^3 + 4*x1*x2^2*x3 + 4*x1*x2^2*x4 + 2*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 3*x1*x3^3 + x1*x3^2*x4 + x1*x3*x4^2 + x1*x4^3 + 4*x2^4 + 2*x2^3*x3 + 4*x2^3*x4 + 4*x2^2*x3^2 + x2^2*x3*x4 + 2*x2^2*x4^2 + 3*x2*x3^3 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 2*x2*x4^3 + 2*x3^4 + 2*x3^3*x4 + 2*x3^2*x4^2 + x3*x4^3 + 4*x4^4
#
#  ftwo2 = x1^4 + x1^3*x2 + 2*x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x3^2 + 4*x1^2*x3*x4 + x1^2*x4^2 + 4*x1*x2^3 + 4*x1*x2^2*x3 + 3*x1*x2^2*x4 + 2*x1*x2*x3^2 + 4*x1*x2*x3*x4 + 3*x1*x2*x4^2 + 4*x1*x3^3 + 2*x1*x4^3 + x2^3*x3 + x2^3*x4 + 4*x2^2*x3^2 + 3*x2^2*x3*x4 + 4*x2*x3^3 + x2*x3^2*x4 + x2*x3*x4^2 + 2*x2*x4^3 + x3^2*x4^2 + x3*x4^3
#
#  fthree1 = 2*x1^4 + x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + x1^2*x2*x3 + 4*x1^2*x2*x4 + x1^2*x3^2 + 4*x1^2*x3*x4 + 3*x1^2*x4^2 + 4*x1*x2^3 + 3*x1*x2^2*x3 + x1*x2^2*x4 + 2*x1*x2*x3^2 + 3*x1*x2*x3*x4 + x1*x3^3 + 4*x1*x3*x4^2 + 2*x1*x4^3 + x2^3*x3 + 3*x2^3*x4 + 4*x2^2*x3^2 + 4*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + 4*x2*x3*x4^2 + 3*x2*x4^3 + 4*x3^4 + 3*x3^3*x4 + 2*x3*x4^3 + 3*x4^4
#
#  ffour1 = 4*x1^4 + 2*x1^3*x3 + 4*x1^3*x4 + 3*x1^2*x2^2 + 3*x1^2*x2*x3 + 4*x1^2*x2*x4 + 2*x1^2*x3*x4 + x1^2*x4^2 + 3*x1*x2^3 + x1*x2^2*x3 + x1*x2^2*x4 + x1*x2*x3^2 + x1*x2*x3*x4 + x1*x2*x4^2 + 2*x1*x3^3 + 2*x1*x3^2*x4 + x1*x3*x4^2 + 2*x1*x4^3 + 4*x2^4 + 3*x2^3*x3 + x2^3*x4 + 3*x2^2*x3^2 + 3*x2^2*x3*x4 + x2^2*x4^2 + 2*x2*x3^3 + 3*x2*x3^2*x4 + x2*x3*x4^2 + 3*x2*x4^3 + 3*x3^4 + 2*x3^3*x4 + 4*x3^2*x4^2 + x3*x4^3
#
#  ffive1 = 4*x1^2*x2^2 + 4*x1^2*x2*x3 + 2*x1^2*x2*x4 + 3*x1^2*x3*x4 + 4*x1*x2^2*x3 + 3*x1*x2*x3*x4 + x1*x3^3 + x1*x3^2*x4 + 3*x1*x3*x4^2 + x1*x4^3 + x2^4 + 2*x2*x3^2*x4 + 4*x2*x3*x4^2
#
#  fsix1 = 4*x1^3*x2 + 4*x1^3*x4 + 4*x2^4 + 2*x2*x3^2*x4 + x2*x3*x4^2 + x3^3*x4
#
#  fseven1 = 4*x1^3*x4 + x1*x2^3 + 2*x1*x3^3 + 2*x2^3*x4 + x2^2*x4^2 + x3^3*x4 + 3*x3^2*x4^2 + 2*x3*x4^3
#
#  feight1 = 2*x1^3*x3 + x1^2*x2^2 + x1^2*x4^2 + x1*x2^2*x4 + 3*x1*x3^2*x4 + x2^4 + 2*x2^3*x3 + 3*x2^2*x4^2 + 4*x2*x3^2*x4 + 3*x3*x4^3
#
#  fnine1 = 3*x1^4 + 3*x1^3*x2 + 3*x1^3*x3 + x1^2*x2^2 + 3*x1^2*x2*x3 + 3*x1^2*x2*x4 + 3*x1^2*x3^2 + 2*x1^2*x3*x4 + 2*x1^2*x4^2 + 4*x1*x2^3 + 2*x1*x2^2*x3 + 4*x1*x2*x3^2 + 2*x1*x2*x3*x4 + 4*x1*x2*x4^2 + x1*x3^3 + 3*x1*x3^2*x4 + 3*x1*x3*x4^2 + x1*x4^3 + 3*x2^3*x3 + 4*x2^3*x4 + 3*x2^2*x3*x4 + x2^2*x4^2 + 4*x2*x3^2*x4 + 4*x2*x3*x4^2 + 4*x2*x4^3 + 3*x3*x4^3 + 4*x4^4
#  
#  ften1 = 2*x1^4 + 4*x1^3*x2 + 3*x1^3*x3 + x1^3*x4 + x1^2*x2^2 + 2*x1^2*x2*x3 + 2*x1^2*x2*x4 + 4*x1^2*x3^2 + 4*x1^2*x3*x4 + 2*x1^2*x4^2 + x1*x2^3 + 4*x1*x2^2*x4 + 3*x1*x2*x3^2 + 3*x1*x2*x4^2 + 2*x1*x3^3 + 3*x1*x3^2*x4 + 2*x1*x3*x4^2 + x1*x4^3 + 3*x2^4 + 2*x2^3*x3 + 2*x2^3*x4 + 4*x2^2*x3^2 + 3*x2^2*x3*x4 + 3*x2^2*x4^2 + x2*x3^3 + 2*x2*x3*x4^2 + 2*x2*x4^3 + 4*x3^4 + x3^3*x4 + 3*x3^2*x4^2 + 4*x3*x4^3 + 3*x4^4
#
#  @test qfs_height_fn(ftwo1) == 2
#  @test qfs_height_fn(ftwo2) == 2
#  @test qfs_height_fn(fthree1) == 3
#  @test qfs_height_fn(ffour1) == 4
#  @test qfs_height_fn(ffive1) == 5
#  @test qfs_height_fn(fsix1) == 6
#  @test qfs_height_fn(fseven1) == 7
#  @test qfs_height_fn(feight1) == 8
#  @test qfs_height_fn(fnine1) == 9
#  @test qfs_height_fn(ften1) == 10
#
#end



#function test_K3_3_gpu()
#    qfs_height_fn = MMPSingularities.quasiFSplitHeight_CY_lift_sort_gpu
#
#    pregen = MMPSingularities.pregen_delta1(4,3)
#    R, (x, y, z, w) = polynomial_ring(GF(3),4)
#
#    f1 = x^4 + y^4 + z^4 + 2w^4 + x^2*y*w + y*z^2*w
#    #2*w^4 + w*x^2*y + w*y*z^2 + x^4 + y^4 + z^4
#    f2 = x^4 + 2y^4 + 2z^4 + 2w^4 + x*y*z^2
#    #2*w^4 + x^4 + x*y*z^2 + 2*y^4 + 2*z^4
#    f3 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2 + z^3*w
#    #w^4 + w*z^3 + x^4 + x^2*z^2 + x*y*z^2 + y^4 + z^4
#    f4 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x*y*z^2
#    #w^4 + x^4 + x^2*z^2 + x*y*z^2 + y^4 + z^4
#    f5 = x^4 + y^4 + z^4 + w^4 + x^3*z + z^3*w + y*z^2*w + y*z*w^2
#    #w^4 + w^2*y*z + w*y*z^2 + w*z^3 + x^4 + x^3*z + y^4 + z^4
#    f6 = x^4 + y^4 + z^4 + w^4 + x^2*z^2 + x^2*y*z + x*z^3
#    #w^4 + x^4 + x^2*y*z + x^2*z^2 + y^4 + z^4
#    f7 = x^4 + y^4 + z^4 + w^4 + x*y^2*z + x*z^2*w + y*z*w^2 + y^2*z*w
#    #w^4 + w^2*y*z + w*x*z^2 + w*y*z^2 + x^4 + x*y^2*z + y^4 + z^4
#    f8 = x^4 + x^2*y*z + x^2*y*w + 2x^2*z^2 + x*y*w^2 + 2*y^4 + y^3*w + z^4 + w^4
#    #w^4 + w^2*x*y + w*x^2*y + w*y^3 + x^4 + x^2*y*z + 2*x^2*z^2 + 2*y^4 + z^4
#    f9 = x^4 + y^4 + z^4 + w^4 + x*y^3 + y^3*w + z^2*w^2 + 2*x*y*z^2 + y*z*w^2
#    #w^4 + w^2*y*z + w^2*z^2 + w*y^3 + x^4 + x*y^3 + 2*x*y*z^2 + y^4 + z^4
#    f10 = x^4 + 2*x^2*y*z + x^2*y*w + x*y^2*w + y^4 + y^3*w + y^2*z^2 + 2*y^2*z*w + y^2*w^2 + y*z^3 + y*z^2*w + y*z*w^2 + z^4 + z*w^3
#    #w^3*z + w^2*y^2 + w^2*y*z + w*x^2*y + w*x*y^2 + w*y^3 + 2*w*y^2*z + w*y*z^2 + x^4 + 2*x^2*y*z + y^4 + y^2*z^2 + y*z^3 + z^4
#    finfty = x^4 + y^4 + z^4 + w^4
#    #w^4 + x^4 + y^4 + z^4
#    fm10 = 2w^4 + x*w^3 + w^2*x^2 + w^2*x*y + w^2*y*z + 2w*x^3 + 2w*x^2*y + w*x*y^2 + w*x*y*z + w*x*z^2 + w*y^3 + 2w*y^2*z + w*y*z^2 + x^4 + 2*x^3*y + x^2*y*z + x^2*z^2 + x*y^3 + x*y^2*z + 2x*y*z^2 + x*z^3 + 2*y^4 + 2*y^3*z + y^2*z^2 + y*z^3 + 2*z^4
#
#    @test qfs_height_fn(3, f1, 10, pregen) == 1
#    @test qfs_height_fn(3, f2, 10, pregen) == 2
#    @test qfs_height_fn(3, f3, 10, pregen) == 3
#    @test qfs_height_fn(3, f4, 10, pregen) == 4
#    @test qfs_height_fn(3, f5, 10, pregen) == 5
#    @test qfs_height_fn(3, f6, 10, pregen) == 6
#    @test qfs_height_fn(3, f7, 10, pregen) == 7
#    @test qfs_height_fn(3, f8, 10, pregen) == 8
#    @test qfs_height_fn(3, f9, 10, pregen) == 9
#    @test qfs_height_fn(3, f10, 10, pregen) == 10
#    infty_substitute = qfs_height_fn(3, finfty, 10, pregen)
#    @test infty_substitute == 11 || infty_substitute == 12
#end
