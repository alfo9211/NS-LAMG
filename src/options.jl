include("../src/smoothers.jl"); 


type options

  ##########################################
  # Options for solve Ax = 0
  ##########################################
  # rsv_smoother = smoother used for finding the right singular vector
  rsv_smoother::Function 
  #w = weight for smoother (if weighted Jacobi)
  rsv_w::Float64 
  # V(nu1,nu2)-cycle for right signular vector
  rsv_nu1::Int64 
  rsv_nu2::Int64 
  # rsv_theta = Strength of Connection matrix parameter for rsv
  rsv_theta::Float64
  # ncoarse = size of matrix for direct solve
  rsv_ncoarse::Int64 
  # rsv_maxLev = max number of levels for finding rsv
  rsv_maxLev::Int64
  # rsv_max_num_vcycles = number of V-Cycles to do for finding rsv
  rsv_max_num_vcycle::Int64 
  # rsv_tol = error tolerance for finding rsv
  rsv_tol::Float64 
  # rsv_conv_cutoff = cutoff for multigrid cycle convergence factor ( should be <= 1 ) for finding rsv
  rsv_conv_cutoff::Float64



  ##########################################
  # Options for solve Ax = b
  ##########################################

	# smoother = vCycle smoother
  smoother::Function 
  # w = weight for smoother (if weighted Jacobi)
  w::Float64 
  # V(nu1,nu2)-cycle
  nu1::Int64 
  nu2::Int64 
  # strength of connection threshold
  theta::Float64
  # ncoarse = size of matrix for direct solve
  ncoarse::Int64 
  #max number of levels
  maxLev::Int64
  # max_num_vcycles = number of V-Cycles to do
  max_num_vcycle::Int64 
  # tol = error tolerance
  tol::Float64 
  # conv_cutoff = cutoff for multigrid cycle convergence factor ( should be <= 1 )
  conv_cutoff::Float64



  ##########################################
  # Options for both solve Ax = 0 and Ax = b
  ##########################################
  #normtype = norm used
  normtype::Int64
   # gamma = Cycle Index
  gamma::Float64
  # g = coarsening work guard, or "guard"
  g::Float64
   # I_max = max aggregation sweeps
  I_max::Int64 
  # d_0 = number of test vectors to use on finest level
  d_0::Int64 
  # d_add = number of additional test vectors to use on each coarser aggregation level
  d_add::Int64
  # d_max = maximum number of test vectors to use
  d_max::Int64
   # nu = number of GS relaxation sweeps on test vectors
  nu::Int64 

  nprev::Int64 




  ##########################################
  # Options Low-degree Elimination
  ##########################################
  minEliminationFraction::Float64 


  options(rsv_smoother,rsv_w, rsv_nu1, rsv_nu2, rsv_theta, rsv_ncoarse, rsv_maxLev, rsv_max_num_vcycle, rsv_tol, rsv_conv_cutoff, smoother, w, nu1, nu2, theta,ncoarse,maxLev, max_num_vcycle,tol,conv_cutoff, normtype,gamma,g,I_max,d_0,d_add,d_max,nu,nprev,minEliminationFraction) = new(rsv_smoother,copy(rsv_w),copy(rsv_nu1), copy(rsv_nu2), copy(rsv_theta), copy(rsv_ncoarse), copy(rsv_maxLev), copy(rsv_max_num_vcycle),copy(rsv_tol), copy(rsv_conv_cutoff), smoother, copy(w), copy(nu1), copy(nu2), copy(theta), copy(ncoarse), copy(maxLev), copy(max_num_vcycle), copy(tol), copy(conv_cutoff), copy(normtype),copy(gamma),copy(g),copy(I_max),copy(d_0),copy(d_add),copy(d_max),copy(nu),copy(nprev),copy(minEliminationFraction))
end

type solveStats
  # rsv_conFacHistory = convergence factor History
  rsv_conFacHistory::AbstractVector
  # Total time for finding rsv
  rsvT::Float64
  #number of vcycles for rsv 
  rsv_numCycles::Int64
  # conFacHistory = convergence factor History
  conFacHistory::AbstractVector
  # setupT = setup time
  setupT::Float64
  #solveT = solve time
  solveT::Float64
converged::Bool
  solveStats(rsv_conFacHistory,rsvT,rsv_numCycles,conFacHistory,setupT,solveT,converged) = new(copy(rsv_conFacHistory),copy(rsvT),copy(rsv_numCycles),copy(conFacHistory),copy(setupT),copy(solveT), copy(converged))
end

function setSolveStats()
  rsv_conFacHistory = Array(Float64,0)
  rsvT = 0.0
  rsv_numCycles = 0
  conFacHistory = Array(Float64,0)
  setupT = 0.0
  solveT = 0.0
  converged = false
  return solveStats(rsv_conFacHistory,rsvT,rsv_numCycles,conFacHistory,setupT,solveT,converged)
end

function setopts()

  ##########################################
  # Options for solve Ax = 0
  ##########################################
  # rsv_smoother = smoother used for finding the right singular vector
  rsv_smoother = wjacobi!
  #w = weight for smoother (if weighted Jacobi)
  rsv_w=2/3
  # V(nu1,nu2)-cycle for right signular vector
  rsv_nu1 = 2
  rsv_nu2 =2
  # rsv_theta = Strength of Connection matrix parameter for rsv
  rsv_theta= .8
  # ncoarse = size of matrix for direct solve
  rsv_ncoarse = 100 
  # rsv_maxLev = max number of levels for finding rsv
  rsv_maxLev = 100
  # rsv_max_num_vcycles = number of V-Cycles to do for finding rsv
  rsv_max_num_vcycle = 100
  # rsv_tol = error tolerance for finding rsv
  rsv_tol = 10.0^(-7) 
  # rsv_conv_cutoff = cutoff for multigrid cycle convergence factor ( should be <= 1 ) for finding rsv
  rsv_conv_cutoff = 1



  ##########################################
  # Options for solve Ax = b
  ##########################################

  # smoother = vCycle smoother
  smoother= wjacobib!
  # smoother = gmres_smoother!
  # smoother = richardson!
  # smoother= gaussseidel!
  # w = weight for smoother (if weighted Jacobi)
  w = 2/3 
  # w =1/20
  # V(nu1,nu2)-cycle
  nu1 = 2
  nu2 =2
  # strength of connection threshold
  theta = .8
  # ncoarse = size of matrix for direct solve
  ncoarse = 100
  #max number of levels
  maxLev = 100
  # max_num_vcycles = number of V-Cycles to do
  max_num_vcycle = 100
  # tol = error tolerance
  tol = 10.0^(-7) 
  # conv_cutoff = cutoff for multigrid cycle convergence factor ( should be <= 1 )
  conv_cutoff = 1



  ##########################################
  # Options for both solve Ax = 0 and Ax = b
  ##########################################
  #normtype = norm used
  normtype = 1

  gamma = 1.0
  # gamma = 0.1
  # g = coarsening work guard, or "guard"
  # g = 0.7
  g = 0.5
  # I_max = max aggregation sweeps
  I_max = 10
  # d_0 = number of test vectors to use on finest level
  d_0 = 1
  # d_add = number of additional test vectors to use on each coarser aggregation level
  d_add = 1
  # d_max = maximum number of test vectors to use
  d_max = 10
  # nu = number of GS relaxation sweeps on test vectors
  nu = 10

nprev = 3

 ##########################################
  # Options Low-degree Elimination
  ##########################################
   minEliminationFraction= .01
	
	return options(rsv_smoother, rsv_w,rsv_nu1, rsv_nu2, rsv_theta, rsv_ncoarse, rsv_maxLev, rsv_max_num_vcycle, rsv_tol,rsv_conv_cutoff, smoother, w,nu1, nu2, theta, ncoarse,maxLev,max_num_vcycle,tol,conv_cutoff,normtype,gamma,g,I_max,d_0,d_add,d_max,nu,nprev,minEliminationFraction)
end
