

type complexHistory
  # store number of MatVecs
  numA::Int64
  numP::Int64
  numR::Int64
  # store number of nonzeros in operators
  nnzA::Int64
  nnzP::Int64
  nnzR::Int64

  complexHistory(numA,numP,numR,nnzA,nnzP,nnzR) = new(copy(numA),copy(numP),copy(numR),copy(nnzA),copy(nnzP),copy(nnzR))
end
# --------------------------------------------------------------
# Function that does the final solve for Ax=b for nonsymmetric graphs
# The hierarchy has been found and is stored in lv, stats will store 
# the solve statistics and opt hold all the options for the method. 
# -------------------------------------------------------------- 

function solveV!(lv::Array{level,1}, stats::solveStats, opts::options)
  println(" lvfun = ", lv[1].lvFunc)

 # if top level is elimination, restrict down, and do V-cycles below top level
  # if string(lv[2].lvFunc) == "elim" 
  #   println(size(lv[1].R)) 
    
  #   topV = 2
  #   # for elimination, restrict RHS
  #   lv[2].b = lv[1].R*lv[1].b
  #   # for elimination, inject x on coarse points
  #   lv[2].x = lv[1].Rinj*lv[1].x
  #   # lv[2].opcount.R+=1

  # # if top level is aggregation, then do full V-cycles
  # else 
  #   topV = 1
  # end
  topV = 1

  # save current residual
  prev_resid = norm(lv[topV].A*lv[topV].x-lv[topV].b,opts.normtype)/norm(lv[topV].b,opts.normtype)
  # lv[1].opcount.L+=1

  # loop and do V-cycles
  for j=1:opts.max_num_vcycle

    
    resetMatvecCount(lv)
    println("Vcycle = ",j)
   
    # Do a V-Cycle
    
    # initialize worst recombination condition number to zero (using a vector so it passes by reference, only use first entry)
    vCycle!(lv[topV:end],opts,length(lv[topV:end]))

    # find norm of residual
    new_resid = norm(lv[topV].A*lv[topV].x-lv[topV].b,opts.normtype)/norm(lv[topV].b,opts.normtype)

    # compute V-cycle convergence factor
    conFac = new_resid/prev_resid
    println("Iteration ",j," V-cycle Residual = ",new_resid)
    println("Iteration ",j," V-cycle Convergence Factor = ",conFac)
    
    # record convergence history
    push!(stats.conFacHistory,conFac)
    # record recombination condition numbers
    # push!(stats.recombCondition,worstRecombCondition[1])
    # record variance of coefficients
    # push!(stats.alphaVar,worstAlphaVar[1])

    # terminate solver if residual is below tolerance, or residual increases
    # if new_resid < opts.tol || conFac < opts.tol  ||  new_resid > prev_resid
    if new_resid < opts.tol 
      # || conFac > opts.conv_cutoff  

      # if top level is elimination (we've been doing V-cycles below this level)
      # then we must bring solution up to top level

      if new_resid < opts.tol
        stats.converged = true
      end

      if string(lv[1].lvFunc) == "elim"

        # extend solution vector to top level
        lv[1].x = lv[1].P*lv[2].x + lv[1].Q*lv[1].b 
      #   alpha = (lv[1].x'*lv[1].rnk)/(lv[1].rnk'*lv[1].rnk) 
      # lv[1].x = lv[1].x + alpha[1]*lv[1].rnk 
        # lv[1].vecs.x = lv[1].vecs.x - mean(lv[1].vecs.x)
        # lv[2].opcount.P+=1
        new_resid = norm(lv[1].A*lv[1].x-lv[1].b,opts.normtype)/norm(lv[1].b,opts.normtype)
        println("new resid = ", new_resid)
        println("new resid = ", norm(lv[2].A*lv[2].x-lv[2].b,opts.normtype)/norm(lv[2].b,opts.normtype))

      end

      println("Final resid = ",new_resid)
      return 
    end

    prev_resid = new_resid

  end

end


function resetMatvecCount(lv::Array{level,1})
  for it in lv 
    it.opcount.L = 0
    it.opcount.R = 0
    it.opcount.P = 0
  end
end
# --------------------------------------------------------------
# Function Prints the statistics of the method.
# -------------------------------------------------------------- 
function PrintStatistics(fileName::AbstractString,lv::Array{level,1},stats::solveStats,opts::options )
  file = open(fileName,"w")
  numNodes = length(lv[1].x)
  numEdges = (countnz(lv[1].A) - numNodes)
  write(file,string("NumNodes: ",numNodes,"\n"))
  write(file,string("NumEdges: ",numEdges,"\n"))  
  write(file,string("NumLevels: ", length(lv),"\n"))
  write(file,string("rsvTime: ", stats.rsvT,"\n"))
  confacgeom = prod(stats.rsv_conFacHistory)^(1/length(stats.rsv_conFacHistory))
  write(file, string("rsv_rhoGeom: ", confacgeom,"\n"))
    write(file, string("rsv_numCycles: ", stats.rsv_numCycles,"\n"))

   error = norm(lv[1].A*lv[1].rnk,opts.normtype)
  write(file,string("|Lrsv|: ", error,"\n"))
  write(file,string("SetupTime: ", stats.setupT,"\n"))
  write(file,string("SolveTime: ", stats.solveT,"\n"))
  error = norm(lv[1].A*lv[1].x-lv[1].b,opts.normtype)/norm(lv[1].b,opts.normtype)
  write(file,string("|Lx-b|/|b|: ", error,"\n"))
  confacgeom = prod(stats.conFacHistory)^(1/length(stats.conFacHistory))
  write(file, string("rhoGeom: ", confacgeom,"\n"))

  history = calculateComplexHistory(lv)
  sums=0.0
  for i in history
    sums+= i.nnzA*i.numA; 
    sums+= i.nnzP*i.numP; 
    sums+= i.nnzR*i.numR; 
  end

  base = countnz(lv[1].A); 
  sumA = 0

  for l in lv
    sumA+= countnz(l.A); 
  end


  #Not sure I am calculating cycle comp right need to ask Geoff
  cycComp = sums/base; 
  write(file, string("cycComp: ", cycComp,"\n"))
  opComp = sumA/base;
  write(file, string("opComp: ", opComp,"\n"))
    write(file, string("efc: ", confacgeom^(1.0/cycComp),"\n"))

  write(file, string("NumVcycles: ", length(stats.conFacHistory),"\n"))
  write(file, string("converged: ", stats.converged,"\n"))

  close(file)

end


function calculateComplexHistory(lv::Array{level,1})
history_counts = Array(complexHistory,0)
  for it =1:length(lv)

      numA = lv[it].opcount.L
    
      numP = lv[it].opcount.P 
      numR = lv[it].opcount.R
    
    nnzA = countnz(lv[it].A);

    if it == 1
      nnzP = 0
      nnzR = 0
    else
    if lv[it-1].P != ones(2,2)
      nnzP = countnz(lv[it-1].P);
    else 
      nnzP = 0; 
    end
    if lv[it-1].R != ones(2,2)
      nnzR = countnz(lv[it-1].R);
    else 
      nnzR=0;
    end
  end
    println("numA = ",numA,"numP = ",numP," numR = ",numR," nnzA = ",nnzA, " nnzP = ",nnzP," nnzR = ",nnzR)

   if lv[it].lvFunc == agg
    # println(lv[it].lvFunc)

    history_counts  = vcat(history_counts, complexHistory(numA,numP,numR,nnzA,nnzP,nnzR))
  end
  end
  return history_counts

end



# ----------------------------------------------------------
# Function to complete a LAMG mu-cycle (aggregation & elimination)  ( opts.mu = 1  ==>  V-cycle,  opts.mu = 2  ==>  W-cycle )
# ----------------------------------------------------------
function vCycle!(lv::Array{level,1}, opts::options, num_levels::Int64)
    # println("lv = ",length(lv))

  #  nu1 relaxations on current grid level
  #  no need for pre-relaxation on elimination levels ( on elimination levels lv[1].lev.nu1 = 0 )
  if opts.nu1 > 0 &&size(lv[1].A,1) >10
    # println("lv[1].func = ",lv[1].lvFunc)
    #  println("size R = ", size(lv[1].A))
    #   println("size x", size(lv[1].x))
    #   println("size b = ",size(lv[1].b))
    opts.smoother(lv[1].A,lv[1].x,lv[1].b,lv[1].rnk,opts.nu1, opts.w)
    alpha = dot(lv[1].x,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk) 
    lv[1].x = lv[1].x - alpha[1]*lv[1].rnk 
    lv[1].opcount.L += opts.nu1

  end

  if length(lv)>1
    
    #  restrict to coarse grid
    if string(lv[1].lvFunc) == "agg"

      # for aggregation, restrict residual for new RHS
      # println("size R = ", size(lv[1].R))
      # println("size x", size(lv[1].x))
      # println("size A = ", size(lv[1].A))
      # println("size b = ",size(lv[1].b))
      lv[2].b = lv[1].R'*(lv[1].b-lv[1].A*lv[1].x)
      lv[1].opcount.L+=1
      lv[2].opcount.R+=1
      # for aggregation, zero initial guess
      lv[2].x = zeros(size(lv[2].b))
      


    elseif string(lv[1].lvFunc) == "elim"
# println("ELim")
#  println("size R = ", size(lv[1].R))
#       println("size x", size(lv[1].x))
#       println("size A = ", size(lv[1].A))
#       println("size b = ",size(lv[1].b))
      # for elimination, restrict RHS
      lv[2].b = lv[1].R*lv[1].b
      # for elimination, inject x on coarse points
      lv[2].x = lv[1].Rinj*lv[1].x
      # lv[2].opcount.R+=1
    
    else 
      println(string(lv[1].lvFunc))
      println("Invalid Level Type")
      lv = []
      return 
    
    end

    #################################################
    # recursive call mu-cycle
    vCycle!(lv[2:end],opts,num_levels)
    # println("lv = ",length(lv))
    #################################################
      
    # 3.) coarse grid correction
    if string(lv[1].lvFunc) == "agg"

      # do coarse grid correction
      lv[1].x = lv[1].x + lv[1].P*lv[2].x
alpha = (lv[1].x'*lv[1].rnk)/(lv[1].rnk'*lv[1].rnk) 
    lv[1].x = lv[1].x - alpha[1]*lv[1].rnk 
          lv[2].opcount.P+=1

    elseif string(lv[1].lvFunc) == "elim"

      # for elimination, project up solution
      lv[1].x = lv[1].P*lv[2].x + lv[1].Q*lv[1].b 
alpha = (lv[1].x'*lv[1].rnk)/(lv[1].rnk'*lv[1].rnk) 
    lv[1].x = lv[1].x - alpha[1]*lv[1].rnk 
          # lv[2].opcount.P+=1

    else 
      println("Invalid Level Type")
      lv = []
      return 

    end

  else


    # solve on coarsest grid

    #  TODO :  account for cost of coarse-grid solve (either direct or by relaxation)

    if size(lv[1].A,1) < opts.ncoarse 
      # Coarse Grid Direct Solve
      lv[1].x=(full(lv[1].A+ 1/size(lv[1].A,1)*ones(size(lv[1].A))))\(lv[1].b)
        # lv[1].x = lv[1].x - mean(lv[1].x) 

        # Removes Kernal
        alpha = dot(lv[1].x,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk) 
    lv[1].x = lv[1].x - alpha[1]*lv[1].rnk 
          lv[1].opcount.L+=1

    else 
      # if matrix is large, solve by relaxation (since ACF is small)
      ######### TODO :  pick nusolve more intelligently
      ######### for now, maxACF^16 = 0.3^16 < 10^-8
      nusolve = 16
      opts.smoother(lv[1].A,lv[1].x,lv[1].b,lv[1].rnk,nusolve, opts.w)
      alpha = dot(lv[1].x,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk) 
    lv[1].x = lv[1].x - alpha[1]*lv[1].rnk 
      lv[1].opcount.L+=nusolve
    end

  end

  # 4.) relax nu2 times on fine grid
  #  no need for post-relaxation on elimination levels ( on elimination levels lv[1].lev.nu2 = 0 )
  if opts.nu2 > 0
    opts.smoother(lv[1].A,lv[1].x,lv[1].b,lv[1].rnk,opts.nu2, opts.w)
aalpha = dot(lv[1].x,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk) 
    lv[1].x = lv[1].x - alpha[1]*lv[1].rnk 
        lv[1].opcount.L+=opts.nu2
    
  end

  # this does iterate recombination
  # if string(lv[1].lev.lvFunc) == "aggLAMG" && length(lv)!=num_levels
  if string(lv[1].lvFunc) == "agg"
    nprev = opts.nprev

    # save current vectors x and Ax into histories
    # NOTE:  we recycle the columns here, overwriting cyclicly as we go
    #  achieves two gains:  only store the number of history vectors we will continue to use & allows full pre-allocation
   
    lv[1].histx[:,rem( lv[1].histcount ,nprev) + 1 ] = lv[1].x 
    lv[1].histAx[:,rem( lv[1].histcount ,nprev) + 1 ] = lv[1].A*lv[1].x 

    lv[1].histcount+=1
    lv[1].opcount.L+=1   
    # solve least-squares problem for recombination coefficients alpha
    # save worst recombination condition number
    # worstRecombCondition[1] = max(worstRecombCondition[1],cond(lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)]))
    # use backslash (Julia uses pivoted QR here)
    alpha = lv[1].histAx[:,1:min(nprev,lv[1].histcount)] \ lv[1].b 
    # if worstAlphaVar[1] == 0.0
    #   worstAlphaVar[1] = var(alpha)
    # else 
    #   worstAlphaVar[1] = min(worstAlphaVar[1],var(alpha))
    # end 
    # try iterative solver
    # alpha = lsqr( lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)], lv[1].vecs.b )[1]

    # solve normal equations
    # alpha = At_mul_B(lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)], lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)]) \ At_mul_B(lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)], lv[1].vecs.b )
    # # let's see if we get a more stable result using SVD
    # U,S,V = svd(lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)])
    # alpha = V*diagm(S.^(-1))*U'*lv[1].vecs.b
    # println("SVs = ",S')
    # println("coefficients = ",alpha')

    # if length(alpha) > 1
    #   if var(alpha) > 10.0^(-16)
    # # compute new x from optimal linear combination of history
    lv[1].x = lv[1].histx[:,1:min(nprev,lv[1].histcount)]*alpha
    alpha = dot(lv[1].x,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk) 
    lv[1].x = lv[1].x - alpha[1]*lv[1].rnk 
    # # change last vectors in history to match recombination result
    lv[1].histx[:,rem( lv[1].histcount-1 ,nprev) + 1 ] = lv[1].x 
    lv[1].histAx[:,rem( lv[1].histcount-1 ,nprev) + 1 ] = lv[1].A*lv[1].x 
    #   end 
    # end 
    # # should we be counting this ?
    # lv[1].opcount.L += 1

  end

  # W-cycle recursive call (2nd recursion)
  #  Aly :  add details on how the W-cycle is performed (e.g. what levels we recurse on)



  return

end







function solve_rsv!(lv::Array{level,1}, stats::solveStats, opts::options)
  println(" lvfun = ", lv[1].lvFunc)

 # if top level is elimination, restrict down, and do V-cycles below top level
  # if string(lv[2].lvFunc) == "elim" 
  #   println(size(lv[1].R)) 
    
  #   topV = 2
  #   # for elimination, restrict RHS
  #   lv[2].b = lv[1].R*lv[1].b
  #   # for elimination, inject x on coarse points
  #   lv[2].x = lv[1].Rinj*lv[1].x
  #   # lv[2].opcount.R+=1

  # # if top level is aggregation, then do full V-cycles
  # else 
  #   topV = 1
  # end
  topV = 1
lv[1].b = zeros(size(lv[1].A,1))
  # save current residual
  prev_resid = norm(lv[topV].A*lv[topV].rnk-lv[topV].b,opts.normtype)
  # lv[1].opcount.L+=1
println("Begining resid = ",prev_resid)
  # loop and do V-cycles
  for j=1:opts.max_num_vcycle

    
    resetMatvecCount(lv)
    println("Vcycle = ",j)
   
    # Do a V-Cycle
    
    # initialize worst recombination condition number to zero (using a vector so it passes by reference, only use first entry)
    vCycle_rsv!(lv[topV:end],opts,length(lv[topV:end]))

    # find norm of residual
    new_resid = norm(lv[topV].A*lv[topV].rnk-lv[topV].b,opts.normtype)

    # compute V-cycle convergence factor
    conFac = new_resid/prev_resid
    println("Iteration ",j," V-cycle Residual = ",new_resid)
    println("Iteration ",j," V-cycle Convergence Factor = ",conFac)
    
    # record convergence history
    push!(stats.conFacHistory,conFac)
    # record recombination condition numbers
    # push!(stats.recombCondition,worstRecombCondition[1])
    # record variance of coefficients
    # push!(stats.alphaVar,worstAlphaVar[1])

    # terminate solver if residual is below tolerance, or residual increases
    # if new_resid < opts.tol || conFac < opts.tol  ||  new_resid > prev_resid
  		 new_rsvtol= dot(lv[1].A*lv[1].rnk,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk)
		println("rsvtol = ", new_rsvtol)

		if abs(new_rsvtol) < opts.rsv_tol  || conFac > opts.conv_cutoff  

      # if top level is elimination (we've been doing V-cycles below this level)
      # then we must bring solution up to top level

      if new_resid < opts.tol
        stats.converged = true
      end

      # if string(lv[1].lvFunc) == "elim"

      #   # extend solution vector to top level
      #   lv[1].x = lv[1].P*lv[2].x + lv[1].Q*lv[1].b 
      # #   alpha = (lv[1].x'*lv[1].rnk)/(lv[1].rnk'*lv[1].rnk) 
      # # lv[1].x = lv[1].x + alpha[1]*lv[1].rnk 
      #   # lv[1].vecs.x = lv[1].vecs.x - mean(lv[1].vecs.x)
      #   # lv[2].opcount.P+=1
      #   new_resid = norm(lv[1].A*lv[1].x-lv[1].b,opts.normtype)/norm(lv[1].b,opts.normtype)
      #   println("new resid = ", new_resid)
      #   println("new resid = ", norm(lv[2].A*lv[2].x-lv[2].b,opts.normtype)/norm(lv[2].b,opts.normtype))

      # end

      println("Final resid = ",new_resid)
      return 
    end

    prev_resid = new_resid

  end

end

function vCycle_rsv!(lv::Array{level,1}, opts::options, num_levels::Int64)
    println("lv = ",length(lv))

  #  nu1 relaxations on current grid level
  #  no need for pre-relaxation on elimination levels ( on elimination levels lv[1].lev.nu1 = 0 )
  if opts.nu1 > 0
    # println("Jacobi")
    #  println("size R = ", size(lv[1].A))
    #   println("size x", size(lv[1].x))
    #   println("size b = ",size(lv[1].b))
	println("residual before smoother = ",norm(lv[1].A*lv[1].rnk-lv[1].b,opts.normtype))


    opts.rsv_smoother(lv[1].A,lv[1].rnk,opts.nu1, opts.w)
    	println("residual after smoother = ",norm(lv[1].A*lv[1].rnk-lv[1].b,opts.normtype))

    lv[1].opcount.L += opts.nu1

  end

  if length(lv)>1

    #  restrict to coarse grid
    if string(lv[1].lvFunc) == "agg"

      # for aggregation, restrict residual for new RHS
      # println("size R = ", size(lv[1].R))
      # println("size x", size(lv[1].x))
      # println("size A = ", size(lv[1].A))
      # println("size b = ",size(lv[1].b))
      	println("norm of P before agg = ",lv[1].P[:,1])
      	# println("norm of P before agg = ",norm([spdiagm(lv[1].rnk)*lv[1].R][:,1]))

      lv[1].P = spdiagm(lv[1].rnk)*lv[1].R

		# scaleing by by diag(P'1)^(-1). This allows the coarser graph to be order 1 along the diagonal. 
		lv[1].P = spdiagm(lv[1].rnk)*lv[1].R*spdiagm(1./((spdiagm(lv[1].rnk)*lv[1].R)'*ones(size(lv[1].P,1))))
      	println("norm of P before agg = ",lv[1].P[:,1])
	# compute coarse-grid operator
		lv[2].A = lv[1].R'*lv[1].A*lv[1].P
      lv[1].opcount.L+=1
      lv[2].opcount.R+=1
      # for aggregation, zero initial guess
      


    elseif string(lv[1].lvFunc) == "elim"
# println("ELim")
#  println("size R = ", size(lv[1].R))
#       println("size x", size(lv[1].x))
#       println("size A = ", size(lv[1].A))
#       println("size b = ",size(lv[1].b))
      # for elimination, restrict RHS
      lv[2].A = lv[1].R*lv[1].A*lv[1].P
      # for elimination, inject x on coarse points
      lv[2].rnk = lv[1].Rinj*lv[1].rnk
      # lv[2].opcount.R+=1
    
    else 
      println(string(lv[1].lvFunc))
      println("Invalid Level Type")
      lv = []
      return 
    
    end

    #################################################
    # recursive call mu-cycle
    vCycle_rsv!(lv[2:end],opts,num_levels)
    # println("lv = ",length(lv))
    #################################################
      
    # 3.) coarse grid correction
    if string(lv[1].lvFunc) == "agg"

      # do coarse grid correction
      println("residual before agg = ",norm(lv[1].A*lv[1].rnk,opts.normtype))

      y = spdiagm(lv[1].R*spdiagm(1./(lv[1].R'*lv[1].rnk))*lv[2].rnk)*lv[1].rnk
			normy = norm(y,1)
			for i=1:length(lv[1].rnk)
				lv[1].rnk[i] = y[i]/normy; 
			end 
		println("residual after agg = ",norm(lv[1].A*lv[1].rnk,opts.normtype))

      lv[2].opcount.P+=1

    elseif string(lv[1].lvFunc) == "elim"

      # for elimination, project up solution
      y = lv[1].P*lv[2].rnk 
		normy = norm(y,1)
			for i=1:length(lv[1].rnk)
				lv[1].rnk[i] = y[i]/normy; 
			end 
          # lv[2].opcount.P+=1

    else 
      println("Invalid Level Type")
      lv = []
      return 

    end

  else


    # solve on coarsest grid

    #  TODO :  account for cost of coarse-grid solve (either direct or by relaxation)

    if size(lv[1].A,1) < opts.ncoarse 
      # Coarse Grid Direct Solve
      x1=(full(vcat([1 zeros(size(v[1].A,1)-1)'],v[1].A[2:end,:])))\([1;zeros(length(v[1].rnk)-1)])
      normy = norm(x1)

		for i=1:length(lv[1].rnk)
			lv[1].rnk[i] = lv[1].rnk[i]/normy; 
		end 
        # lv[1].x = lv[1].x - mean(lv[1].x) 

        # Removes Kernal
        
          lv[1].opcount.L+=1

    else 
      # if matrix is large, solve by relaxation (since ACF is small)
      ######### TODO :  pick nusolve more intelligently
      ######### for now, maxACF^16 = 0.3^16 < 10^-8
      nusolve = 16
       opts.rsv_smoother(lv[1].A,lv[1].rnk,nusolve, opts.w)
      lv[1].opcount.L+=nusolve
    end

  end

  # 4.) relax nu2 times on fine grid
  #  no need for post-relaxation on elimination levels ( on elimination levels lv[1].lev.nu2 = 0 )
  if opts.nu2 > 0
     opts.rsv_smoother(lv[1].A,lv[1].rnk,opts.nu2, opts.w)
        lv[1].opcount.L+=opts.nu2
    
  end

  # this does iterate recombination
  # if string(lv[1].lev.lvFunc) == "aggLAMG" && length(lv)!=num_levels
  if string(lv[1].lvFunc) == "agg"

    # save current vectors x and Ax into histories
    # NOTE:  we recycle the columns here, overwriting cyclicly as we go
    #  achieves two gains:  only store the number of history vectors we will continue to use & allows full pre-allocation
   
    lv[1].histx[:,rem( lv[1].histcount ,4) + 1 ] = lv[1].x 
    lv[1].histAx[:,rem( lv[1].histcount ,4) + 1 ] = lv[1].A*lv[1].x 

    lv[1].histcount+=1
    lv[1].opcount.L+=1   
    # solve least-squares problem for recombination coefficients alpha
    nprev = 4
    # save worst recombination condition number
    # worstRecombCondition[1] = max(worstRecombCondition[1],cond(lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)]))
    # use backslash (Julia uses pivoted QR here)
    alpha = lv[1].histAx[:,1:min(nprev,lv[1].histcount)] \ lv[1].b 
    # if worstAlphaVar[1] == 0.0
    #   worstAlphaVar[1] = var(alpha)
    # else 
    #   worstAlphaVar[1] = min(worstAlphaVar[1],var(alpha))
    # end 
    # try iterative solver
    # alpha = lsqr( lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)], lv[1].vecs.b )[1]

    # solve normal equations
    # alpha = At_mul_B(lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)], lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)]) \ At_mul_B(lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)], lv[1].vecs.b )
    # # let's see if we get a more stable result using SVD
    # U,S,V = svd(lv[1].history.Ax[:,1:min(nprev,lv[1].history.count)])
    # alpha = V*diagm(S.^(-1))*U'*lv[1].vecs.b
    # println("SVs = ",S')
    # println("coefficients = ",alpha')

    # if length(alpha) > 1
    #   if var(alpha) > 10.0^(-16)
    # # compute new x from optimal linear combination of history
    lv[1].x = lv[1].histx[:,1:min(nprev,lv[1].histcount)]*alpha
    alpha = dot(lv[1].x,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk) 
    lv[1].x = lv[1].x - alpha[1]*lv[1].rnk 
    # # change last vectors in history to match recombination result
    lv[1].histx[:,rem( lv[1].histcount-1 ,4) + 1 ] = lv[1].x 
    lv[1].histAx[:,rem( lv[1].histcount-1 ,4) + 1 ] = lv[1].A*lv[1].x 
    #   end 
    # end 
    # # should we be counting this ?
    # lv[1].opcount.L += 1

  end

  # W-cycle recursive call (2nd recursion)
  #  Aly :  add details on how the W-cycle is performed (e.g. what levels we recurse on)



  return

end
