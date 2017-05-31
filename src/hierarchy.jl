include("../src/level.jl")
include("../src/prolongation.jl")
include("../src/aggregation.jl")
include("../src/pruning.jl")
# using MATLAB


function bottomLevel()
	return
end



function hierarchy_test(A,At,rnk,x,b,opts::options,aggFunc)

	i = 1
	numTV = 1
	opcount = OPCOUNT(0,0,0)
	lv = [level(A,spzeros(2,2),spzeros(2,2),spzeros(2,2),spzeros(2,2),zeros(2),rnk,x,b,spzeros(2,2),spzeros(2,2),0,opcount,bottomLevel)]
	AA = copy(A)
	rnkc = copy(rnk)

		while  i< opts.maxLev && size(AA,1)>opts.ncoarse

		# for lvFunc in Any[aggFunc]
		for lvFunc in [elim, aggFunc]


			# if lvFunc == aggFunc
			# 	# # test relaxation speed and truncate elimination/aggregation if it is fast
			# 	#THIS IS THE WAY TOBY HAS IT< WE SHOULD USE TVTEST AS A TEST VECTOR!!!!
			# 	# THIS SHOULD ONLY BE DONT ON AGGREGATION STAGES!+
			# 	#This is an option in Toby, change later

			# 	# Minimun times we must smooth before calculating ACF
			# 	relaxAcfMinSweeps = 7
			# 	level_index = i	
			# 	# number of times we sweep
			# 	nu = relaxAcfMinSweeps+2*level_index

			# 	# how much we smooth for the test vectors. We will save the test vector that
			# 	# is created here for 
			# 	tvNu = opts.nu

			# 	# Discard these many initial iterations from ACF estimates
			# 	initial = 3

			# 	# random initial vector
			# 	xtest = rand(size(AA,1))
				
			# 	# zero right hand side
			# 	rtest = zeros(size(AA,1))

			# 	# smooth initial vector the minimal amount of time

			# 	# opts.smoother(AA,SGF, xtest, rtest, initial, w)
			# 	opts.smoother(AA,xtest, rtest, rnkc,initial, opts.w)

			# 	alpha = (xtest'*rnkc)/(rnkc'*rnkc) 
			# 	xtest = xtest- alpha[1]*rnkc

			# 	# make a copy of currect smoothed vector for our test vector
			# 	tvtest = copy(xtest)
			# 	# smooth test vector to get the right amount of sweeps so we can save
			# 	# opts.smoother(AA,SGF, tvtest, rtest, tvNu -initial,w)
			# 	opts.smoother(AA, tvtest, rtest, rnkc,tvNu -initial,opts.w)

			# 	alpha = (tvtest'*rnkc)/(rnkc'*rnkc) 
			# 	tvtest = tvtest- alpha[1]*rnkc


			# 	# we now smooth our current vector further to test for ACF. 
			# 	# smooth a total of nu = relaxAcfMinSweeps+2*level_index times. 
			# 	# Note that we have already smoothed opts.nu times. 
			# 	ytest = copy(tvtest)
			# 	# opts.smoother(AA,SGF, ytest, rtest, nu-tvNu, w)
			# 	opts.smoother(AA, ytest, rtest,  rnkc,nu-tvNu, opts.w)

			# 	alpha = (ytest'*rnkc)/(rnkc'*rnkc) 
			# 	ytest = ytest- alpha[1]*rnkc

			# 	# find the ratio between the more smoothed vector ytest 
			# 	# verse the initally smoothed vector xtest
			# 	ratio = norm(ytest,2)/norm(xtest,2)

			# 	#Calculate the ACF. 
			# 	relaxACF = ratio^(1/(nu-initial))

			# 	# Determine if relaxation is fast
			# 	# if relaxation is fast, truncate coarsening.  Solve by relaxation on coarsest (current) level
			# 	maxACF = 0.75
			# 	println("relaxACF = ", relaxACF)
			# 	if relaxACF<  maxACF
			# 		println("Relaxation is fast")
			# 		return lv 
			# 	end

	
			# end

			println("function  = ", lvFunc)
			if i < opts.maxLev && size(AA,1)>opts.ncoarse

				# println("size = ",size(AA,1))

			P,R,Q,Rinj,Ac,bc,xc,rnkc,ZFC = lvFunc(AA,sparse(transpose(AA)),zeros(size(AA,1)),lv[end].rnk,opts,numTV)
			if lvFunc == aggFunc
				numTV+=1
			end
			println("size A = ",size(AA,1))
				println("size Ac = ",size(Ac,1))


			i+=1
			if size(AA,1)>size(Ac,1) 
							# println("size A = " ,size(AA))

			# 	println("size Ac = " ,size(Ac))
			# println("lvfunc = " ,string(lvFunc))
			# 			println("size P = " ,size(P))

			# println("size R = " ,size(R))

					if size(P,1) > size(P,2)

						#set the coarse grid rsv
						# //rnkc = R'*lv[end].rnk

						lv[end].P = P
						lv[end].R = R
						lv[end].Rinj = Rinj
						lv[end].Q = Q
						if lvFunc == agg_rsv || lvFunc==agg_CR
							lv[end].lvFunc = agg
						else
							lv[end].lvFunc = lvFunc
						end

						lv[end].histx = spzeros(size(AA,1),opts.nprev)
						lv[end].histAx = spzeros(size(AA,1),opts.nprev)
						lv[end].ZFC = ZFC

						
					# else 
					# 	return lv
					end
				
					opcount = OPCOUNT(0,0,0)
					lv = vcat(lv, level(Ac,spzeros(2,2),spzeros(2,2),spzeros(2,2),spzeros(2,2),zeros(2),rnkc,xc,bc,spzeros(2,2),spzeros(2,2),0, opcount, bottomLevel))
					AA = copy(Ac)
			end
			else
				return lv
			end
			
			end
		
	
	end
	 return lv

end


function hierarchy_noElim(A,At,rnk,x,b,opts,aggFunc)

	i = 1
	numTV = 1
	opcount = OPCOUNT(0,0,0)
	lv = [level(A,spzeros(2,2),spzeros(2,2),spzeros(2,2),spzeros(2,2),zeros(2),rnk,x,b,spzeros(2,2),spzeros(2,2),0,opcount,bottomLevel)]
	AA = copy(A)
	rnkc = copy(rnk)

		while  i< opts.maxLev && size(AA,1)>opts.ncoarse

		# for lvFunc in Any[aggFunc]
		for lvFunc in [aggFunc]


			# if lvFunc == aggFunc
			# 	# # test relaxation speed and truncate elimination/aggregation if it is fast
			# 	#THIS IS THE WAY TOBY HAS IT< WE SHOULD USE TVTEST AS A TEST VECTOR!!!!
			# 	# THIS SHOULD ONLY BE DONT ON AGGREGATION STAGES!+
			# 	#This is an option in Toby, change later

			# 	# Minimun times we must smooth before calculating ACF
			# 	relaxAcfMinSweeps = 7
			# 	level_index = i	
			# 	# number of times we sweep
			# 	nu = relaxAcfMinSweeps+2*level_index

			# 	# how much we smooth for the test vectors. We will save the test vector that
			# 	# is created here for 
			# 	tvNu = opts.nu

			# 	# Discard these many initial iterations from ACF estimates
			# 	initial = 3

			# 	# random initial vector
			# 	xtest = rand(size(AA,1))
				
			# 	# zero right hand side
			# 	rtest = zeros(size(AA,1))

			# 	# smooth initial vector the minimal amount of time

			# 	# opts.smoother(AA,SGF, xtest, rtest, initial, w)
			# 	opts.smoother(AA,xtest, rtest, rnkc,initial, opts.w)

			# 	alpha = (xtest'*rnkc)/(rnkc'*rnkc) 
			# 	xtest = xtest- alpha[1]*rnkc

			# 	# make a copy of currect smoothed vector for our test vector
			# 	tvtest = copy(xtest)
			# 	# smooth test vector to get the right amount of sweeps so we can save
			# 	# opts.smoother(AA,SGF, tvtest, rtest, tvNu -initial,w)
			# 	opts.smoother(AA, tvtest, rtest, rnkc,tvNu -initial,opts.w)

			# 	alpha = (tvtest'*rnkc)/(rnkc'*rnkc) 
			# 	tvtest = tvtest- alpha[1]*rnkc


			# 	# we now smooth our current vector further to test for ACF. 
			# 	# smooth a total of nu = relaxAcfMinSweeps+2*level_index times. 
			# 	# Note that we have already smoothed opts.nu times. 
			# 	ytest = copy(tvtest)
			# 	# opts.smoother(AA,SGF, ytest, rtest, nu-tvNu, w)
			# 	opts.smoother(AA, ytest, rtest,  rnkc,nu-tvNu, opts.w)

			# 	alpha = (ytest'*rnkc)/(rnkc'*rnkc) 
			# 	ytest = ytest- alpha[1]*rnkc

			# 	# find the ratio between the more smoothed vector ytest 
			# 	# verse the initally smoothed vector xtest
			# 	ratio = norm(ytest,2)/norm(xtest,2)

			# 	#Calculate the ACF. 
			# 	relaxACF = ratio^(1/(nu-initial))

			# 	# Determine if relaxation is fast
			# 	# if relaxation is fast, truncate coarsening.  Solve by relaxation on coarsest (current) level
			# 	maxACF = 0.75
			# 	println("relaxACF = ", relaxACF)
			# 	if relaxACF<  maxACF
			# 		println("Relaxation is fast")
			# 		return lv 
			# 	end

	
			# end

			println("function  = ", lvFunc)
			if i < opts.maxLev && size(AA,1)>opts.ncoarse

				# println("size = ",size(AA,1))

			P,R,Q,Rinj,Ac,bc,xc,rnkc,ZFC = lvFunc(AA,AA',zeros(size(AA,1)),lv[end].rnk,opts,numTV)




			if lvFunc == aggFunc
				numTV+=1
			end
			println("size A = ",size(AA,1))
				println("size Ac = ",size(Ac,1))


			i+=1
			if size(AA,1)>size(Ac,1) 
							# println("size A = " ,size(AA))

			# 	println("size Ac = " ,size(Ac))
			# println("lvfunc = " ,string(lvFunc))
			# 			println("size P = " ,size(P))

			# println("size R = " ,size(R))

					if size(P,1) > size(P,2)

						#set the coarse grid rsv
						# //rnkc = R'*lv[end].rnk

						lv[end].P = P
						lv[end].R = R
						lv[end].Rinj = Rinj
						lv[end].Q = Q
						if lvFunc == agg_rsv || lvFunc==agg_CR
							lv[end].lvFunc = agg
						else
							lv[end].lvFunc = lvFunc
						end

						lv[end].histx = spzeros(size(AA,1),opts.nprev)
						lv[end].histAx = spzeros(size(AA,1),opts.nprev)
						lv[end].ZFC = ZFC

						
					# else 
					# 	return lv
					end
				
					opcount = OPCOUNT(0,0,0)
					lv = vcat(lv, level(Ac,spzeros(2,2),spzeros(2,2),spzeros(2,2),spzeros(2,2),zeros(2),rnkc,xc,bc,spzeros(2,2),spzeros(2,2),0, opcount, bottomLevel))
					AA = copy(Ac)
			end
			else
				return lv
			end
			
			end
		
	
	end
	 return lv

end



function find_rsv_no_save_hierarcy_RS(A::SparseMatrixCSC,At::SparseMatrixCSC, rnk::Vector, opts::options,stats::solveStats)
	# this will store the Tentative prolongation operaters that are found
	# adptively. Once a good rnk vector is found the final operators are 
	# store here so that we can use them for the solve of Ax = b

	prev_resid = norm(A*rnk,opts.normtype)
	 Ac,P,R,Q,Rinj = elimination(A,At,opts)



	# runs through vcycles
	for j=1:opts.max_num_vcycle

		# one adaptive vcycle updates tent. 
		rnkc = Rinj*rnk
		# rsv_no_save_hierarcy_RS!(Ac,rnkc,opts,1,1,4)
		 rsv_no_save_hierarcy_without_elim_RS!(Ac,rnkc,opts,1,1,4)
		   # println("size of tent = ", length(tent))
   
		rnk = P*rnkc
		new_resid = norm(A*rnk,opts.normtype)

		conFac = new_resid/prev_resid
		push!(stats.rsv_conFacHistory,conFac)

		println("residual = ", new_resid)
		println("Iteration ",j," V-cycle Convergence Factor = ",conFac)
		# checks if we have reach tolerance or the error is growing. 
		if new_resid < opts.tol  || conFac > opts.conv_cutoff 
			println("Final resid = ",new_resid)
			break
      	end

      	prev_resid = new_resid
	end


end


function elim(A, B,x,rnk,opts, i)
	# println("Elim size of A = ", size(A,1))

	Ae,P,R,Q,Rinj,ZFC = elimination(A,B,opts)
	println("size of R = ", size(Rinj))
	println("size of rnk = ", size(rnk))
	rnkc = Rinj*rnk
		xe = zeros(size(rnkc))

	be = zeros(size(xe))

	return P,R,Q,Rinj,Ae,be,xe,rnkc,ZFC

end



function agg(A,At,x,rnk,opts,i)
	
end


function agg_CR(A,At,x,rnk,opts,i)
	println("level = ",i)
	# println("Agg size of A = ", size(A,1))
	# println("size A = ", size(A,1))
	# println("time for aggregation")
	# @time	status = aggregateLAMG(A,opts, 1.0,2)
	# if (2+i-1>10)
	# 	d = 10
	# else
	# 	d = 2+i-1
	# end
	d=4

	gamma_l = 1.0
	if (2+i-1>10)
		d = 10
	else
		d = 2+i-1
	end
	status, rnk = aggregate_Affinity_RS_new(A,rnk,d,opts,i)
		# status, rnk = aggregate_degree_aware(A,rnk,d,opts,i)

	# println("norm  = ", norm(rnk))

	# status,rnk = v7_aggregate_CR(A,opts.theta,1,i,opts)

	# status = aggregate_Affinity_RS_new(A,rnk,d,opts)
	# status = aggregation_RS_new(A,rnk,opts.theta,1)
	# status = aggregation_RS(A,rnk,opts.theta,1)
	# Builds Tentative prolongation operator
	R = prolongation(status)

	# P,rnk = smoothP(A,R,rnk)
		println("size of R = ", size(R))

	# println("size of status = ", size(status))

	# println("zero  = ", length(find(x->x==0,rnk)))
	# rnk = spdiagm(vec(rnk))
	P = spdiagm(rnk)*R

	# scaleing by by diag(P'1)^(-1). This allows the coarser graph to be order 1 along the diagonal. 
	P = P*spdiagm(1./(P'*ones(size(P,1))))

	# compute coarse-grid operator
	Ac = R'*A*P
	# Ac =T'*A*spdiagm(x)*T*spdiagm(1./(T'*x))

	# println("min element = ", minimum(minimum(Ac)))
	# println("Ac = ", Ac)
	normxc= norm(R'*rnk,opts.normtype)

	rnkc = R'*rnk/normxc
	bc = zeros(size(R,1))
	xc = zeros(size(R,1))
	ZFC = zeros(size(xc))

	return P,R,spzeros(2,2), spzeros(2,2), Ac, bc,xc,rnkc,ZFC
end


function agg_rsv(A,At,x,rnk,opts,i)
	println("level = ",i)
	# println("Agg size of A = ", size(A,1))
	# println("size A = ", size(A,1))
	# println("time for aggregation")
	# @time	status = aggregateLAMG(A,opts, 1.0,2)
	if (2+i-1>10)
		d = 10
	else
		d = 2+i-1
	end
	# status = aggregate_Affinity_RS_new(A,rnk,d,opts)
	status = aggregation_RS_new(A,rnk,opts.theta,1)
	# status = aggregation_RS(A,rnk,opts.theta,1)
	# Builds Tentative prolongation operator
	R = prolongation(status)
	# println("R =",full(R))

	P = spdiagm(rnk)*R

	# scaleing by by diag(P'1)^(-1). This allows the coarser graph to be order 1 along the diagonal. 
	P = P*spdiagm(1./(P'*ones(size(P,1))))

	# P = (speye(size(A,1)) - .7*spdiagm(1./diag(A))*A)*spdiagm(rnk)*R
	# R = (speye(size(A,1)) - .7*spdiagm(1./diag(A))*A)*R

	# S = R'*spdiagm(diag(A))*P 
	# # println("R =",full(R))
	# # println("P =",full(P))
	# # println("S =",(S))

	# G  =R'*(A-spdiagm(diag(A)))*P 
	# # println("G =",(G))

	# # compute coarse-grid operator
	Ac = R'*A*P
	# 	# println("Ac =",(Ac))
	# # compute affinities directly sweeping though CSC format
	# for col = 1:size(Ac,1)
	# 	# pre-store column components
	# 	for j in Ac.colptr[col]:(Ac.colptr[col+1]-1)

	# 		if Ac.nzval[j]>0 && Ac.rowval[j]!=col
	# 			println(col)
	# 			println(Ac.rowval[j])

	# 			beta1 = (-1)*G[Ac.rowval[j],col]+ S[Ac.rowval[j],col]
	# 			beta2 = (-1)*G[col,Ac.rowval[j]]+ S[col,Ac.rowval[j]]
				
	# 			beta = max(beta1,beta2)
	# 			# println("beta ", beta)
	# 			Ac.nzval[j] += -beta
	# 			Ac[col,Ac.rowval[j]] += -beta
	# 			Ac[Ac.rowval[j],Ac.rowval[j]] += beta
	# 			Ac[col,col] += beta
	# 		end
	# 	end
	# end

	# 		println("Ac =",(Ac))
	# 		println(sum(Ac,1))
	# 					println(sum(Ac,2))



	# Ac =T'*A*spdiagm(x)*T*spdiagm(1./(T'*x))

	# println("min element = ", minimum(minimum(Ac)))
	# println("Ac = ", Ac)
	normxc= norm(R'*rnk,opts.normtype)

	rnkc = R'*rnk/normxc
	bc = zeros(size(R,1))
	xc = zeros(size(R,1))
	ZFC = zeros(2)


	return P,R,spzeros(2,2), spzeros(2,2), Ac, bc,xc,rnkc,ZFC
end




function find_rsv_RS!(A::SparseMatrixCSC, rnk::Vector, opts::options,stats::solveStats)
	# this will store the Tentative prolongation operaters that are found
	# adptively. Once a good rnk vector is found the final operators are 
	# store here so that we can use them for the solve of Ax = b

	println("in find")
	#P,R,Q,Rinj,Ac,bc,xc,rnkc = elim(A,A',zeros(size(A,1)),rnk,opts,1)

	prev_resid = norm(A*rnk,opts.normtype)

	prev_rsvtol= dot(A*rnk,rnk)/dot(rnk,rnk)
	# runs through vcycles
	for j=1:opts.rsv_max_num_vcycle
	# for j = 1:1

		# one adaptive vcycle updates tent. 
		# tent = [tentitive(zeros(2,2))]
		# println(" time for cycle = ", j)
		println("cycle = ", j)

		rsv_RS!(A,rnk,opts,1)
		   # println("size of tent = ", length(tent))
   

		new_resid = norm(A*rnk,opts.normtype)

		conFac = new_resid/prev_resid
		push!(stats.rsv_conFacHistory,conFac)

		println("residual = ", new_resid)
		# println("Iteration ",j," V-cycle Convergence Factor = ",conFac)
		# println("Opts.tol = ", opts.rsv_tol)
		# println("Is new_resid < opts.rsv = ", new_resid<opts.rsv_tol)
		# checks if we have reach tolerance or the error is growing. 
		new_rsvtol= dot(A*rnk,rnk)/dot(rnk,rnk)
		println("rsvtol = ", new_rsvtol)

		if abs(new_resid) < opts.rsv_tol  
			# || conFac > opts.conv_cutoff 
			println("in if")
			# println("Final resid = ",new_resid)

			#rnk = P*rnkc
			stats.rsv_numCycles = j; 
			return rnk

      	end

      	prev_resid = new_resid
	end

	# rnk = P*rnkc
	stats.rsv_numCycles = opts.rsv_max_num_vcycle; 
	return rnk

end

function rsv_RS!(A,rnk,opts, k::Int64)
	# smoother nu1 times
	# println("smoother k =", k)
	# opts.rsv_smoother(A,x,opts.rsv_nu1,opts.w)

	# n1 = size(A.nzval,1)
	println("size A = " ,opts.rsv_ncoarse )
	if size(A,1)>opts.rsv_ncoarse && k<opts.rsv_maxLev
		println("in rsv")
		# println("create status k= ", k)
		# if mod(k,2)==1
			# println("agg")
			opts.rsv_smoother(A,rnk,opts.rsv_nu1,opts.w)

			
			P,R,Q,Rinj,Ac,bc,xc,rnkc = agg_rsv(A,A',zeros(size(A,1)),rnk,opts,k)
			k+=1
			if size(Ac,1) !=size(A,1)
				

			rsv_RS!(Ac,rnkc,opts,k)
				    		# println("After solve = ", norm(Ac*xc))

			y = spdiagm(R*spdiagm(1./(R'*rnk))*rnkc)*rnk
			normy = norm(y,opts.normtype)
			for i=1:length(rnk)
				rnk[i] = y[i]/normy; 
			end 

			opts.rsv_smoother(A,rnk,opts.rsv_nu2,opts.w)
			# println("norm error  after agg = ", norm(A*rnk))
		# else

		# 	println("elim")
		# 	P,R,Q,Rinj,Ac,bc,xc,rnkc = elim(A,zeros(size(A,1)),rnk,opts,k)
		# 	k+=1
		# 	rsv_RS!(Ac,rnkc,opts,k)
		# 	rnk = P*rnkc
		# 	println("norm error  after agg = ", norm(A*rnk))



		end
	else

		# println("Residual before direct solve = ", norm(A*x))
		# Should be a direct Solve!!! what is the kernal tho? Fix this!
		# println("size A= ",size(A[2:end,:]))
		  nusolve = 16
    	opts.rsv_smoother(A,rnk,nusolve,opts.w)

  #   	x1=(full(vcat([1 zeros(size(A,1)-1)'],A[2:end,:])))\([1;zeros(length(rnk)-1)])

		# # println("Residual after direct solve = ", norm(x1))
		# normy = norm(x1,opts.normtype)

		# for i=1:length(rnk)
		# 	rnk[i] = rnk[i]/normy; 
		# end 

	 #    # lv[1].x = lv[1].x - mean(lv[1].x) 

	 #    # Removes Kernal
	    # alpha = (x1'*x)/(x'*x) 
	    # x = x1 - alpha[1]*x

		# opts.rsv_smoother(A,x,20,opts.w)
		    		# println("Residual after direct solve = ", norm(A*x))

	end		

end


# --------------------------------------------------------------
# Function that adaptively finds the Tentative prologation operator 
# for each level. Work is based from the paper for Adaptive Markov Chains. 
# -------------------------------------------------------------- 
function hierarchy_RS(tent,A,rnk,x,b,opts)
	i = 1
	opcount = OPCOUNT(0,0,0)
	lv = [level(A,spzeros(2,2),spzeros(2,2),rnk,x,b,spzeros(2,2),spzeros(2,2),0,opcount)]
	AA = copy(A)
	while i < opts.maxLev
		if size(AA,1)>opts.ncoarse && i < opts.maxLev

			# set R = T
			R = tent[i].T 
			# set P = daig(rsv)T
			P = spdiagm(lv[i].rnk)*tent[i].T

			# scaleing by by diag(P'1)^(-1). This allows the coarser graph to be order 1 along the diagonal. 
			P = P*spdiagm(1./(P'*ones(size(P,1))))
		
			# compute coarse-grid operator
			Ac = R'*AA*P

			# coarse grid right hand side will be replaced by restricted residual during solve
			# initialize with zeros
			bc = zeros(size(R,1))
			# start with zero initial guess on coarse levels
			xc = zeros(size(bc))

			if size(P,1) > size(P,2)

				#set the coarse grid rsv
				rnkc = R'*lv[end].rnk

				lv[end].P = P
				lv[end].R = R
				lv[end].histx = zeros(size(AA,1),4)
				lv[end].histAx = zeros(size(AA,1),4)

				
			end
			i +=1
		
			opcount = OPCOUNT(0,0,0)
			lv = vcat(lv, level(Ac,spzeros(2,2),spzeros(2,2),rnkc,xc,bc,spzeros(2,2),spzeros(2,2),0, opcount))
			AA = copy(Ac)
		else
			return lv
		end
	end
	 return lv

end


function find_rsv_RS_unweighted!(lv, rnk, opts::options,stats::solveStats)
		println("in Find")

	# this will store the Tentative prolongation operaters that are found
	# adptively. Once a good rnk vector is found the final operators are 
	# store here so that we can use them for the solve of Ax = b


	# P,R,Q,Rinj,Ac,bc,xc,rnkc = elim(A,A',zeros(size(A,1)),rnk,opts,1)

	prev_resid = norm(lv[1].A*rnk,opts.normtype)

	prev_rsvtol= abs(dot(lv[1].A*lv[1].rnk,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk))
	U = Vector{Any}
	rnkprev = copy(lv[1].rnk)
	# runs through vcycles
	for j=1:opts.rsv_max_num_vcycle
	# for j = 1:1

		# one adaptive vcycle updates tent. 
		# tent = [tentitive(zeros(2,2))]
		# println(" time for cycle = ", j)
		println("cycle = ", j)

		rsv_RS_unweighted!(lv,rnk,opts,1)
		# U = push!(U,lv[1].rnk-rnkprev)
		
		println("norm = ", norm(lv[1].A*lv[1].rnk,opts.normtype))
#  nprev = opts.nprev
     lv[1].histx[:,rem( lv[1].histcount ,opts.nprev) + 1 ] = lv[1].rnk
    # lv[1].histAx[:,rem( lv[1].histcount ,opts.nprev) + 1 ] = lv[1].A*lv[1].rnk - rnkprev
#     println("rem = ",nprev)
# println("x1  = ",norm( lv[1].A*lv[1].histx[:,1],1))
# println("x2  = ",norm( lv[1].A*lv[1].histx[:,2],1))

	
# 	lv[1].histcount+=1

# #     lv[1].opcount.L+=1 

# 	if   lv[1].histcount>1
# 		qin,rin = qr(full(lv[1].histx[:,1:min(opts.nprev,lv[1].histcount)]))
# 		qout,rout = qr(full(lv[1].A*qin))

# 		e,v = eig(rout'*rout)
# 		alpha = qin*v[:,findmin(e)[2]]
# 		alpha = alpha/norm(alpha,opts.normtype)
# 		if minimum(alpha)>0
# 			lv[1].rnk = alpha
# 		else
# 			println("Backing up")
# 			D = lv[1].histx[:,1:min(opts.nprev,lv[1].histcount)]
# 			D[:,rem( lv[1].histcount ,opts.nprev) + 1] = zeros(length(lv[1].rnk))
# 			qin,rin = qr(full(D))
# 			qout,rout = qr(full(lv[1].A*qin))

# 			e,v = eig(rout'*rout)
# 			alpha = qin*v[:,findmin(e)[2]]
# 			alpha = alpha/norm(alpha,opts.normtype)
# 			if minimum(alpha)>0
# 				lv[1].rnk = alpha
# 			end

# 		end
# 				println("norm = ", norm(lv[1].A*lv[1].rnk,opts.normtype))

# 	end


	# if   lv[1].histcount>1

	# 	x = Variable(min(opts.nprev,lv[1].histcount))
	# 	# The problem is to minimize ||Ax - b||^2 subject to x >= 0
	# 	# This can be done by: minimize(objective, constraints)
	# 	problem = minimize(sumsquares(lv[1].A*lv[1].histx[:,1:min(opts.nprev,lv[1].histcount)]*x), [sum(x) == 1.0]);
	# 	solve!(problem);
	# 	println(x.value)
	# 	alpha = lv[1].histx[:,1:min(opts.nprev,lv[1].histcount)]*vec(x.value)
	# 					println("norm = ", norm(lv[1].A*alpha,opts.normtype))

	# end



#     if   lv[1].histcount>1

#     	# for k = 1:2
#     		m = minimum(lv[1].histx)
#     	#   lv[1].histx[:,k][lv[1].histx[:,k] .<.1*m] =.1*m
#     	#    lv[1].histx[:,k] =  lv[1].histx[:,k]./norm( lv[1].histx[:,k],1)
#     	# end

#     	# B = At_mul_B(lv[1].histAx[:,1:min(nprev,lv[1].histcount)],lv[1].histA[:,1:min(nprev,lv[1].histcount)])
#     	# X = At_mul_B(lv[1].hist[:,1:min(nprev,lv[1].histcount)],lv[1].histA[:,1:min(nprev,lv[1].histcount)])

#     	x = Variable(min(nprev,lv[1].histcount))

# 		# The problem is to minimize ||Ax - b||^2 subject to x >= 0
# 		# This can be done by: minimize(objective, constraints)
# 		problem = minimize(sumsquares(lv[1].histAx[:,1:min(nprev,lv[1].histcount)]*x), [lv[1].histx[1:2,1:min(nprev,lv[1].histcount)]*x >0, sum(x) == 1.0]);
# 		solve!(problem);
# 		println(x.value)
# # x.value = x.value./norm(x.value,1)
# 		println("x before = ",norm(lv[1].A*lv[1].rnk,1))
# 		alpha = lv[1].histx[:,1:min(nprev,lv[1].histcount)]*vec(x.value)
# 		alpha = alpha./norm(alpha,1)

# 		println("alpha = ",norm(lv[1].A*alpha,1))
# 		if norm(lv[1].A*lv[1].rnk,1) >=norm(lv[1].A*alpha,1)
# 			lv[1].rnk = abs(lv[1].histx[:,1:min(nprev,lv[1].histcount)]*vec(x.value))
# 			# lv[1].histx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].rnk
# 			# lv[1].histAx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].A*lv[1].rnk
# 			println("x = ",norm(lv[1].A*lv[1].rnk,1))

# 		# else
# 		# 	println("in if")
# 		# 	println( lv[1].histx[find(x->x<0,alpha),1:2])
# 		# 	println( lv[1].histx[find(x->x<0,alpha),1:2]*vec(x.value))
# 		end

#     end
#     	 println("resid before = ", abs(dot(lv[1].A*lv[1].rnk,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk)))
#     	X = lv[1].histx[:,1:min(nprev,lv[1].histcount)]
#     	qin,rin = qr(full(X))
#     	qout,rout = qr(full(A*qin))
#     	B = At_mul_B(rout,rout,)
#     	f = zeros(min(nprev,lv[1].histcount),1)
#     	 @mput B 
# 		   @mput f
# # 		   @mput X
# # 		   @mput b
# # 		   @mput Aeq
# 		   		   @matlab alpha,flag, exitflag = quadprog(B,f)
# 		   @mget alpha
# 		   		   x1 = qin*alpha
# 		   		   if length(find(x->x<0, x1  ))==0
# 		   		   lv[1].rnk = x1/norm(x1,1)
# 		   		   println("BAD  = ", length(find(x->x<0,lv[1].rnk)))
# 		   		#    lv[1].histx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].rnk
# 			    # lv[1].histAx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].A*lv[1].rnk
# 			     println("resid after = ", norm(lv[1].A*lv[1].rnk,1))
# lv[1].rnk[lv[1].rnk .<=0] = .01
# opts.rsv_smoother(lv[1].A,lv[1].rnk,2,opts.w)
 # println("BAD  = ", length(find(x->x<0,lv[1].rnk)))
 # if  length(find(x->x<0,lv[1].rnk)) >0
 # 	println("STOP")
 # 	break
 # end
	# 	   		       	#  println("resid after = ", norm(lv[1].A*lv[1].rnk,1))

		   		       	# end


# 		   # @matlab alpha,flag, exitflag = quadprog(c,f,-X,b,Aeq,1.0);
# 		    @mget alpha
# 		    @mget exitflag
    # end
    # # solve least-squares problem for recombination coefficients alpha
#     if lv[1].histcount>1
#  			opts.rsv_smoother(lv[1].A,lv[1].rnk,1,opts.w)

# 		  println("resid before = ", norm(lv[1].A*lv[1].rnk,1))

# 		  B = At_mul_B(lv[1].A*lv[1].histx[:,1:min(nprev,lv[1].histcount)],lv[1].A*lv[1].histx[:,1:min(nprev,lv[1].histcount)])
# 		    println("b = ", size(B))
# 		  f = zeros(min(nprev,lv[1].histcount),1)
# 		  		  println("Aefq = ", size(f))

# 		  X = lv[1].histx[:,1:min(nprev,lv[1].histcount)]
# 		  b = zeros(length(lv[1].rnk),1)
# 		  # b = zeros(10,1)
# 		  C = inv(full(At_mul_B(X,X)))*B
# 		  Aeq = ones(1,min(nprev,lv[1].histcount))

# 		  println("b = ", size(b))
# 		  		  println("Aeq = ", size(X))
# @mput C
# 		   @mput B 
# 		   @mput f
# 		   @mput X
# 		   @mput b
# 		   @mput Aeq
# 		   		   @matlab alpha,flag, exitflag = quadprog(B,f)

# 		   # @matlab alpha,flag, exitflag = quadprog(c,f,-X,b,Aeq,1.0);
# 		    @mget alpha
# 		    @mget exitflag

# 		  alpha = alpha./norm(alpha,1)
# 		   	if length(find(x->x<0, lv[1].histx[:,1:min(nprev,lv[1].histcount)]*alpha  ))==0

# 		   		println("alapha = ", alpha)

# 			   	lv[1].rnk=lv[1].histx[:,1:min(nprev,lv[1].histcount)]*alpha
# 			   	lv[1].rnk = lv[1].rnk./norm(lv[1].rnk,1) 
# 		    	# # change last vectors in history to match recombination result
# 			    # lv[1].histx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].rnk
# 			    # lv[1].histAx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].A*lv[1].rnk

# 				println("resid after = ", norm(lv[1].A*lv[1].rnk,1))
# 			   	println("BAD  = ", length(find(x->x<0,lv[1].rnk)))
# 			   	itflag = 0
# 			   	 opts.rsv_smoother(lv[1].A,lv[1].rnk,1,opts.w)
# 			   	 println("BAD  = ", length(find(x->x<0,lv[1].rnk)))
# 			end
# 			# count-=1
# 		# end
# 	end
	


	   	# println(lv[1].rnk[find(x->x<0,lv[1].rnk)])
# 
		new_resid = norm(lv[1].A*lv[1].rnk,opts.normtype)

		conFac = new_resid/prev_resid
		# conFac = new_rsvtol /prev_rsvtol
		println("Confac = ", conFac)

		push!(stats.rsv_conFacHistory,conFac)

		println("residual = ", new_resid)
		# println("Iteration ",j," V-cycle Convergence Factor = ",conFac)
		# println("Opts.tol = ", opts.rsv_tol)
		# println("Is new_resid < opts.rsv = ", new_resid<opts.rsv_tol)
		# checks if we have reach tolerance or the error is growing. 
		new_rsvtol= abs(dot(lv[1].A*lv[1].rnk,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk))
		println("rsvtol = ", new_rsvtol)

		if abs(new_resid) < opts.rsv_tol  
			# || conFac > opts.conv_cutoff 
			println("in if")
			# println("Final resid = ",new_resid)

			# rnk = P*rnkc
			stats.rsv_numCycles = j; 
			return lv

      	end
      	
      	prev_resid = new_resid
	end

	# rnk = P *rnkc
	stats.rsv_numCycles = opts.max_num_vcycle; 
	return lv

end

function rsv_RS_unweighted!(lv::Array{level,1}, rnk::Vector,opts::options, num_levels::Int64)
    # println("lv = ",length(lv))

  #  nu1 relaxations on current grid level
  #  no need for pre-relaxation on elimination levels ( on elimination levels lv[1].lev.nu1 = 0 )
  if opts.nu1 > 0
  
  	# println("residual before = ", norm(lv[1].A*lv[1].rnk))
    opts.rsv_smoother(lv[1].A,lv[1].rnk,opts.rsv_nu1,opts.w)
            # gmres_smoother!(lv[1].A,lv[1].rnk,zeros(legnth(lv[1].rnk)),opts.rsv_nu1,opts.w)

      	# println("residual after = ", norm(lv[1].A*lv[1].rnk))

    # alpha = dot(lv[1].x,lv[1].rnk)/dot(lv[1].rnk,lv[1].rnk) 
    # lv[1].x = lv[1].x - alpha[1]*lv[1].rnk 
    # lv[1].opcount.L += opts.nu1

  end

  if length(lv)>1
    
    #  restrict to coarse grid
    if string(lv[1].lvFunc) == "agg"

		# P = spdiagm(rnk)*R
	# P = P*spdiagm(1./(P'*ones(size(P,1))))

		# scaleing by by diag(P'1)^(-1). This allows the coarser graph to be order 1 along the diagonal. 
		lv[1].P = spdiagm(lv[1].rnk)*lv[1].R*spdiagm(1./((spdiagm(lv[1].rnk)*lv[1].R)'*ones(size(lv[1].R,1))))

    	lv[2].A = lv[1].R'*lv[1].A*lv[1].P
      # for aggregation, zero initial guess
      lv[2].rnk= lv[1].R'*lv[1].rnk
      	normxc= norm(lv[2].rnk,opts.normtype)
      lv[2].rnk=lv[2].rnk./normxc
  #     normy = norm(lv[2].rnk,1)
		# for i=1:length(lv[2].rnk)
		# 	lv[2].rnk[i] = lv[2].rnk[i]/normy; 
		# end
      


    elseif string(lv[1].lvFunc) == "elim"

    # println("Size Ac = ", size(lv[2].A))

    nZ = sum(lv[1].ZFC .== 2)
    nF = sum(lv[1].ZFC .== 0)
    nC = size(lv[1].A,1)- nZ - nF

    # println("size nc = ", nC )
    perm = sortperm(lv[1].ZFC)
    # permmat = speye(length(perm))
    # permmat.rowval = perm
    AFFinv = diag(lv[1].A[lv[1].ZFC .== 0,lv[1].ZFC .== 0]).^(-1) 

 #    println("size AFF =",size(AFFinv))
	# println("size permmat =",size(permmat))
	# println("size nZ =",nZ)
	# println("size nC =",nC)
	# println("size nF =",nF)



   #  lv[1].Q = permute!(spdiagm( [ zeros( nZ ); AFFinv; zeros( nC ) ]),perm,perm
   #  # push!(Qii,Qi)
   #  lv[1].R = A_mul_Bt([ spzeros(nC,nZ) -lv[1].A[lv[1].ZFC .== 1, lv[1].ZFC  .== 0]*spdiagm(AFFinv) speye(nC) ],permmat)
    
   #  lv[1].Rinj = A_mul_Bt( [spzeros(nC,nZ+nF) speye(nC)], permmat)

   #  # println(size(lv[1].A))
   #  # println("size of P = ", size(lv[1].A[lv[1].ZFC .== 1, lv[1].ZFC.== 2]))
  	# lv[1].P = permmat*[ spzeros(nZ,nC );-spdiagm(AFFinv)*lv[1].A[lv[1].ZFC .== 0, lv[1].ZFC.== 1]; speye(nC) ]



    lv[1].Q = permute!(spdiagm( [ zeros( nZ ); AFFinv; zeros( nC ) ]),perm,perm)
    # push!(Qii,Qi)

   
    lv[1].R = permute!([ spzeros(nC,nZ) -lv[1].A[lv[1].ZFC .== 1, lv[1].ZFC .== 0]*spdiagm(AFFinv) speye(nC) ],round(Int64,[i for i in 1:nC]),perm)
    
    lv[1].Rinj = permute!( [spzeros(nC,nZ+nF) speye(nC)],round(Int64,[i for i in 1:nC]),perm)

    # println(size(lv[1].A))
    # println("size of P = ", size(lv[1].A[lv[1].ZFC .== 1, lv[1].ZFC.== 2]))
  	lv[1].P =permute!( [ spzeros(nZ,nC );-spdiagm(AFFinv)*lv[1].A[lv[1].ZFC .== 0, lv[1].ZFC.== 1]; speye(nC) ],perm,round(Int64,[i for i in 1:nC]))
    # R = A_mul_Bt( [spzeros(nC,nZ+nF) speye(nC)], permmat)
    # println("diff = ", norm(lv[1].R*lv[1].A*lv[1].P - lv[2].A,1))
    lv[2].A =lv[1].R*lv[1].A*lv[1].P

    	lv[2].rnk= lv[1].Rinj*lv[1].rnk
      normxc= norm(lv[2].rnk,opts.normtype)
      lv[2].rnk=lv[2].rnk./normxc
     	
      # lv[2].opcount.R+=1
    
    else 
      println(string(lv[1].lvFunc))
      println("Invalid Level Type")
      lv = []
      return 
    
    end

    #################################################
    # recursive call mu-cycle
    rsv_RS_unweighted!(lv[2:end],lv[2].rnk,opts,num_levels)
    # println("lv = ",length(lv))
    #################################################
      
    # 3.) coarse grid correction
    if string(lv[1].lvFunc) == "agg"
    	  	# println("residual before = ", norm(lv[1].A*lv[1].rnk,opts.normtype))

		y = spdiagm(lv[1].R*spdiagm(1./(lv[1].R'*lv[1].rnk))*lv[2].rnk)*lv[1].rnk
		normy = norm(y,opts.normtype)
		for i=1:length(lv[1].rnk)
			lv[1].rnk[i] = y[i]/normy; 
		end 

	# println("residual after = ", norm(lv[1].A*lv[1].rnk,opts.normtype))

      # do coarse grid correction
    elseif string(lv[1].lvFunc) == "elim"
    	  	# println("elim residual before = ", norm(lv[1].A*lv[1].rnk,opts.normtype))

      # for elimination, project up solution
     y = lv[1].P*lv[2].rnk 
		normy = norm(y,opts.normtype)
		for i=1:length(lv[1].rnk)
			lv[1].rnk[i] = y[i]/normy; 
		end 
		    	  	# println("elim residual after = ", norm(lv[1].A*lv[1].rnk,opts.normtype))

    else 
      println("Invalid Level Type")
      lv = []
      return 

    end

  else


    # solve on coarsest grid

    #  TODO :  account for cost of coarse-grid solve (either direct or by relaxation)
    	  	# println("solve residual before = ", norm(lv[1].A*lv[1].rnk,opts.normtype))

    if size(lv[1].A,1) < opts.ncoarse 
      # Coarse Grid Direct Solve
        nusolve = 16
    	opts.rsv_smoother(lv[1].A,lv[1].rnk,nusolve,opts.w)
        # lv[1].rnk=(full(vcat([1 zeros(size(lv[1].A,1)-1)'],lv[1].A[2:end,:])))\([1;zeros(length(lv[1].rnk)-1)])

    else 
      # if matrix is large, solve by relaxation (since ACF is small)
      ######### TODO :  pick nusolve more intelligently
      ######### for now, maxACF^16 = 0.3^16 < 10^-8
      nusolve = 16
    	opts.rsv_smoother(lv[1].A,lv[1].rnk,nusolve,opts.w)
   
    end
    	  	# println("solve residual before = ", norm(lv[1].A*lv[1].rnk,opts.normtype))

  end

  # 4.) relax nu2 times on fine grid
  #  no need for post-relaxation on elimination levels ( on elimination levels lv[1].lev.nu2 = 0 )
  if opts.nu2 > 0
    opts.rsv_smoother(lv[1].A,lv[1].rnk,opts.rsv_nu2,opts.w)
        # gmres_smoother!(lv[1].A,lv[1].rnk,zeros(legnth(lv[1].rnk)),opts.rsv_nu1,opts.w)


  end

  # this does iterate recombination
  # if string(lv[1].lev.lvFunc) == "aggLAMG" && length(lv)!=num_levels
#   if string(lv[1].lvFunc) == "agg"


#   	 nprev = opts.nprev
#      lv[1].histx[:,rem( lv[1].histcount ,nprev) + 1 ] = lv[1].rnk
#     lv[1].histAx[:,rem( lv[1].histcount ,nprev) + 1 ] = lv[1].A*lv[1].rnk
#     println("rem = ",rem( lv[1].histcount ,nprev))

#   	    if   lv[1].histcount>1

#     	# for k = 1:2
#     		m = minimum(lv[1].histx)
#     	#   lv[1].histx[:,k][lv[1].histx[:,k] .<.1*m] =.1*m
#     	#    lv[1].histx[:,k] =  lv[1].histx[:,k]./norm( lv[1].histx[:,k],1)
#     	# end

#     	# B = At_mul_B(lv[1].histAx[:,1:min(nprev,lv[1].histcount)],lv[1].histA[:,1:min(nprev,lv[1].histcount)])
#     	# X = At_mul_B(lv[1].hist[:,1:min(nprev,lv[1].histcount)],lv[1].histA[:,1:min(nprev,lv[1].histcount)])

#     	x = Variable(min(nprev,lv[1].histcount))

# 		# The problem is to minimize ||Ax - b||^2 subject to x >= 0
# 		# This can be done by: minimize(objective, constraints)
# 		problem = minimize(sumsquares(lv[1].histAx[:,1:min(nprev,lv[1].histcount)]*x), [lv[1].histx[:,1:min(nprev,lv[1].histcount)]*x >=.1*m,sum(abs(x))==1]);
# 		solve!(problem);
# 		println(x.value)

# 		println("x before = ",norm(lv[1].A*lv[1].rnk,1))
# 		alpha = lv[1].histx[:,1:min(nprev,lv[1].histcount)]*vec(x.value)

# 		println("alpha = ",norm(lv[1].A*alpha,1))
# 		if norm(lv[1].A*lv[1].rnk,1) >=norm(lv[1].A*alpha,1)
# 			lv[1].rnk = abs(lv[1].histx[:,1:min(nprev,lv[1].histcount)]*vec(x.value))
# 			lv[1].histx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].rnk
# 			lv[1].histAx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].A*lv[1].rnk
# 			println("x = ",norm(lv[1].A*lv[1].rnk,1))

# 		# else
# 		# 	println("in if")
# 		# 	println( lv[1].histx[find(x->x<0,alpha),1:2])
# 		# 	println( lv[1].histx[find(x->x<0,alpha),1:2]*vec(x.value))
# 		end

#     end
# end



#   	 nprev = opts.nprev
#      lv[1].histx[:,rem( lv[1].histcount ,nprev) + 1 ] = lv[1].rnk
#     lv[1].histAx[:,rem( lv[1].histcount ,nprev) + 1 ] = lv[1].A*lv[1].rnk
#     println("rem = ",rem( lv[1].histcount ,nprev))




#     lv[1].histcount+=1
#     lv[1].opcount.L+=1 

#   	  if   lv[1].histcount>1
# 			     println("resid before = ", norm(lv[1].A*lv[1].rnk,1))
#     	X = lv[1].histx[:,1:min(nprev,lv[1].histcount)]
#     	qin,rin = qr(full(X))
#     	qout,rout = qr(full(lv[1].A*qin))
#     	B = At_mul_B(rout,rout,)
#     	f = zeros(min(nprev,lv[1].histcount),1)
#     	 @mput B 
# 		   @mput f
# # 		   @mput X
# # 		   @mput b
# # 		   @mput Aeq
# 		   		   @matlab alpha,flag, exitflag = quadprog(B,f)
# 		   @mget alpha
# 		   		   x1 = qin*alpha
# 		   		   if length(find(x->x<0, x1  ))==0
# 		   		   lv[1].rnk = x1/norm(x1,1)
# 		   		   # println("BAD  = ", length(find(x->x<0,lv[1].rnk)))
# 		   		   lv[1].histx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].rnk
# 			    lv[1].histAx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].A*lv[1].rnk
# 			     println("resid after = ", norm(lv[1].A*lv[1].rnk,1))
# # opts.rsv_smoother(lv[1].A,lv[1].rnk,2,opts.w)
# # 		   		       	 println("resid after = ", norm(lv[1].A*lv[1].rnk,1))

# 		   		       	end


# # 		   # @matlab alpha,flag, exitflag = quadprog(c,f,-X,b,Aeq,1.0);
# # 		    @mget alpha
# # 		    @mget exitflag
#     end
 # nprev = opts.nprev
 #     lv[1].histx[:,rem( lv[1].histcount ,nprev) + 1 ] = lv[1].rnk
 #    lv[1].histAx[:,rem( lv[1].histcount ,nprev) + 1 ] = lv[1].A*lv[1].rnk
 #    println("rem = ",rem( lv[1].histcount ,nprev))




 #    lv[1].histcount+=1
 #    lv[1].opcount.L+=1   
 #    # # solve least-squares problem for recombination coefficients alpha
 #    if lv[1].histcount>1
 # 			opts.rsv_smoother(lv[1].A,lv[1].rnk,1,opts.w)

	# 	  println("resid before = ", norm(lv[1].A*lv[1].rnk,1))

	# 	  B = At_mul_B(lv[1].A*lv[1].histx[:,1:min(nprev,lv[1].histcount)],lv[1].A*lv[1].histx[:,1:min(nprev,lv[1].histcount)])
	# 	    println("b = ", size(B))
	# 	  f = zeros(min(nprev,lv[1].histcount),1)
	# 	  		  println("Aefq = ", size(f))

	# 	  X = lv[1].histx[1:2,1:min(nprev,lv[1].histcount)]
	# 	  # X = zeros(length(lv[1].rnk),min(nprev,lv[1].histcount))
	# 	  b = zeros(length(lv[1].rnk),1)
	# 	  b = zeros(2,1)
	# 	  # C = inv(full(At_mul_B(X,X)))*B
	# 	  Aeq = ones(1,min(nprev,lv[1].histcount))
	# 	  lb = zeros(min(nprev,lv[1].histcount),1)
	# 	ub = ones(min(nprev,lv[1].histcount),1)
	# 	  println("b = ", size(b))
	# 	  		  println("Aeq = ", size(X))

	# 	   @mput B 
	# 	   @mput f
	# 	   @mput X
	# 	   @mput b
	# 	   @mput Aeq
	# 	   @mput ub
	# 	   @mput lb
	# 	   		   # @matlab alpha,flag, exitflag = quadprog(2*B,f)

	# 	   @matlab alpha,flag, exitflag = quadprog(B,f,-X,b,Aeq,1.0);
	# 	    @mget alpha
	# 	    @mget exitflag

	# 	  # alpha = 1*alpha./norm(alpha,1)
	# 	   	if length(find(x->x<0, lv[1].histx[:,1:min(nprev,lv[1].histcount)]*alpha  ))==0 && sum(alpha)!=0

	# 	   		println("alapha = ", alpha)

	# 		   	lv[1].rnk=lv[1].histx[:,1:min(nprev,lv[1].histcount)]*alpha
	# 		   	lv[1].rnk = lv[1].rnk./norm(lv[1].rnk,1) 
	# 	    	# # change last vectors in history to match recombination result
	# 		    # lv[1].histx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].rnk
	# 		    # lv[1].histAx[:,rem( lv[1].histcount-1 ,opts.nprev) + 1 ] = lv[1].A*lv[1].rnk

	# 			println("resid after = ", norm(lv[1].A*lv[1].rnk,1))
	# 		   	println("BAD  = ", length(find(x->x<0,lv[1].rnk)))
	# 		   	itflag = 0
	# 		   	 opts.rsv_smoother(lv[1].A,lv[1].rnk,1,opts.w)
	# 		   	 println("BAD  = ", length(find(x->x<0,lv[1].rnk)))
	# 		end
	# 		# count-=1
	# 	# end
	# end
	
  # end

  # W-cycle recursive call (2nd recursion)
  #  Aly :  add details on how the W-cycle is performed (e.g. what levels we recurse on)



  return

end
