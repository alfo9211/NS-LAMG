
# function to compute test vectors
# function testVectorsLAMG(A::SparseMatrixCSC, d::Int64,)
function testVectorsLAMG(A::SparseMatrixCSC, d::Int64, opts::options,status::Array{Float64,1})

	# do nu Gauss-Seidel relaxation sweeps on AX=0
	# X is a size-d vector of test vector (this avoids unnecessary copies of TVs in aggregation routine)
	n=size(A,1);
	# X = Vector{Any}[.5*ones(n)]	# start each vector as random +1/-1
	# X[1][bitrand(n)] = -.5
	X = Vector{Float64}[ones(n)]	# start each vector as random +1/-1
	X[1][bitrand(n)] = -1
	X[1][status.==0] = 0
	# X[1] = X[1]/norm(X[1])

	# println(X)
	# do relaxation sweeps
	wjacobi!(A, X[1], 5, 2/3)
	# alpha = (X[1]'*rnk)/(rnk'*rnk) 
	# 			X[1] = X[1]- alpha[1]*rnk
	
	# println(X)
	for i = 2:d 

		# add vectors until you have d of them
		# X = push!(X,.5*ones(n))
		# X[i][bitrand(n)] = -.5
		X = push!(X,ones(n))
		X[i][bitrand(n)] = -1
		X[i][status.==0] = 0

	end
		for i = 2:d 
			wjacobi!(A, X[i], 5, 2/3)
				# alpha = (X[i]'*rnk)/(rnk'*rnk) 
				# X[i] = X[i]- alpha[1]*rnk
				# X[i] = X[i]/norm(X[i])
				# X[i] = X[i] - mean(X[i])
			end
			
			# println("rho = ",rho)
		# end
	# println("number of TV = ", length(X))

	
	# for i = 1:d
	# 	wjacobi!(A, X[i], 1, 2/3)
	# 	alpha = (X[i]'*rnk)/(rnk'*rnk) 
	# 	X[i] = X[i]- alpha[1]*rnk
		
	# end
	# rho = [norm(X[i]) for i=1:d]
	# println("rho = ",rho)

	return X
end

function aggregateLAMG(lv::level, opts::options, gamma_l::Float64, d::Int64)
	# build test vectors
	X = testVectorsLAMG(lv, d, opts)

	# alpha_max := g/gamma (see LAMG long paper section 3.4.3)
	a_max = opts.g/gamma_l

	# save an off diagonal only version of A 
	# TODO:  do this more elegantly
	Aod = lv.ops.A - spdiagm(diag(lv.ops.A))
	# run aggregation routine
	status=aggregate(lv.ops.A, Aod, X, a_max, opts.I_max);
	# return status vector which defines aggregates (same value => aggregate together)
	return status
end

function aggregate(A::SparseMatrixCSC, At::SparseMatrixCSC,Aod::SparseMatrixCSC, a_max::Float64, I_max::Int64,l,opts)
	# 1:  initialize variables
	# I_max=2; # doesn't seem like I_max should be hard coded
	B = Inf*ones(I_max); n=size(A,1);
	status = -1*ones(n)

	# 2:  initialize variables
	n_c = n; a = 1; i = 0; delta = 0.9;
	d = 3+l
	if d>10
		d = 10
	end

	
	println("d = ",d)
	A_v = A.colptr[2:end]-A.colptr[1:end-1]-1;
	A_u = At.colptr[2:end]-At.colptr[1:end-1]-1;
	A_v = A_v+A_u;	
	status[A_v .>= 6*median(A_v)] = 0;
	X = testVectorsLAMG(A, d, opts,status)
	E = [maximum([abs(X[i][j]) for i = 1:length(X)]) for j = 1:n]
		x = zeros(n)
	for j = 1:n
		x[j] = E[j]
	end
	x = x/norm(x,opts.normtype)
	opts.rsv_smoother(A,x,2,opts.w)

	
	# 3:  compute LAMG affinities (C) from test vectors (X)
	C,maxC = affinityLAMG(Aod, X);
	maxC = vec(maximum(C, 2));	
	D = transpose(C)
	maxD = vec(maximum(D, 2));


	# 5:  compute neighborhood sizes
	# A_v = number of neighbors of each vertex
	
	# status[E .>= (maximum(E) - median(E))/2+median(E)] = 0;

	println("number of seeds = ", length(find(x->x == 0, status)))

	# 6:  initialize status array
	S = zeros(size(status,1),0);
	# 1:  initialize variables
	# I_max=2; # doesn't seem like I_max should be hard coded
	aggregateSize=ones(n);
	
	
	println(" number of unclaimed nodes = ", length(find(x->x == -1,status)))

	a = n_c/n
	while (a >= a_max) && (i < I_max)
		# 7:  initialize aggregation sweep
		i += 1; 			
	

		# delta=0.6*delta;

		# 8:  main aggregation step
		status, n_c, aggregateSize, X = aggregationStage!(status, Aod, A, n_c,C,maxC,D,maxD,X,aggregateSize,delta,A_v, x,x);
		println(" number of unclaimed nodes = ", length(find(x->x == -1,status)))
		# X = testVectorsLAMG(A, d, opts,status)
		# 	# X = push!(X,rnk./norm(rnk))

		# E = [maximum([abs(X[i][j]) for i = 1:length(X)]) for j = 1:n]

		# C,maxC = affinityLAMG(Aod, X);
		# maxC = vec(maximum(C, 2));	

		# D = transpose(C)
		# maxD = vec(maximum(D, 2));

		println("n_c = ",n_c)
		delta=0.8*delta;
		# 9:  a = aggregation ratio (B to identify most aggressive aggregation which remains below maximum)
		a = n_c/n
		println("a = ",a)

		B[i]= 1-a + 2*a*(a>a_max)
		println("B = ", B)


		# 10:  concatentate status vector to status array
		S = [S status];

		# 11:  Set seed status to their own node number (to match associate nodes)
		S[status .== 0, i] = (1:length(status))[status .== 0];
		# Nate's addition:  set unaggregated nodes to their own aggregate
		S[status .== -1, i] = (1:length(status))[status .== -1];

		# println("Status = ", status)

	# 12:
	end

	i=findmin(B)[2];
	
	
	# 14:  return status
	return status,x
end
# function to compute affinity matrix C
function affinityLAMG(Aod::SparseMatrixCSC, X::Array{Array{Float64,1},1})
	C = copy(Aod)
	n = size(Aod,1)
	maxC = zeros(n)

	# compile vectors stored in X as a matrix Y
	Y = zeros(n,length(X))
	for i = 1:length(X)
		Y[:,i] = X[i]
	end

	# compute affinities directly sweeping though CSC format
	for col = 1:n 
		# pre-store column components
		Ycolvec = vec(Y[col,:])
		normYcolvecsq = dot(Ycolvec,Ycolvec)
		for j in C.colptr[col]:(C.colptr[col+1]-1)
			# pre-store row components
			Yrowvec = vec(Y[C.rowval[j],:])
			normYrowvecsq = dot(Yrowvec,Yrowvec)
			# compute LAMG affinities
			C.nzval[j] = ( dot( Yrowvec, Ycolvec ) )^2 / ( normYrowvecsq * normYcolvecsq )
			# generate maxC as we go to save time
			if C.nzval[j] > maxC[col]
				maxC[col] = C.nzval[j]
			end
		end
	end

	return C,maxC
end



function aggregate2(A::SparseMatrixCSC, At::SparseMatrixCSC,Aod::SparseMatrixCSC,  a_max::Float64, I_max::Int64,l,opts)
	println("in agg2")
	# 1:  initialize variables
	# I_max=2; # doesn't seem like I_max should be hard coded
	B = Inf*ones(I_max); n=size(A,1);
	aggregateSize=ones(n);
	status = -1*ones(n)

	# 2:  initialize variables
	n_c = n; a = 1; i = 0; delta = 0.9;

	# 3:  compute LAMG affinities (C) from test vectors (X)
	# d=10

	d = 3+l
	
	println("d = ", d)
	if d>10
		d = 10
	end
		
	
		
	A_v = A.colptr[2:end]-A.colptr[1:end-1]-1;
	A_u = At.colptr[2:end]-At.colptr[1:end-1]-1;
	A_v = A_v+A_u;
	X = testVectorsLAMG(A, d, opts,status)
	E = [maximum([abs(X[i][j]) for i = 1:length(X)]) for j = 1:n]
	x = zeros(n)
	for j = 1:n
		x[j] = E[j]
	end
	x = x/norm(x)


	status[A_v .>= 8*median(A_v)] = 0;

	println("sizex =", size(x))

		C,maxC = affinityLAMG(Aod, X);
		maxC = vec(maximum(C, 2));

		D = transpose(C)
		maxD = vec(maximum(D, 2));
		strong_neighbors = ones(Bool,length(C.nzval))
		# maxC = maximum(C, 1);

		for k = 1:n
	   		for j in C.colptr[k]:(C.colptr[k+1]-1)
	   			if C.nzval[j] >= delta*max(maxC[C.rowval[j]], maxC[k])
	   				# strong_neighbors == true for strong connections
	   				# thus, no change necessary here 
	   			else 
	   				# strong_neighbors == false for weak connections
	   				strong_neighbors[j] = false
	   			end
	   		end
		end

		strong_neighbors_row = ones(Bool,length(C.nzval))
			# maxC = maximum(C, 1);

			for k = 1:n
		   		for j in D.colptr[k]:(D.colptr[k+1]-1)
		   			if D.nzval[j] >= delta*max(maxD[D.rowval[j]], maxD[k])
		   				# strong_neighbors == true for strong connections
		   				# thus, no change necessary here 
		   			else 
		   				# strong_neighbors == false for weak connections
		   				strong_neighbors_row[j] = false
		   			end
		   		end
			end



	alphaInc = .1*n
	i=0
	while (n_c/n>a_max)	&& (i < 3)
		println( "i = ",i )
# while  (i < 3)
		i+=1	
		
		
			U = (status .!=0);
			ncInc = 0
		

			for u=sortperm(E,rev=true)
						if U[u] 

							rows_all = D.rowval[ D.colptr[u]:D.colptr[u+1]-1 ] 
							strong_mask = strong_neighbors_row[ D.colptr[u]:D.colptr[u+1]-1 ] 
							Sval=D.nzval[D.colptr[u]:D.colptr[u+1]-1][strong_mask];
							N_u = rows_all[strong_mask];
							t = N_u
							t=t[status[t] .== -1]

							rows_all   = C.rowval[ C.colptr[u]:C.colptr[u+1]-1 ] 
							strong_mask = strong_neighbors[ C.colptr[u]:C.colptr[u+1]-1 ] 
							Sval=C.nzval[C.colptr[u]:C.colptr[u+1]-1][strong_mask];
							N_u = rows_all[strong_mask];
							t1 = N_u

					
							t1=t1[status[t1] .== -1]
							s = unique(vcat(t,t1))

							
							flag = 0
													# 4-6:
							if length(s)>0 && status[u] > -1

									for k = s
										# if E[u]>=E[k]
											status[k] = u
											# # status[t[findmin(Sval)[2]]]  = u
											# n_c -=1
											# aggregateSize[u] +=length(t)
											# status[u] = 0
											U[k] = false
											# status[N_u[ status[N_u].==-1][findmin(s)[2]]] = u
										# end

									end
									status[u] = 0
									ncInc+=1

							elseif length(s)>0
								status[u] = 0
								ncInc+=1
								n_c+=1
								for k = s
										# if E[u]>=E[k]
											status[k] = u
											# # status[t[findmin(Sval)[2]]]  = u
											# n_c -=1
											# aggregateSize[u] +=length(t)
											# status[u] = 0
											U[k] = false
											# status[N_u[ status[N_u].==-1][findmin(s)[2]]] = u
										# end

									end


							end

							if ncInc >= alphaInc
								# println("break if ")
								break
							end
							# println("ncInc = ", ncInc)
							
						end
						if ncInc >alphaInc
							# println("break inside for ")
							break
						end

					
			end
			# status[status.!=0]=-1
			# delta = .9*delta
			X = testVectorsLAMG(A, d, opts,status)
			E = [maximum([abs(X[i][j]) for i = 1:length(X)]) for j = 1:n]
		# 		C,maxC = affinityLAMG(Aod, X);
		# maxC = vec(maximum(C, 2));

		# D = transpose(C)
		# maxD = vec(maximum(D, 2));
		# strong_neighbors = ones(Bool,length(C.nzval))
		# # maxC = maximum(C, 1);

		# for k = 1:n
	 #   		for j in C.colptr[k]:(C.colptr[k+1]-1)
	 #   			if C.nzval[j] >= delta*max(maxC[C.rowval[j]], maxC[k])
	 #   				# strong_neighbors == true for strong connections
	 #   				# thus, no change necessary here 
	 #   			else 
	 #   				# strong_neighbors == false for weak connections
	 #   				strong_neighbors[j] = false
	 #   			end
	 #   		end
		# end

		# strong_neighbors_row = ones(Bool,length(C.nzval))
		# 	# maxC = maximum(C, 1);

		# 	for k = 1:n
		#    		for j in D.colptr[k]:(D.colptr[k+1]-1)
		#    			if D.nzval[j] >= delta*max(maxD[D.rowval[j]], maxD[k])
		#    				# strong_neighbors == true for strong connections
		#    				# thus, no change necessary here 
		#    			else 
		#    				# strong_neighbors == false for weak connections
		#    				strong_neighbors_row[j] = false
		#    			end
		#    		end
		# 	end


	end
	status[status.!=0] = -1

	println("alpha = ", n_c/n)
	# n_c = n
	a = n_c/n
	i = 0
	while (a >= a_max) && (i < I_max)
		i+=1
		# 7:  initialize aggregation sweep
		delta=0.8*delta;
		status, n_c, aggregateSize, X = aggregationStage!(status, Aod, A, n_c,C,maxC,D,maxD,X,aggregateSize,delta,A_v, x,rnk);
		println(" number of unclaimed nodes = ", length(find(x->x == -1,status)))
		a = n_c/n
		# if (a>=a_max)
			# delta=0.8*delta;
			# X = testVectorsLAMG(A, d, opts,status)
			# E = [maximum([abs(X[i][j]) for i = 1:length(X)]) for j = 1:n]
			# C,maxC = affinityLAMG(Aod, X);	
		# 	maxC =  vec(maximum(C, 2));
		# end
		println("n_c = ",n_c)
		# 9:  a = aggregation ratio (B to identify most aggressive aggregation which remains below maximum)
		a = n_c/n
		println("a = ",a)

		# B[i]= 1-a + 2*a*(a>a_max)
		println("B = ", B)


		# 10:  concatentate status vector to status array


		# println("Status = ", status)

	# 12:
	end

	i=findmin(B)[2];

	println("number of seeds = ", length(find(x->x == 0, status)))
	
	status[status .== 0] = (1:length(status))[status .== 0];
		# Nate's addition:  set unaggregated nodes to their own aggregate
		status[status .== -1] = (1:length(status))[status .== -1];
	# 14:  return status
	return status, x
end
# function to compute affinity matrix C
function affinityLAMG(Aod::SparseMatrixCSC, X::Array{Array{Float64,1},1})
	C = copy(Aod)
	n = size(Aod,1)
	maxC = zeros(n)

	# compile vectors stored in X as a matrix Y
	Y = zeros(n,length(X))
	for i = 1:length(X)
		Y[:,i] = X[i]
	end

	# compute affinities directly sweeping though CSC format
	for col = 1:n 
		# pre-store column components
		Ycolvec = vec(Y[col,:])
		normYcolvecsq = dot(Ycolvec,Ycolvec)
		for j in C.colptr[col]:(C.colptr[col+1]-1)
			# pre-store row components
			Yrowvec = vec(Y[C.rowval[j],:])
			normYrowvecsq = dot(Yrowvec,Yrowvec)
			# compute LAMG affinities
			C.nzval[j] = ( dot( Yrowvec, Ycolvec ) )^2 / ( normYrowvecsq * normYcolvecsq )
			# generate maxC as we go to save time
			if C.nzval[j] > maxC[col]
				maxC[col] = C.nzval[j]
			end
		end
	end

	return C,maxC
end


# Algorithm 4
function aggregationStage!(status::Vector{Int64}, Aod::SparseMatrixCSC, A::SparseMatrixCSC, n_c::Int64,C::SparseMatrixCSC,maxC::Vector{Float64},D::SparseMatrixCSC,maxD::Vector{Float64},X::Array{Array{Float64,1},1},aggregateSize::Vector{Int64},delta::Float64,A_v::Vector{Float64}, E,rnk)
n=size(A,1);
medE = median(E)

	# 1:  find strong neighbors 
	strong_neighbors = ones(Bool,length(C.nzval))
	# maxC = maximum(C, 1);

	for i = 1:n
   		for j in C.colptr[i]:(C.colptr[i+1]-1)
   			if C.nzval[j] >= delta*max(maxC[C.rowval[j]], maxC[i])
   				# strong_neighbors == true for strong connections
   				# thus, no change necessary here 
   			else 
   				# strong_neighbors == false for weak connections
   				strong_neighbors[j] = false
   			end
   		end
	end


	strong_neighbors_row = ones(Bool,length(C.nzval))
	# maxC = maximum(C, 1);

	for i = 1:n
   		for j in D.colptr[i]:(D.colptr[i+1]-1)
   			if D.nzval[j] >= delta*max(maxD[D.rowval[j]], maxD[i])
   				# strong_neighbors == true for strong connections
   				# thus, no change necessary here 
   			else 
   				# strong_neighbors == false for weak connections
   				strong_neighbors_row[j] = false
   			end
   		end
	end

	# 2:  identify undecided nodes
	U = (status .== -1);

	# 3:  find undecided nodes with strong neighbors
	for u = 1:n
		if U[u]
			# purge zero values (weak neighbors) manually here
			#   do we need intersection with A_v ???
			rows_all   = C.rowval[ C.colptr[u]:C.colptr[u+1]-1 ] 
			strong_mask = strong_neighbors[ C.colptr[u]:C.colptr[u+1]-1 ] 
			Sval=C.nzval[C.colptr[u]:C.colptr[u+1]-1][strong_mask];
			N_u = rows_all[strong_mask];


			rows_all_row   = D.rowval[ D.colptr[u]:D.colptr[u+1]-1 ] 
			strong_mask_row = strong_neighbors_row[ D.colptr[u]:D.colptr[u+1]-1 ] 
			Sval_row=D.nzval[D.colptr[u]:D.colptr[u+1]-1][strong_mask_row];
			N_u_row= rows_all_row[strong_mask_row];
			# Sval=Sval[status[S] .<= 0]
			# S=S[status[S] .<= 0]

			# 4-6:
			if length(N_u)>0 && status[u] == -1
				# 7: find best Seed using fancy energy inflation minimizing search
				# s = bestSeedSimple(A, Aod, X, aggregateSize, N_u, C, u, status);
				s=bestSeedSimple(A,Sval,Sval_row,X,aggregateSize,N_u,N_u_row,C,u,status,E,medE)
				# s_min=bestSeedSimple_row(A,Sval,Sval_row,X,aggregateSize,N_u,N_u_row,C,u,status,rnk,medE)

				# 8: if a seed is found, then do some aggregation
					if s > 0.5 
						# || s_min>.5
						#9: update status vector
						# if s_min==0
							# if E[s]>E[u]
								status[s]=0; status[u]=s; n_c=n_c-1;
								U[s] = false
								U[u] = false
								for i = 1:length(X)
									X[i][u] =  X[i][s]
									E[u] = E[s]
								end
							# elseif status[s]!=0
							# 	status[u]=0; status[s]=u; n_c=n_c-1;
							# 	U[s] = false
							# 	U[u] = false

							# 	for i = 1:length(X)
							# 		X[i][u] =  X[i][s]
							# 		E[u] = E[s]
							# 	end

							# end
						# elseif s==0
						
						# 	if E[s_min]<=E[u] && status[s_min]!=0
								
						# 		status[u]=0; status[s_min]=u; n_c=n_c-1;
						# 		U[s] = false
						# 		for i = 1:length(X)
						# 			X[i][s_min] = X[i][u]
						# 			E[s_min] = E[u]
						# 		end
							# end
						# else
						# 	if E[s_min]<=E[u] && status[s_min]!=0
						# 		U[s] = false
						# 		status[u]=0; status[s_min]=u; n_c=n_c-1;
						# 		for i = 1:length(X)
						# 			X[i][s_min] = X[i][u]
						# 			E[s_min] = E[u]
						# 		end
						# 	elseif E[s]>E[u]
						# 		U[s] = false
						# 		status[s]=0; status[u]=s; n_c=n_c-1;
						# 		for i = 1:length(X)
						# 			X[i][u] =  X[i][s]
						# 			E[u] = E[s]
						# 		end
						# 	end

						# end
					end

						 # status[u]=s; n_c=n_c-1;
					# else 
						# status[u] = status[s]; n_c=n_c-1;
					# end
					# println(" In stage n_c = ",n_c)

					#10: update TV values on new aggregate
					# for i = 1:length(X)
					# 	X[i][u] = X[i][s]
					# 	E[u] = E[s]
					# end

					#11: update aggregate size
					# aggregateSize[u] += 1; aggregateSize[s] += 1;
				
			end 
		end
	end


return status, n_c, aggregateSize, X


end

function bestSeedSimple(A::SparseMatrixCSC,Sval::Vector,Sval_row::Vector,X::Array{Array{Float64,1},1},aggregateSize::Vector{Float64},N_u::Vector{Float64},N_u_row::Vector{Float64},C::SparseMatrixCSC,u::Int64,status::Vector{Float64},E,medE)
	# 1:  identify strong neighbors which are seeds or undecided
	S= N_u 
	Eval=E[S]
	Sval = Eval[status[S] .<= 0]
	S=S[status[S] .<= 0]

	S_row= N_u_row 
	Eval_row =E[S_row]
	Sval_row = Eval_row[status[S_row ] .<= 0]
	S_row =S_row[status[S_row ] .<= 0]

	
	# 2:
	if isempty(S)  && isempty(S_row)
		# let 0 indicate NotFound
		return 0
	elseif isempty(S)
		return  S_row[findmax(Sval_row)[2]]
	elseif isempty(S_row)
		return S[findmax(Sval)[2]]
	else
		rowS = S_row[findmax(Sval_row)[2]]
		colS = S[findmax(Sval)[2]]
		# # println( "rowS = ", rowS)
		# # println( "colS = ", colS)
		# if E[rowS] > E[colS]
		if findmax(Sval_row)[1] >findmax(Sval)[1]
			return rowS
		else
			return colS
		end 
	end
end


function bestSeedSimple_row(A::SparseMatrixCSC,Sval::Vector,Sval_row::Vector,X::Array{Array{Float64,1},1},aggregateSize::Vector,N_u::Vector{Float64},N_u_row::Vector{Float64},C::SparseMatrixCSC,u::Int64,status::Vector{Float64},E,medE)
	# 1:  identify strong neighbors which are seeds or undecided
	S= N_u 
	Eval=E[S]
	Sval = Eval[status[S] .<= 0]
	S=S[status[S] .<= 0]

	S_row= N_u_row 
	Eval_row =E[S_row]
	Sval_row = Eval_row[status[S_row ] .<= 0]
	S_row =S_row[status[S_row ] .<= 0]

	
	# 2:
	if isempty(S)  && isempty(S_row)
		# let 0 indicate NotFound
		return 0
	elseif isempty(S)
		return  S_row[findmin(Sval_row)[2]]
	elseif isempty(S_row)
		return S[findmax(Sval)[2]]
	else
		rowS = S_row[findmin(Sval_row)[2]]
		colS = S[findmin(Sval)[2]]
		# println( "rowS = ", rowS)
		# println( "colS = ", colS)
		if E[rowS] < E[colS]
			return rowS
		else
			return colS
		end 
	end
end



# Algorithm 6 - Better bestSeed algorithm
function bestSeed(A::SparseMatrixCSC,Aod::SparseMatrixCSC,X::Array{Array{Float64,1},1},aggregateSize::Vector{Int64},N_u::Vector{Float64},C::SparseMatrixCSC,u::Int64,status::Vector{Int64})

	# 1:  identify strong neighbors which are seeds or undecided
	# S = N_u[status[N_u] .<= 0]
	# the following double loop is faster than the elementwise compare above
	S_count = 0
	for i = 1:length(N_u)
		if status[N_u[i]] <= 0
			S_count += 1
		end
	end

	# 2:
	# if isempty(S)
	# skip second loop if it won't find anything to put in S
	if S_count == 0
		# let 0 indicate NotFound
		return 0
	else 

		# second loop
		S = zeros(Int,S_count)
		S_count = 0
		for i = 1:length(N_u)
			if status[N_u[i]] <= 0
				S_count += 1
				S[S_count] = N_u[i]
			end
		end

		# 3: compute TV energy ratios
		K = length(X)
		auu = Auu(A,u)
		# F = zeros(K,1)
		# Ctk = zeros(length(S),K)
		# Ctk_divF = zeros(length(S),K)
		qt = zeros(length(S))

		####### merged loops for steps 4 & 5, and switched order for loop #5
		for k = 1:K 
			# pre-compute Bux() and Cux() and pass them to Eu() several times (since they don't depend on y)
			# Bu = Bux( Aod, u, X[k] )[1]
			# Cu = Cux( Aod, u, X[k] )[1]
			Bu,Cu = BuxCux( Aod, u, X[k] )
			
			# 4:
			# F[k] = Eu( A, Aod, u, X[k], Bu/auu, Bu, Cu)
			F = Eu( A, Aod, u, X[k], Bu/auu, Bu, Cu)
			# compute estimated energy inflation ratios for each TV (k) and neighbor-agg (t) pair
			for t = 1:length(S)
				# Ctk[t,k] = Eu( A, Aod, u, X[k], X[k][S[t]], Bu, Cu )
				# Ctk_divF[t,k] = Eu( A, Aod, u, X[k], X[k][S[t]], Bu, Cu ) / F
				Ctk_temp = Eu( A, Aod, u, X[k], X[k][S[t]], Bu, Cu ) / F
				if Ctk_temp > qt[t]
					qt[t] = Ctk_temp
				end
			end
		end

		S_count = 0
		for i = 1:length(qt)
			if qt[i] <= 2.5
				S_count += 1
			end
		end
		
		# if isempty(tildeS)
		# skip second loop if it won't find anything to put in tildeS
		if S_count == 0
			# let 0 indicate NotFound
			return 0
		else 

			# second loop
			tildeS = zeros(Int,S_count)
			S_count = 0
			for i = 1:length(qt)
				if qt[i] <= 2.5
					S_count += 1
					tildeS[S_count] = S[i]
				end
			end

			return tildeS[findmin(aggregateSize[tildeS])[2]] 
		end
	end
end

# define E_u(xx,y) [equation (3.20)]
function Eu(A::SparseMatrixCSC,Aod::SparseMatrixCSC,u::Int64,xx::Vector{Float64},y::Float64,Bu::Float64,Cu::Float64)
	a_uu = Auu(A,u)
	return ((1/2*a_uu*y-Bu)*y+Cu)[1]
end

# a_{uu} [equation (3.20)]
function Auu(A::SparseMatrixCSC,u::Int64)
	return A[u,u]
end

# B_u(xx) [equation (3.20)]
function Bux(Aod::SparseMatrixCSC,u::Int64,xx::Vector{Float64})
	# rows = Aod.colptr[u]:Aod.colptr[u+1]-1
	# return At_mul_B( xx[ Aod.rowval[ rows ] ], -Aod.nzval[ rows ] )
	return At_mul_B( xx[ Aod.rowval[ Aod.colptr[u]:Aod.colptr[u+1]-1 ] ], -Aod.nzval[ Aod.colptr[u]:Aod.colptr[u+1]-1 ] )
end

# C_u(xx) [equation (3.20)]
function Cux(Aod::SparseMatrixCSC,u::Int64,xx::Vector{Float64})
	# rows = Aod.colptr[u]:Aod.colptr[u+1]-1
	# return 1/2*At_mul_B( (xx[ Aod.rowval[ rows ] ]).^2, -Aod.nzval[ rows ] )
	return 1/2*At_mul_B( (xx[ Aod.rowval[ Aod.colptr[u]:Aod.colptr[u+1]-1 ] ]).^2, -Aod.nzval[ Aod.colptr[u]:Aod.colptr[u+1]-1 ] )
end

# B_u(xx) [equation (3.20)]
# C_u(xx) [equation (3.20)]
function BuxCux(Aod::SparseMatrixCSC,u::Int64,xx::Vector{Float64})
	a = xx[ Aod.rowval[ Aod.colptr[u]:Aod.colptr[u+1]-1 ] ]
	b = -Aod.nzval[ Aod.colptr[u]:Aod.colptr[u+1]-1 ]
	return At_mul_B(a, b)[1], 1/2*At_mul_B(a.^2, b)[1]
end

function aggregate_Affinity_RS_new(A::SparseMatrixCSC,x::Vector{Float64},d::Int64,opts::options,l)
	d =3+l
	if d>10
		d = 10
	end
	println("d = ",d)
	n = size(A,1)
	Aod = A - spdiagm(diag(A))
	At = transpose(A)
# A_v = A.colptr[2:end]-A.colptr[1:end-1]-1;
# 	A_u = At.colptr[2:end]-At.colptr[1:end-1]-1;
# 	A_v = A_v+A_u;	
# 	highD = -1*ones(n)
# 	highD[A_v .>= 6*median(A_v)] = 0;


	theta = .95
	status=-1*ones(n)
	# println("status = ", length(find(x->x== -1, status)))
	X = testVectorsLAMG(A,d,opts,status)
	E = [maximum([abs(X[i][j]) for i = 1:length(X)]) for j = 1:n]
	x = zeros(n)
		C,maxC = affinityLAMG(Aod, X);
	maxC = maximum(C,2)
		for j = 1:n
			x[j] = E[j]
		end
		# x = x/norm(x)

		opts.rsv_smoother(A,x,2,opts.w)
			# X = testVectorsLAMG(A, d, opts)

			 # C = -Aod*spdiagm(x)
			 # maxC = maximum(C,2)

			 D = transpose(C)
			 maxD = maximum(D,2)

		status=-1*ones(n)
		strong_neighbors = ones(Bool,length(C.nzval))
		for i = 1:n
				for j in C.colptr[i]:(C.colptr[i+1]-1)
					if C.nzval[j] >= theta*maxC[C.rowval[j]]
						# strong_neighbors == true for strong connections
						# thus, no change necessary here 
					else 
						# strong_neighbors == false for weak connections
						strong_neighbors[j] = false
					end
				end
		end

		strong_neighbors_row = ones(Bool,length(C.nzval))
		for i = 1:n
				for j in D.colptr[i]:(D.colptr[i+1]-1)
					if D.nzval[j] >= theta*maxD[D.rowval[j]]
						# strong_neighbors == true for strong connections
						# thus, no change necessary here 
					else 
						# strong_neighbors == false for weak connections
						strong_neighbors_row[j] = false
					end
				end
		end

	# println("strong  = ",strong_neighbors)


	# U = (highD .== 0);
	# 	# println(U)

	# 	# 3:  find undecided nodes with strong neighbors
	# 	for u=sortperm(x,rev = true)
	# 		# println("u = ",u)
	# 		if U[u]
	# 			# purge zero values (weak neighbors) manually here
	# 			#   do we need intersection with A_v ???
	# 			rows_all   = C.rowval[ C.colptr[u]:C.colptr[u+1]-1 ] 
	# 			Sval   = C.nzval[ C.colptr[u]:C.colptr[u+1]-1 ] 
	# 			# println(rows_all)
	# 			strong_mask = strong_neighbors[ C.colptr[u]:C.colptr[u+1]-1 ] 
	# 			Sval = Sval[strong_mask]
	# 			# println("sval = ", Sval )

	# 			N_u = rows_all[strong_mask];
	# 			# println(N_u)
	# 			U[u] = false
	# 			s = N_u[status[N_u].==-1]
	# 			Sval = x[s]
	# 			# Sval = Sval[status[N_u].==-1]
	# 			# println("s = ", Sval )
	# 			status[u]=u
				


	# 			rows_all   = D.rowval[ D.colptr[u]:D.colptr[u+1]-1 ] 
	# 			Sval_row   = D.nzval[ D.colptr[u]:D.colptr[u+1]-1 ] 

	# 			# println(rows_all)
	# 			strong_mask_row = strong_neighbors_row[ D.colptr[u]:D.colptr[u+1]-1 ] 
	# 			Sval_row = Sval_row[strong_mask_row]
	# 			N_u = rows_all[strong_mask_row];
	# 			# println(N_u)
	# 			s = [s;N_u[status[N_u].==-1] ]
	# 			Sval = [Sval;Sval_row[status[N_u].==-1]]
	# 			# sortperm(s[x])
	# 			# println("s = ", s )
	# 			if !isempty(s)
	# 				count = 0
	# 				for i in sortperm(x[s])
	# 				# for i in sortperm(x[s])
	# 					if status[s[i]] == -1 && highD[s[i]] == -1
	# 						status[s[i]] = u
						
	# 						U[s[i]] = false
	# 						count+=1

	# 					end
	# 					if count>=4
	# 						break
	# 					end
	# 				end
	# 			end
	# 					# println(U)

	# 		end
	# 	end  

		U = (status .== -1);
		# println(U)

		# 3:  find undecided nodes with strong neighbors
		for u=sortperm(x,rev = true)
			# println("u = ",u)
			if U[u]
				# purge zero values (weak neighbors) manually here
				#   do we need intersection with A_v ???
				rows_all   = C.rowval[ C.colptr[u]:C.colptr[u+1]-1 ] 
				Sval   = C.nzval[ C.colptr[u]:C.colptr[u+1]-1 ] 
				# println(rows_all)
				strong_mask = strong_neighbors[ C.colptr[u]:C.colptr[u+1]-1 ] 
				Sval = Sval[strong_mask]
				# println("sval = ", Sval )

				N_u = rows_all[strong_mask];
				# println(N_u)
				U[u] = false
				s = N_u[status[N_u].==-1]
				Sval = x[s]
				# Sval = Sval[status[N_u].==-1]
				# println("s = ", Sval )
				status[u]=u
				


				rows_all   = D.rowval[ D.colptr[u]:D.colptr[u+1]-1 ] 
				Sval_row   = D.nzval[ D.colptr[u]:D.colptr[u+1]-1 ] 

				# println(rows_all)
				strong_mask_row = strong_neighbors_row[ D.colptr[u]:D.colptr[u+1]-1 ] 
				Sval_row = Sval_row[strong_mask_row]
				N_u = rows_all[strong_mask_row];
				# println(N_u)
				s = [s;N_u[status[N_u].==-1] ]
				Sval = [Sval;Sval_row[status[N_u].==-1]]
				# sortperm(s[x])
				# println("s = ", s )
				if !isempty(s)
					count = 0
					for i in sortperm(x[s])
					# for i in sortperm(x[s])
						if status[s[i]] == -1
							status[s[i]] = u
						
							U[s[i]] = false
							count+=1

						end
						if count>=1
							break
						end
					end
				end
						# println(U)

			end
		end  



	return status, x
end


function aggregate_Affinity_RS_new_rnk(A::SparseMatrixCSC,x::Vector{Float64},d::Int64,opts::options,l, rnk)
	d =3+l
	if d>10
		d = 10
	end
	println("d = ",d)
	n = size(A,1)
	Aod = A - spdiagm(diag(A))
	At = transpose(A)
	A_v = A.colptr[2:end]-A.colptr[1:end-1]-1;
	A_u = At.colptr[2:end]-At.colptr[1:end-1]-1;
	A_v = A_v+A_u;	
	highD = -1*ones(n)
	highD[A_v .>= 6*median(A_v)] = 0;


	
		U = (status .== -1);
		# println(U)

		# 3:  find undecided nodes with strong neighbors
		for u=sortperm(rnk,rev = true)
			# println("u = ",u)
			if U[u]
				# purge zero values (weak neighbors) manually here
				#   do we need intersection with A_v ???
				rows_all   = C.rowval[ C.colptr[u]:C.colptr[u+1]-1 ] 
				Sval   = C.nzval[ C.colptr[u]:C.colptr[u+1]-1 ] 
				# println(rows_all)
				strong_mask = strong_neighbors[ C.colptr[u]:C.colptr[u+1]-1 ] 
				Sval = Sval[strong_mask]
				# println("sval = ", Sval )

				N_u = rows_all[strong_mask];
				# println(N_u)
				U[u] = false
				s = N_u[status[N_u].==-1]
				Sval = x[s]
				# Sval = Sval[status[N_u].==-1]
				# println("s = ", Sval )
				status[u]=u
				


				rows_all   = D.rowval[ D.colptr[u]:D.colptr[u+1]-1 ] 
				Sval_row   = D.nzval[ D.colptr[u]:D.colptr[u+1]-1 ] 

				# println(rows_all)
				strong_mask_row = strong_neighbors_row[ D.colptr[u]:D.colptr[u+1]-1 ] 
				Sval_row = Sval_row[strong_mask_row]
				N_u = rows_all[strong_mask_row];
				# println(N_u)
				s = [s;N_u[status[N_u].==-1] ]
				Sval = [Sval;Sval_row[status[N_u].==-1]]
				# sortperm(s[x])
				# println("s = ", s )
				if !isempty(s)
					count = 0
					for i in sortperm(rnk[s])
					# for i in sortperm(x[s])
						if status[s[i]] == -1
							status[s[i]] = u
						
							U[s[i]] = false
							count+=1

						end
						if count>=4
							break
						end
					end
				end
						# println(U)

			end
		end  



	return status, x
end



function aggregate_degree_aware(A::SparseMatrixCSC,x::Vector{Float64},d::Int64,opts::options,l)
	d =2+l
	# println("d = ",d)
	n = size(A,1)
	Aod = A - spdiagm(diag(A))

	Aod.nzval = ones(length(Aod.nzval))
	dout = findDegrees(A)
    din = findDegrees(A')
    B = speye(n) - Aod*spdiagm(1./dout)


	theta = .9
	status=-1*ones(n)
	# # println("status = ", length(find(x->x== -1, status)))
	X = testVectorsLAMG(A,d,opts,status)
	E = [maximum([abs(X[i][j]) for i = 1:length(X)]) for j = 1:n]
	x = zeros(n)
		C,maxC = affinityLAMG(Aod, X);
		for j = 1:n
			x[j] = E[j]
		end
		x = x/norm(x,1)

		opts.rsv_smoother(A,x,8,opts.w)
	# 		# X = testVectorsLAMG(A, d, opts)

	# 		 # C = -Aod*spdiagm(x)
	# 		 maxC = maximum(C,2)

			 D = transpose(C)
		# C = abs(Aod)

		# D = transpose(abs(Aod))
	maxC = maximum(C,2)
			 maxD = maximum(D,2)


		status=-1*ones(n)
		strong_neighbors = ones(Bool,length(C.nzval))
		for i = 1:n
				for j in C.colptr[i]:(C.colptr[i+1]-1)
					if C.nzval[j] >= theta*maxC[C.rowval[j]]
						# strong_neighbors == true for strong connections
						# thus, no change necessary here 
					else 
						# strong_neighbors == false for weak connections
						strong_neighbors[j] = false
					end
				end
		end

		strong_neighbors_row = ones(Bool,length(C.nzval))
		for i = 1:n
				for j in D.colptr[i]:(D.colptr[i+1]-1)
					if D.nzval[j] >= theta*maxD[D.rowval[j]]
						# strong_neighbors == true for strong connections
						# thus, no change necessary here 
					else 
						# strong_neighbors == false for weak connections
						strong_neighbors_row[j] = false
					end
				end
		end

	# println("strong  = ",strong_neighbors)

	
    dt = dout+din
		U = (status .== -1);
		# println(U)

		# 3:  find undecided nodes with strong neighbors
		for u=sortperm(dt,rev = true)
			# println("u = ",u)
			if U[u]
				# purge zero values (weak neighbors) manually here
				#   do we need intersection with A_v ???
				rows_all   = C.rowval[ C.colptr[u]:C.colptr[u+1]-1 ] 
				Sval   = C.nzval[ C.colptr[u]:C.colptr[u+1]-1 ] 
				# println(rows_all)
				strong_mask = strong_neighbors[ C.colptr[u]:C.colptr[u+1]-1 ] 
				Sval = Sval[strong_mask]
				# println("sval = ", Sval )

				N_u = rows_all[strong_mask];
				# println(N_u)
				U[u] = false
				s = N_u[status[N_u].==-1]
				# Sval = x[s]
				Sval = Sval[status[N_u].==-1]
				# println("s = ", Sval )
				status[u]=u
				


				rows_all   = D.rowval[ D.colptr[u]:D.colptr[u+1]-1 ] 
				Sval_row   = D.nzval[ D.colptr[u]:D.colptr[u+1]-1 ] 

				# println(rows_all)
				strong_mask_row = strong_neighbors_row[ D.colptr[u]:D.colptr[u+1]-1 ] 
				Sval_row = Sval_row[strong_mask_row]
				N_u = rows_all[strong_mask_row];
				# println(N_u)
				s = [s;N_u[status[N_u].==-1]]
				Sval = [Sval;Sval_row[[status[N_u].==-1]]]
				# sortperm(s[x])
				# println("s = ", s )
				if !isempty(s)
					count = 0
					for i in sortperm(Sval, rev = true)
					# for i in sortperm(x[s])
						if status[s[i]] == -1
						status[s[i]] = u
					
						U[s[i]] = false
												count+=1

						end
						if count>=8
							break
						end
					end
				end
						# println(U)

			end
		end  



	return status, x
end

function aggregation_RS_new(A::SparseMatrixCSC, x::Vector{Float64},theta, mix::Int64)
	n = size(A,1)
	println("size = ",n)

	# status=-1*ones(n)
	# A_v = A.colptr[2:end]-A.colptr[1:end-1]-1;
	# status[A_v .>= 8*median(A_v)] = 0
	# println("median(A_v)", median(A_v))
	# println("max(A_v)", maximum(A_v))
	# println("number of high = ", n - countnz(status))
	# Aod = -A+spdiagm(diag(A))
	 C = -A*spdiagm(x)
	status=-1*ones(n)
	strong_neighbors = ones(Bool,length(C.nzval))
	maxC = maximum(C, 2);
	for i = 1:n
			for j in C.colptr[i]:(C.colptr[i+1]-1)
				if C.nzval[j] >= theta*maxC[C.rowval[j]]
					# strong_neighbors == true for strong connections
					# thus, no change necessary here 
				else 
					# strong_neighbors == false for weak connections
					strong_neighbors[j] = false
				end
			end
	end

	# println("strong  = ",strong_neighbors)

		U = (status .== -1);
		# println(U)

		# 3:  find undecided nodes with strong neighbors
		for u=sortperm(x,rev = true)
			# println("u = ",u)
			if U[u]
				# purge zero values (weak neighbors) manually here
				#   do we need intersection with A_v ???
				rows_all   = C.rowval[ C.colptr[u]:C.colptr[u+1]-1 ] 
				# println(rows_all)
				strong_mask = strong_neighbors[ C.colptr[u]:C.colptr[u+1]-1 ] 

#
				N_u = rows_all[strong_mask];
				# println(N_u)
				U[u] = false
				s = N_u[status[N_u].==-1]
				status[u]=u
				if !isempty(s)
				status[N_u[status[N_u].==-1]] = u
				

				U[N_u] = false
				end
						# println(U)

			end
		end  



	return status
end

function aggregation_RS_CR_new(A::SparseMatrixCSC, x::Vector{Float64},theta, mix::Int64,E)
	n = size(A,1)
	println("size = ",n)

	# status=-1*ones(n)
	# A_v = A.colptr[2:end]-A.colptr[1:end-1]-1;
	# status[A_v .>= 8*median(A_v)] = 0
	# println("median(A_v)", median(A_v))
	# println("max(A_v)", maximum(A_v))
	# println("number of high = ", n - countnz(status))
	# Aod = -A+spdiagm(diag(A))
	 C = -A*spdiagm(x)
	status=-1*ones(n)
	strong_neighbors = ones(Bool,length(C.nzval))
	maxC = maximum(C, 2);
	for i = 1:n
			for j in C.colptr[i]:(C.colptr[i+1]-1)
				if C.nzval[j] >= theta*maxC[C.rowval[j]]
					# strong_neighbors == true for strong connections
					# thus, no change necessary here 
				else 
					# strong_neighbors == false for weak connections
					strong_neighbors[j] = false
				end
			end
	end

	# println("strong  = ",strong_neighbors)

		U = (status .== -1);
		# println(U)

		# 3:  find undecided nodes with strong neighbors
		for u=sortperm(E)
			# println("u = ",u)
			if U[u]
				# purge zero values (weak neighbors) manually here
				#   do we need intersection with A_v ???
				rows_all   = C.rowval[ C.colptr[u]:C.colptr[u+1]-1 ] 
				# println(rows_all)
				strong_mask = strong_neighbors[ C.colptr[u]:C.colptr[u+1]-1 ] 

#
				N_u = rows_all[strong_mask];
				# println(N_u)
				U[u] = false
				s = N_u[status[N_u].==-1]
				status[u]=u
				if !isempty(s)
				status[N_u[status[N_u].==-1]] = u
				

				U[N_u] = false
				end
						# println(U)

			end
		end  



	return status
end






function aggregation_DNR(A::SparseMatrixCSC,mix::Int64)
	n = size(A,1)
	d = diag(A)

	status=zeros(n); 

	xold = copy(d)

	# A_v = A.colptr[2:end]-A.colptr[1:end-1]-1;

	count = 0
	while count<length(d) 
		# g = countnz(status)
		# println(countnz(xold) )

	# maxC = maximum(C, 1);
	i = indmax(xold)
	aggregateSize=0;

	xold[i] = 0
   	status[i]=i
   	count+=1
	# println(xold)
	s  = [i]
   		for j in collect(A.colptr[i]:(A.colptr[i+1]-1))
   			if status[A.rowval[j]] ==0
	   				xold[A.rowval[j]] = 0
					status[A.rowval[j]] = i
					count+=1
					aggregateSize+=1; 

					# S.nzval[find(x->x ==S.rowval[j], S.rowval)] = 0
					# S.nzval[j] = 0
					push!(s,A.rowval[j])

   			end
   		end

   		if mix ==2
	   		for h in s
		   		for j in collect(A.colptr[h]:(A.colptr[h+1]-1))
		   			if status[A.rowval[j]] ==0
			   				xold[A.rowval[j]] = 0
							status[A.rowval[j]] = i
							count+=1
							aggregateSize+=1; 

							# S.nzval[find(x->x ==S.rowval[j], S.rowval)] = 0
							# S.nzval[j] = 0
		   			end
		   		end
	   		end
	   	end


	   	if mix ==3
	   		if aggregateSize<=6
		   		for h in s
			   		for j in collect(A.colptr[h]:(A.colptr[h+1]-1))
			   			if status[A.rowval[j]] ==0
				   				xold[A.rowval[j]] = 0
								status[A.rowval[j]] = i
								count+=1
								aggregateSize+=1; 

								# S.nzval[find(x->x ==S.rowval[j], S.rowval)] = 0
								# S.nzval[j] = 0
			   			end
			   		end
		   		end
	   		end
	   	end


   	



	end
		
	

	return status
end


############## CR CODE ##############################



function v7_aggregate_CR(A::SparseMatrixCSC, theta::Float64, mix::Int64,l::Int64,opts::options)
	n = size(A,1)
	println("size = ",n)
	println()

	
	a_max = .7

	# save an off diagonal only version of A 
	# TODO:  do this more elegantly
	Aod = A - spdiagm(diag(A))
	# run aggregation routine
	At = transpose(A)
	status,rnk =aggregate(A,At, Aod, a_max, opts.I_max,l,opts);
	

	status[status .== 0]= (1:length(status))[status .== 0];
		# Nate's addition:  set unaggregated nodes to their own aggregate
	status[status .== -1] = (1:length(status))[status .== -1];


	return status,rnk
end



