function findDegrees(Ac::SparseMatrixCSC)
  degrees = zeros(Int,length(Ac.colptr)-1)
    for i = 1:length(degrees)
      # minus 1 to remove count for diagonal entry
      #  Note:  degrees[i] == -1  implies i has degree zero
      degrees[i] = Ac.colptr[i+1]-Ac.colptr[i]-1
    end
  return degrees
end 

function elimination(A::SparseMatrixCSC,At::SparseMatrixCSC,opts::options)
  # 1:  initialize
  n = size(A,1)
  # TODO:  DO WE REALLY NEED TO COPY A HERE?
  # Ac = copy(A)
  Ac = A 
  Atc = At
  nc_prev = n 
  P = speye(n)
    R = speye(n)

  Q = spzeros(n,n)
  q = 0
  Rinj = speye(n)
  # Pii = {}
  # Qii = {}

  # 2:  loop-de-doop
  # while nc_prev > 1
    
    # store Z,F,C identifiers in one vector
    #  ZFC[i] = 0  ==>  i \in Z
    #         = 1  ==>  i \in F
    #         = 2  ==>  i \in C
    # start all nodes in C
    ZFC = 2*ones(nc_prev)

    # 3:  find disconnected and low degree nodes
    # note:  this will return degrees = -1 for degree zero nodes
    # degrees = Ac.colptr[2:end]-Ac.colptr[1:end-1]-1
    dout = findDegrees(Ac)
    din = findDegrees(Atc)

    # println("find Degree = ",dout)
    #  println("find Degree = ",din)

    # node set Z is ZFC .== 0
    # ZFC[ dout .== -1 ] = 0
    #     ZFC[ din .== -1 ] = 0

    nZ = sum(ZFC .== 0)
    # node set F is ZFC .== 1
    ZFC[ lowDegreeNodes(Ac,Atc,1,dout,din) ] = 1

    # println("out of lowdegree")
    nF = sum(ZFC .== 1)
    nC = nc_prev - nZ - nF

    
    # 4:  terminate pruning if no disconnected commponents and few low degree nodes
    # if (sum(ZFC .== 0) == 0 ) && (sum(ZFC .== 1 )<0.01*nc_prev)
    if (nC == 0 ) ||  (nF<opts.minEliminationFraction*nc_prev) 

      # || nC == nF
    # if (nC == 0 ) || nC == nF
    # println("nC = ", nC)
    # println("nF = ", nF)

    	println(opts.minEliminationFraction*nc_prev)
      # 5:
      # return Ac,P,Q,Pii,Qii

      ZFC = vec(sum(Rinj,1))
      return Ac,P,R,Q,Rinj,ZFC
      # 6:
    end
    
    # 7:  coarse points
    # node set C is ZFC .== 2
    q += 1

    # 8:  build permutation order for [Z,F,C]
    # this strategy appears to check out:  transpose(permmat)*ZFC should be sorted in ascending order
    perm = sortperm(ZFC)
    # println("perm = ",perm)
    # permmat = speye(length(perm))
    # permmat.rowval = perm

    # permmatT = transpose(permmat)

    # 9:  Build Qi, and concatenate onto Qii (skipping this)
    #  appears to be a typo in paper:  third term along diagonal of Qi should be O_C
    #  also need another transpose(permutation) on RHS of Qi
    AFFinv = diag(Ac[ZFC .== 1,ZFC .== 1]).^(-1) 

    # Qi = A_mul_Bt( permmat*spdiagm( [ zeros( nZ ); AFFinv; zeros( nC ) ] ), permmat)
 Qi= permute!(spdiagm( [ zeros( nZ ); AFFinv; zeros( nC ) ]),perm,perm)

    # push!(Qii,Qi)
    # Ri = A_mul_Bt([ spzeros(nC,nZ) -Ac[ZFC .== 2, ZFC .== 1]*spdiagm(AFFinv) speye(nC) ],permmat)
    println("size mat = ", size([ spzeros(nC,nZ) -Ac[ZFC .== 2, ZFC .== 1]*spdiagm(AFFinv) speye(nC) ]))
    println("length perm = ", length(perm))
    println("nF = ", nF)
     println("nC = ", nC)

     Ri= permute([ spzeros(nC,nZ) -Ac[ZFC .== 2, ZFC .== 1]*spdiagm(AFFinv) speye(nC) ],round(Int64,[i for i in 1:nC]),perm)
     println("size mat = ", size([spzeros(nC,nZ+nF) speye(nC)]))
    println("length perm = ", length(perm))
    println("nF = ", nF)
     println("nC = ", nC)
    # Rinj = A_mul_Bt( [spzeros(nC,nZ+nF) speye(nC)], permmat)*Rinj
     Rinj= permute( [spzeros(nC,nZ+nF) speye(nC)],round(Int64,[i for i in 1:nC]),perm)*Rinj
println("size mat = ", size( [ spzeros(nZ,nC );-spdiagm(AFFinv)*Ac[ZFC .== 1, ZFC.== 2]; speye(nC) ]))
    println("length perm = ", length(perm))
    println("nF = ", nF)
     println("nC = ", nC)
  	# Pi = permmat*[ spzeros(nZ,nC );-spdiagm(AFFinv)*Ac[ZFC .== 1, ZFC.== 2]; speye(nC) ]
    Pi= permute( [ spzeros(nZ,nC );-spdiagm(AFFinv)*Ac[ZFC .== 1, ZFC.== 2]; speye(nC) ],perm,round(Int64,[i for i in 1:nC]))
    # R = A_mul_Bt( [spzeros(nC,nZ+nF) speye(nC)], permmat)


    # 11:  store new Q 
    #    typo in pseudocode here:  no need for R
    #    should be:  Q <- Q + P*Qi*transpose(P)
    # println("dim of Qi",  size(Qi))
    # println("dim of Ri" , size(Ri))
    # # println("dim of Rinj" , size(Qi*Rinj))
    # println("dim of Pi", size(Pi))
    #     println("dim of P", size(P))
    #             println("dim of R", size(R))


        # Rinj = A_mul_Bt( [spzeros(nC,nZ+nF) speye(nC)], permmat)*Rinj

    Q = Q + P*Qi*R



    # println("sizeQ ", size(Qi))


    #   store new P 
    #   typo in pseudocode here:  multiplication order is wrong
    #   should be:  P <- P*Pi
    P = P*Pi
    R = Ri*R
    #   store coarse grid Laplacian Ac
    #   typo in pseudocode
    #   should be:   Ac <- Pi'*Ac*Pi 
    # Ac = At_mul_B(Pi, Ac*Pi)
    Ac = Ri*Ac*Pi
    Atc =Ac'

    # println(Ac)
    nc_prev = size(Ac,1)
    # println("nc_prev = ", nc_prev)

    # 12:
  # end
  # if we've coarsened to exactly 1 node...
  ZFC = vec(sum(Rinj,1))
  return Ac,P,R,Q,Rinj,ZFC
end



function lowDegreeNodes(A::SparseMatrixCSC,At::SparseMatrixCSC,d::Int64,dout::Vector,din::Vector)
  # 1: find low degree nodes
  n = size(A,1)
  println("size n = ", n)
  # degrees = A.colptr[2:end]-A.colptr[1:end-1] -1

  # U[u] == true  ==> u is low degree
  #hard coding for now the degrees
  # U = min(dout.>= 1 , dout .<= 1 )
  U = collect(dout.==1)
   V = collect(din.==1)
  # V = min(din .>= 1 , din .<= 1 )
 Z = min((dout+din) .>= 1 , (dout+din) .<= 8 )
  # println("U = ", length(find(x->x == true , U)))
  # println("V = ",length(find(x->x == true , U)))

# println(degrees)
  # println(U)
  # 2:  visited = 0  ==> NotVisited
  #             = 1  ==> FNode
  #             = 2  ==> NotEliminated
  visited = zeros(length(U))
  # 3:  loop-de-doop
  for u = 1:n 
  	# println("u = ", u)

  	# println("visited = ", visited)
    if Z[u] 
      # 4:  if u has not been visited ...
      if visited[u] == 0
        # 5:  check for Fnodes in neighborhood
        #     if none, u may be eliminated
        Au = A.rowval[ A.colptr[u]:A.colptr[u+1]-1 ]
        Av = At.rowval[ At.colptr[u]:At.colptr[u+1]-1 ]

        # println("Au = ", Au)
        # println("visited[Au] = ", visited[Au])

        # if sum( visited[Au] .== 1 ) == 0
        #   # 7:  set neighbors of u as NotEliminated
        #   visited[Au] = 2
        #   # 6:  set u as Fnode
        #   visited[u] = 1
        #   # 8:  else u has an F-neighbor
        # else 
        #   # 9:  do not eliminate u
        #   visited[u] = 2
        #   # 10:
        # end
        if 1 in visited[Au] || 1 in visited[Av] 
          # 8:  then u has an F-neighbor
          # 9:  do not eliminate u
          visited[u] = 2
        else 
          # else u does not have an F-neighbor
          # 7:  set neighbors of u as NotEliminated
          visited[Au] = 2
          # 6:  set u as Fnode
          visited[u] = 1
          # 10:
        end
        # 11:
      end
    end



    if V[u] 
      # 4:  if u has not been visited ...
      if visited[u] == 0
        # 5:  check for Fnodes in neighborhood
        #     if none, u may be eliminated
                Au = A.rowval[ A.colptr[u]:A.colptr[u+1]-1 ]

        Av = At.rowval[ At.colptr[u]:At.colptr[u+1]-1 ]
        # println("Av = ", Av)
        # println("visited[Av] = ", visited[Av])

        # if sum( visited[Au] .== 1 ) == 0
        #   # 7:  set neighbors of u as NotEliminated
        #   visited[Au] = 2
        #   # 6:  set u as Fnode
        #   visited[u] = 1
        #   # 8:  else u has an F-neighbor
        # else 
        #   # 9:  do not eliminate u
        #   visited[u] = 2
        #   # 10:
        # end
        if 1 in visited[Au] || 1 in visited[Av] 
          # 8:  then u has an F-neighbor
          # 9:  do not eliminate u
          visited[u] = 2
        else 
          # else u does not have an F-neighbor
          # 7:  set neighbors of u as NotEliminated
          visited[Av] = 2
          # 6:  set u as Fnode
          visited[u] = 1
          # 10:
        end
        # 11:
      end
    end

      if U[u] 
      # 4:  if u has not been visited ...
      if visited[u] == 0
        # 5:  check for Fnodes in neighborhood
        #     if none, u may be eliminated
                Au = A.rowval[ A.colptr[u]:A.colptr[u+1]-1 ]

        Av = At.rowval[ At.colptr[u]:At.colptr[u+1]-1 ]
        # println("Av = ", Av)
        # println("visited[Av] = ", visited[Av])

        # if sum( visited[Au] .== 1 ) == 0
        #   # 7:  set neighbors of u as NotEliminated
        #   visited[Au] = 2
        #   # 6:  set u as Fnode
        #   visited[u] = 1
        #   # 8:  else u has an F-neighbor
        # else 
        #   # 9:  do not eliminate u
        #   visited[u] = 2
        #   # 10:
        # end
        if 1 in visited[Au] || 1 in visited[Av] 
          # 8:  then u has an F-neighbor
          # 9:  do not eliminate u
          visited[u] = 2
        else 
          # else u does not have an F-neighbor
          # 7:  set neighbors of u as NotEliminated
          visited[Av] = 2
          # 6:  set u as Fnode
          visited[u] = 1
          # 10:
        end
        # 11:
      end
    end
  #   # 12:
  end
  # 13:  return Fnodes
  return visited .== 1
end