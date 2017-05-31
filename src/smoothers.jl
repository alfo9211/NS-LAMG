function richardson!(A::SparseMatrixCSC,x::Array{Float64,1},b::Array{Float64,1},rnk::Array{Float64,1},numIters::Int64,w::Float64)
  numRows, numCols = size(A); 

  if (numRows != numCols)
    error("The matrix A is not square")
  end

  if (numRows != size(b,1))
    error("The right hand side is inconsisitent with the size of the matrix")
  end

  if (numCols != size(x,1))
    error("The initial guess is inconsistent with the size of the matrix")
  end

for k=1:numIters 


    xhat = b-A*x; 
 

    for i=1:numCols
      x[i] = x[i] + w*xhat[i]; 
    end 

end

  y = copy(x)
  alpha = dot(y,vec(rnk))/dot(vec(rnk),vec(rnk)) 
      
    for i=1:length(x)
      x[i] = y[i]-alpha[1]*rnk[i]; 
    end
end

function gaussseidel!(A::SparseMatrixCSC,x::Vector, b::Vector,numIters::Int64,w::Float64)

  numRows,numCols = size(A); 
  
  if (numRows != numCols)
    error("The matrix A is not square")
  end

  if (numRows != size(b,1))
    error("The right hand side is inconsisitent with the size of the matrix")
  end

  if (numCols != size(x,1))
    error("The initial guess is inconsistent with the size of the matrix")
  end
  y = copy(x)

  for kk=1:numIters 

    dElem = 0.0;
    for jj in 1:size(A,1);  
      xhat = b[jj];
      for ii = A.colptr[jj]:A.colptr[jj+1]-1
        row = A.rowval[ii]; 
        val = A.nzval[ii]; 
        if (row == jj) 
          dElem = val; 
        else 
          xhat -= val*x[row]; 
        end 
      end 
      y[jj] = xhat/dElem; 
    end 

    for i=1:numCols
      x[i] = (1-w)*x[i] + w*y[i]; 
    end 


  end 

end 

# --------------------------------------------------------------
# Function that performs numIters iterations of weighted jacobi on 
# Ax = b with x as the initial guess.  
# Note: not garnateed to converge for nonsymmetric matrices.
# -------------------------------------------------------------- 
function wjacobib!(A::SparseMatrixCSC,x::Vector,b::Vector,rnk::Vector,numIters::Int64,w::Float64)
  numRows, numCols = size(A); 

  if (numRows != numCols)
    error("The matrix A is not square")
  end

  if (numRows != size(b,1))
    error("The right hand side is inconsisitent with the size of the matrix")
  end

  if (numCols != size(x,1))
    error("The initial guess is inconsistent with the size of the matrix")
  end

for k=1:numIters 

    xhat = copy(b); 

    d= zeros(numRows); 
      
    for j = 1:numCols 
      for i = A.colptr[j]:A.colptr[j+1]-1
        row = A.rowval[i]; 
        val = A.nzval[i];
        if (row != j)
          xhat[row] -= val*x[j]; 
        else
          d[row] = val; 
        end
      end 
    end

    for i = 1:numCols 
      # if d[i]>10^(-4.0)
        xhat[i] = xhat[i]/d[i];  
      # end
    end 

    for i=1:numCols
      x[i] = (1-w)*x[i] + w*xhat[i]; 
    end 

  end 

  # y = copy(x)
  # alpha = dot(y,vec(rnk))/dot(vec(rnk),vec(rnk)) 
      
  #   for i=1:length(x)
  #     x[i] = y[i]-alpha[1]*rnk[i]; 
  #   end
end


# --------------------------------------------------------------
# Function that performs numIters iterations of weighted jacobi on 
# Ax = 0 with x as the initial guess. It normalizes the current 
# guess x so that it converges for nonsymmetric matrics. 
# -------------------------------------------------------------- 
function wjacobi!(A::SparseMatrixCSC,x::Vector,numIters::Int64,w::Float64)
  numRows, numCols = size(A); 

  if (numRows != numCols)
    error("The matrix A is not square")
  end


  if (numCols != size(x,1))
    error("The initial guess is inconsistent with the size of the matrix")
  end

for k=1:numIters 

    xhat = zeros(numRows); 

    d= zeros(numRows); 
      
    for j = 1:numCols 
      for i = A.colptr[j]:A.colptr[j+1]-1
        row = A.rowval[i]; 
        val = A.nzval[i];
        if (row != j)
          xhat[row] -= val*x[j]; 
        else
          if (val ==0.0)
            println("VAL IS 0")
          end
          d[row] = val; 
        end
      end 
    end




    for i = 1:numCols 
      xhat[i] = xhat[i]/d[i];  
    end 

    for i=1:numCols
      x[i] = (1-w)*x[i] + w*xhat[i]; 
    end 

  end 
  y = copy(x)
  normy= norm(y,1)
    # println("normy = ", y)

    for i=1:length(x)
      x[i] = y[i]/normy; 
    end

end


# --------------------------------------------------------------
# Function that performs number of iterations (its) of weighted Gauss Seidel
# on Ax = b with x as the initial guess.  
# Note: not garnateed to converge for nonsymmetric matrices.
# -------------------------------------------------------------- 
function wgausseidel!(A::SparseMatrixCSC,x::Vector,its::Int64,w::Float64)
  DL = tril(A)
  for i = 1:its
    y = DL\(A*x)
    for i=1:length(x)
      y[i] =x[i]-w*y[i]; 
    end 
    normy= norm(y,1)
    for i=1:length(x)
      x[i] = y[i]/normy; 
    end 

  end

end



function usymqr!(A::SparseMatrixCSC, x::Vector,b::Vector, its::Int64 )
	r = b-A*x
	beta1 = norm(r,2)
	beta = beta1
	gamma = norm(r,2)
	p = [zeros(length(r)) r./beta]
	q = [zeros(length(r)) r./gamma]
 	for i = 2:its+1
 		u = A*q[:,i]-gamma*p[:,i-1]
 		v = A'*p[:,i]-beta*p[:,i-1]
 		alpha = (p[:,i]'*u)[1]
 		# println(alpha*p[:,i])
 		u = u - alpha*p[:,i]
 		v = v -alpha*q[:,i]
 		beta = norm(u,2)
 		gamma = norm(v,2)
 		if beta == 0 ||gamma == 0
 			return 
 		else
 			p = [p u./beta]
 			q = [q v./gamma]
 		end
 	end

 		T = p[:,2:end]'*A*q[:,2:end]

		e = zeros(size(T,1))
		e[1] = 1
		S = vcat(T,beta*e')
		h = S\(beta1*(vcat(e,0)))
		y = q[:,2:end]*h

		for i= 1:length(x)
			x[i]+=y[i]
		end
end


