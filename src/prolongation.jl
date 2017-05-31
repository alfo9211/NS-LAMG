include("../src/aggregation.jl");
include("../src/smoothers.jl");

# --------------------------------------------------------------
# Function finds the Prolongation operator P from aggregation vector. 
# --------------------------------------------------------------
function prolongation(status::Vector)
	n = length(status)

	# build P directly using sparse(I,J,V) with J = status
	 P = sparse(sparse(collect(1:n), vec(round(Int64,status)), ones(n) ))
	 println("Size P = ",size(P))
	# need to push out all-zero columns of P
	# P.colptr = unique(P.colptr)
	# P.n = length(P.colptr) - 1

	return SparseMatrixCSC(n,length(unique(status)),unique(P.colptr),P.rowval,P.nzval)
end





function createPrologationOpt_DNR(A::AbstractMatrix, opts::options, gamma_l::Float64, d::Int64)
	
	# run LAMG aggregation routine
	#  status indicates which vertices should be aggregated:  if status[i] == status[j], then aggregate i and j together
	status = aggregate_DNR(A, opts, gamma_l, d)
	n = length(status)

	# build P directly using sparse(I,J,V) with J = status
	P = sparse([1:n], int64(status), ones(n) )
	# need to push out all-zero columns of P
	P.colptr = unique(P.colptr)
	P.n = length(P.colptr) - 1

	return P 

end