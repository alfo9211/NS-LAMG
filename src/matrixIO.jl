# using MatrixMarket

# use built-in MatrixMarket Read function
# requires adding MatrixMarket package to Julia
function MatrixMarketRead(filename)
	 A = readdlm(filename)
	 A = A[4:end,1:3]
	return sparse( round(Int64,A[:,1]), round(Int64,A[:,2]), convert(Array{Float64},A[:,3]))
end

# custom MatrixMarket Write function
function MatrixMarketWrite(A::SparseMatrixCSC, filename, comment="")
	# open desired file to write-to
	f_h = open(filename, "w")
	# write matrix market headers
	write(f_h, string("%%MatrixMarket matrix coordinate real general","\n") )
	write(f_h, string("%",comment,"\n") )
	# find the nonzero values in A, and their indices
	r,c,v = findnz(A)
	# write size and nnz's
	(nx, ny) = size(A)
	number_vals = nnz(A)
	write(f_h, "$nx $ny $number_vals\n")
	# write row,col,vals
	for i in 1:length(r)
		write(f_h, "$(r[i]) $(c[i]) $(v[i])\n")
	end
	# close file
	close(f_h)
end

