# operation counter sub-type for level structure
type OPCOUNT
	# countMatVecs for L
	L::Int64
	# countMatVecs for P
	P::Int64
	# countMatVecs for L
	R::Int64

	OPCOUNT(L,P,R) = new(copy(L),copy(P),copy(R))
end 


type level
	A::SparseMatrixCSC
	P::SparseMatrixCSC
	R::SparseMatrixCSC
	Rinj::SparseMatrixCSC
	Q::SparseMatrixCSC
	ZFC::Vector
	rnk::Vector
	x::Vector
	b::Vector
	histx::SparseMatrixCSC
	histAx::SparseMatrixCSC
	histcount::Int64
	opcount::OPCOUNT
	lvFunc::Function
	level(A,P,R,Rinj,Q,ZFC,rnk,x,b, histx,histAx,histcount,opcount,lvFunc) = new(copy(A),copy(P),copy(R),copy(Rinj),copy(Q),copy(ZFC),copy(rnk),copy(x),copy(b),copy(histx), copy(histAx), copy(histcount), opcount,lvFunc)
end 



type tentitive
	T::SparseMatrixCSC
	
	tentitive(T) = new(copy(T))
end 