include("../src/matrixIO.jl")
include("../src/options.jl"); 
include("../src/hierarchy.jl")
include("../src/solve.jl"); 
# include("MarkovChainExamples.jl")

# using IterativeSolvers

dir  = "directedMTX/"
dir_results = "results/"

# ln = "amazon0601.mtx"
ln = "advogato.mtx"
# ln = "citeseer.mtx"
# ln = "cit-HepTh.mtx"

# ln = "web-Stanford.mtx"
# ln = "wikipedia_link_ja.mtx"
ln = "cfinder-google.mtx"
# ln = "moreno_cattle_cattle.mtx"
# ln = "facebook-wosn-wall.mtx"
# ln ="cit-HepPh.mtx"
# ln = "link-dynamic-nlwiki.mtx"
# ln = "dbpedia-all.mtx"
# ln = "uniformchain1.mtx"
# MTXFiles = open(string(dir,"MtxFiles.txt"))



flag = 0
# for ln in eachline(MTXFiles)

	println("Testing NS-LAMG on graph: ", string(strip(ln)))
	# Upload Graph and insure it is a normalized graph Laplacain, i.e L = I -AD^{-1}
	A = MatrixMarketRead(string(dir,strip(ln)));
	A = A- spdiagm(diag(A));
	n = size(A,1);
	A =  speye(n)-A*spdiagm(1./vec(sum(A,1)));


	# Run test twice for the first graph as Julia has a just in time complier
	if flag ==0
		i = 2
		flag = 1
	else
		i =1
	end

	for l = 1:i


			rnk = one(n);
		
			#Upload Solver Options and Solver Statistics	
			opts = setopts();
			stats = setSolveStats();


			srand(1234)
			xexact = rand(n);
			xexact = xexact/norm(xexact);
			#Ensure right-hand side is in the Range of L 
			b = A*xexact;
			x = rand(n);


			# for aggFunc in {agg_CR,agg_rsv}
				aggFunc = agg_CR
				stats.conFacHistory = Array(Float64,0)
				stats.setupT = 0.0
				stats.solveT = 0.0
				stats.converged = false
				
				srand(1234)
				rnk = rand(n)
				time_0 = time()
				lv = hierarchy_test(A,A',vec(rnk),x,b,opts,aggFunc)
				stats.setupT = (time() - time_0)

				time_0 = time();
				lv = find_rsv_RS_unweighted!(lv, lv[1].rnk, opts,stats)
				stats.rsvT = (time() - time_0);
				time_0 = time()
				solveV!(lv,stats,opts)
				stats.solveT = (time() - time_0)

				PrintStatistics(string(dir_results,string(aggFunc),"/",ln[1:search(ln,".mtx")[1]-1],".txt"),lv,stats,opts)

# 
				lv = 0; 


# 				# aggFunc = agg_rsv
# 				# stats.conFacHistory = Array(Float64,0)
# 				# stats.setupT = 0.0
# 				# stats.solveT = 0.0
# 				# stats.converged = false
# 				# b1 = copy(b)
# 				# x1 = copy(x)
# 				# rnk = ones(n)
# 				# time_0 = time()
# 				# lv = hierarchy_test(A,A',vec(rnk),x1,b1,opts,aggFunc)
# 				# stats.setupT = (time() - time_0)

# 				# time_0 = time();
# 				# lv = find_rsv_RS_unweighted!(lv, lv[1].rnk, opts,stats)
# 				# stats.rsvT = (time() - time_0);
# 				# time_0 = time()
# 				# solve!(lv,stats,opts)
# 				# stats.solveT = (time() - time_0)

# 				# println(string("results_RS_vs_LAMG_mix/",ln[1:search(ln,".mtx")[1]-1],".txt"))
# 				# PrintStatistics(string(dir_results,string(aggFunc),"/",ln[1:search(ln,".mtx")[1]-1],".txt"),lv,stats,opts)
# 				# lv = 0; 

# 				# aggFunc = agg_CR
# 				# stats.conFacHistory = Array(Float64,0)
# 				# stats.setupT = 0.0
# 				# stats.solveT = 0.0
# 				# stats.converged = false
# 				# b1 = copy(b)
# 				# x1 = copy(x)
# 				# rnk = ones(n)

# 				# time_0 = time();
# 				# rnk = find_rsv_RS!( A,rnk,opts ,stats);
# 				# stats.rsvT = (time() - time_0);

# 				# time_0 = time()
# 				# lv = hierarchy_test(A,A',vec(rnk),x1,b1,opts,aggFunc)
# 				# stats.setupT = (time() - time_0)
# 				# time_0 = time()
# 				# solve!(lv,stats,opts)
# 				# stats.solveT = (time() - time_0)

# 				# println(string("results_RS_vs_LAMG_mix/",ln[1:search(ln,".mtx")[1]-1],".txt"))
# 				# PrintStatistics(string(dir_results,string(aggFunc),"/",ln[1:search(ln,".mtx")[1]-1],".txt"),lv,stats,opts)
# 				# lv = 0; 
# 				stats = setSolveStats();

# 				aggFunc = agg_rsv
# 				stats.conFacHistory = Array(Float64,0)
# 				stats.setupT = 0.0
# 				stats.solveT = 0.0
# 				stats.converged = false
# 				b1 = copy(b)
# 				x1 = copy(x)
# 				srand(1234)
# 				rnk = rand(n)

# 				time_0 = time();
# 				rnk = find_rsv_RS!( A,rnk,opts ,stats);
# 				stats.rsvT = (time() - time_0);

# 				time_0 = time()
# 				lv = hierarchy_noElim(A,A',vec(rnk),x1,b1,opts,aggFunc)
# 				stats.setupT = (time() - time_0)
# 				time_0 = time()
# 				solveV!(lv,stats,opts)
# 				stats.solveT = (time() - time_0)

# 				println(string(dir_results,string(aggFunc),"/","uniformchain1.txt"))
# 				PrintStatistics(string(dir_results,string(aggFunc),"/",ln[1:search(ln,".mtx")[1]-1],".txt"),lv,stats,opts)
# 									# PrintStatistics(string(dir_results,string(aggFunc),"/","uniformchain1.txt"),lv,stats,opts)

# 				lv = 0; 

				# ######## GMRES #######
				# b1 = copy(b)
				# x1 = copy(x)
				# 	GMRES(A,b1,x1,dir_results,opts,1000)

	end
# end



