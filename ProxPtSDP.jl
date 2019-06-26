#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

  # define data
  lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
  nbuses, nlines, ngens = length(buses), length(lines), length(generators)
  N = 1:nbuses; L = 1:nlines; G = 1:ngens
  # build a dictionary between buses ids and their indexes
  busIdx = mapBusIdToIdx(buses)
  # set up the fromLines and toLines for each bus
  fromLines, toLines = mapLinesToBuses(buses, lines, busIdx)
  fromBus=zeros(Int,nlines); toBus=zeros(Int,nlines)
  for l in L
    fromBus[l] = busIdx[lines[l].from]; toBus[l] = busIdx[lines[l].to] 
  end
  # generators at each bus
  BusGeners = mapGenersToBuses(buses, generators, busIdx)

  Y = Dict()  # Admittances
  Y["ffR"] = opfdata.admittancesAC.YffR; Y["ffI"] = opfdata.admittancesAC.YffI;
  Y["ttR"] = opfdata.admittancesAC.YttR; Y["ttI"] = opfdata.admittancesAC.YttI;
  Y["ftR"] = opfdata.admittancesAC.YftR; Y["ftI"] = opfdata.admittancesAC.YftI;
  Y["tfR"] = opfdata.admittancesAC.YtfR; Y["tfI"] = opfdata.admittancesAC.YtfI;
  Y["shR"] = opfdata.admittancesAC.YshR; Y["shI"] = opfdata.admittancesAC.YshI;

  PD = zeros(nbuses); QD = zeros(nbuses)
  Wmin = zeros(nbuses); Wmax = zeros(nbuses)
  for i in N
	PD[i] = buses[i].Pd / baseMVA; QD[i] = buses[i].Qd / baseMVA
	Wmin[i] = (buses[i].Vmin)^2; Wmax[i] = (buses[i].Vmax)^2
  end
  Pmin = zeros(ngens); Pmax = zeros(ngens)
  Qmin = zeros(ngens); Qmax = zeros(ngens)
  for g in G
	Pmin[g] = generators[g].Pmin; Pmax[g] = generators[g].Pmax
	Qmin[g] = generators[g].Qmin; Qmax[g] = generators[g].Qmax
  end

  println("Done with initial setup.")

MAX_TIME=24*3600
MAX_N_CUTS = 2
MAX_N_ITERS = 30000
M=100

α_val = zeros(nbuses)
β_val = zeros(nbuses)
γ_val = zeros(nbuses)
δ_val = zeros(nbuses)
ζpUB_val = zeros(nbuses)
ζpLB_val = zeros(nbuses)
ζqUB_val = zeros(nbuses)
ζqLB_val = zeros(nbuses)
α_ctr = zeros(nbuses)
β_ctr = zeros(nbuses)
γ_ctr = zeros(nbuses)
δ_ctr = zeros(nbuses)
ζpUB_ctr = zeros(nbuses)
ζpLB_ctr = zeros(nbuses)
ζqUB_ctr = zeros(nbuses)
ζqLB_ctr = zeros(nbuses)
α_trl = zeros(nbuses,MAX_N_CUTS)
β_trl = zeros(nbuses,MAX_N_CUTS)
γ_trl = zeros(nbuses,MAX_N_CUTS)
δ_trl = zeros(nbuses,MAX_N_CUTS)
ζpUB_trl = zeros(nbuses,MAX_N_CUTS)
ζpLB_trl = zeros(nbuses,MAX_N_CUTS)
ζqUB_trl = zeros(nbuses,MAX_N_CUTS)
ζqLB_trl = zeros(nbuses,MAX_N_CUTS)
zK = zeros(MAX_N_ITERS)
sigK = zeros(MAX_N_ITERS)
x_val = zeros(Int,nlines)
W_val = zeros(nbuses)
Wr_val = zeros(nlines)
Wi_val = zeros(nlines)
e_val = zeros(nbuses)
f_val = zeros(nbuses)
objval = 0
linobjval = 0
bestBd = 0
η0Val = 0
etaCtr = 0
etaTrl = zeros(MAX_N_CUTS)
nCuts = 0
new_cut = false
maxNSG = 1
SSC = 0.05
m2 = 0.9
TOL = 1e-8
sg_α = zeros(nbuses,MAX_N_CUTS)
sg_β = zeros(nbuses,MAX_N_CUTS)
sg_γ = zeros(nbuses,MAX_N_CUTS) 
sg_δ = zeros(nbuses,MAX_N_CUTS) 
sg_α_ctr = zeros(nbuses,MAX_N_CUTS)
sg_β_ctr = zeros(nbuses,MAX_N_CUTS)
sg_γ_ctr = zeros(nbuses,MAX_N_CUTS) 
sg_δ_ctr = zeros(nbuses,MAX_N_CUTS) 
cutDuals = zeros(MAX_N_CUTS)
linErrors = zeros(MAX_N_CUTS)
tMAX = 2.5
tMIN = 0.5
tVal = 1
tCntr = 0

function solveMP(firstIt=false)
  global α_val, β_val, γ_val, δ_val
  global α_trl, β_trl, γ_trl, δ_trl
  global x_val
  global nCuts, linErrors
  global η0Val, etaCtr, etaTrl
  global sg_α, sg_β, sg_γ, sg_δ
  global objval, linobjval, iter
  global sigK, zK

# The master problem MP
  #mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_THREADS=1,))
  mMP = Model(solver=IpoptSolver())
  #mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=4))
  @variable(mMP, -1 <= α[i=N] <= 1, start=0)
  @variable(mMP, -1 <= β[i=N] <= 1, start=0)
  @variable(mMP, γp[i=N] >= 0, start=0)
  @variable(mMP, γm[i=N] >= 0, start=0)
  @constraint(mMP, [i=N], γp[i]+γm[i] <= 1)
  @variable(mMP, 0 <= ζpUB[g=G] <= 1, start=0)
  @variable(mMP, 0 <= ζpLB[g=G] <= 1, start=0)
  @variable(mMP, 0 <= ζqUB[g=G] <= 1, start=0)
  @variable(mMP, 0 <= ζqLB[g=G] <= 1, start=0)
  #@variable(mMP, 0 <= x[l=L] <= 1)
  #@constraint(mMP, sum(x[l] for l in L) <= K)
  #x_val[208]=1
  #x_val[183]=1
  x_val[41]=1
  x_val[80]=1
#=
  x_val[60]=1
  x_val[65]=1
  x_val[66]=1
  x_val[72]=1
=#
#=
  for l in L
    setlowerbound(x[l],fixX[l])
    setupperbound(x[l],fixX[l])
  end
=#
  
  for i in N
   for g in BusGeners[i]
    @constraint(mMP, -α[i] + ζpUB[g] - ζpLB[g] == 0 )
    @constraint(mMP, -β[i] + ζqUB[g] - ζqLB[g] == 0 ) 
   end
  end

#=
  if firstIt
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γm[i]*Wmin[i]-γp[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N)
    )
  else
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γm[i]*Wmin[i]-γp[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N)
	- (0.5/tVal)*(sum( (α_ctr[i] -α[i])^2 + (β_ctr[i]- β[i])^2 + (γ_ctr[i] - γm[i])^2 + (δ_ctr[i] - γp[i])^2 for i in N) 
	#+ sum( (ζpUB_ctr[g] - ζpUB[g])^2 + (ζpLB_ctr[g] - ζpLB[g])^2 for g in G) 
	#+ sum( (ζqUB_ctr[g] - ζqUB[g])^2 + (ζqLB_ctr[g] - ζqLB[g])^2 for g in G) 
			) 
    )
  end
=#


#=
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γm[i]*Wmin[i]-γp[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N)
    )
=#
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γm[i]*Wmin[i]-γp[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N)
	- (0.5/tVal)*(sum( (α_ctr[i] -α[i])^2 + (β_ctr[i]- β[i])^2 + (γ_ctr[i] - γm[i])^2 + (δ_ctr[i] - γp[i])^2 for i in N) 
	+ sum( (ζpUB_ctr[g] - ζpUB[g])^2 + (ζpLB_ctr[g] - ζpLB[g])^2 for g in G) 
	+ sum( (ζqUB_ctr[g] - ζqUB[g])^2 + (ζqLB_ctr[g] - ζqLB[g])^2 for g in G) 
			) 
    )
    @expression(mMP, linObjExpr, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γm[i]*Wmin[i]-γp[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N)
    )

#=
  @constraint(mMP, apriori[i=N], (YshR[i] + sum( YffR[l] for l in fromLines[i]) + sum( YttR[l] for l in toLines[i])  )*α[i] 
		- (YshI[i] + sum(YffI[l] for l in fromLines[i]) + sum(YttI[l] for l in toLines[i]  ))*β[i] + γp[i] - γm[i] >= 0  
	     )
=#

  if nCuts > 0
#=
@show sg_α[1:nbuses,nCuts]
@show sg_β[1:nbuses,nCuts]
@show sg_γ[1:nbuses,nCuts]
@show sg_δ[1:nbuses,nCuts]
=#
	@constraint(mMP, cutConstr[s=1:nCuts],  sum( sg_α[i,s] * α[i] for i in N) + sum(sg_β[i,s] * β[i] for i in N )  
		+ sum(sg_γ[i,s] * γm[i] for i in N)  + sum(sg_δ[i,s] * γp[i] for i in N) >= 0)
  end


  status=solve(mMP)
  if status == :Optimal || status == :CPX_STAT_NUM_BEST
	# get the variable values
	for i in N
		α_val[i] = getvalue(α[i])
		β_val[i] = getvalue(β[i])
		γ_val[i] = getvalue(γm[i])
		δ_val[i] = getvalue(γp[i]) 
	end
	for i in G
		ζpUB_val[i] = getvalue(ζpUB[i])
		ζpLB_val[i] = getvalue(ζpLB[i])
		ζqUB_val[i] = getvalue(ζqUB[i])
		ζqLB_val[i] = getvalue(ζqLB[i])
	end
#=
	for l in L
    		x_val[l] = round(getvalue(x[l]))
	end
=#
        objval=getobjectivevalue(mMP)
        linobjval=getvalue(linObjExpr)
        #getsolvetime(mMP)
	if status == :Optimal
if nCuts > 0
	  nMax=nCuts
	  #for n=nMax:-1:1
	  for n=1:nCuts
	    cutDuals[n] = -getdual(cutConstr[n])
	    linErrors[n] =  etaCtr - (etaTrl[n] + sum( sg_α[i,n]*(α_ctr[i]-α_trl[i,n])  for i in N) + sum(sg_β[i,n] * (β_ctr[i] - β_trl[i,n]) for i in N )  
			+ sum(sg_γ[i,n] * (γ_ctr[i]-γ_trl[i,n]) for i in N)  + sum(sg_δ[i,n] * (δ_ctr[i]-δ_trl[i,n]) for i in N)   ) 
#=
	    if abs(cutDuals[n]) < 1e-6
		etaTrl[n] = etaTrl[nCuts]
		for i in N
	    	    α_trl[i,n] = α_trl[i,nCuts]
	    	    β_trl[i,n] = β_trl[i,nCuts]
	            γ_trl[i,n] = γ_trl[i,nCuts]
	            δ_trl[i,n] = δ_trl[i,nCuts]
		    sg_α[i,n]=sg_α[i,nCuts]
		    sg_β[i,n]=sg_β[i,nCuts]
		    sg_γ[i,n]=sg_γ[i,nCuts]
		    sg_δ[i,n]=sg_δ[i,nCuts]
	    	    cutDuals[n] = cutDuals[nCuts]
		    linErrors[n] = linErrors[nCuts]
		end
		nCuts -= 1
	    end
=#
  	  end
#=
	    sigK[iter] = sum( cutDuals[n]*linErrors[n]  for n in 1:nCuts)
	    zK[iter] = 0.0
	    zK[iter] += sum( sum( (cutDuals[n]*sg_α[i,n])  for n in 1:nCuts)^1   for i in N)
	    zK[iter] += sum( sum( (cutDuals[n]*sg_β[i,n])  for n in 1:nCuts)^1   for i in N)
	    zK[iter] += sum( sum( (cutDuals[n]*sg_γ[i,n])  for n in 1:nCuts)^1   for i in N)
	    zK[iter] += sum( sum( (cutDuals[n]*sg_δ[i,n])  for n in 1:nCuts)^1   for i in N)
	    zK[iter] = sqrt(zK[iter])
=#
	  else
	    sigK[iter] = 0.0
	    zK[iter] = 0.0
	  end
#@show cutDuals[1:nCuts]
#@show linErrors[1:nCuts]
	end
        #println(" with optimal value ",solninfo[nlines+1])
        if status == :Stall
	  println("solveNodeAC: Return status $status")
        end	
  else
	println("solveNodeAC: Return status $status")
  end
  return status
end #end of function

function testECP()
   global objval, linobjval, nCuts, new_cut, SSC, m2, tVal, tCntr, etaCtr, iter, tMAX, tUB, tLB
   global α_ctr, β_ctr, γ_ctr, δ_ctr
   global ζpUB_ctr, ζpLB_ctr, ζqUB_ctr, ζqLB_ctr
   global α_val, β_val, γ_val, δ_val
   global ζpUB_val, ζpLB_val, ζqUB_val, ζqLB_val
   start_time = time_ns()
   for kk=1:MAX_N_ITERS
        iter=kk 
	status = solveMP()
        solveEta0Eigs(kk==1)
	dNorm = sqrt(norm(α_val-α_ctr)^2 + norm(β_val-β_ctr)^2 + norm(γ_val-γ_ctr)^2 + norm(δ_val-δ_ctr)^2 
		+ norm(ζpUB_val-ζpUB_ctr)^2 + norm(ζpLB_val-ζpLB_ctr)^2 + norm(ζqUB_val-ζqUB_ctr)^2 + norm(ζqLB_val-ζqLB_ctr)^2
		    )
        #if ( (1.0/tVal)*dNorm <= TOL && -(1.0/tVal)*dNorm^2 - η0Val <= TOL ) 
        if ( (1.0/tVal)*dNorm <= TOL && -η0Val <= TOL ) 
	   break
        end

   end
@show objval,linobjval,η0Val
   end_time = time_ns()
   runtime = (end_time-start_time)/1e9
@show runtime
end

# Define callback function for generating and adding cuts
function solveEta0Eigs(firstIt=false)

	global nSG, iter
	global η0Val, etaCtr, tVal, tCntr, tMAX, tMIN, objval, linobjval
	global α_val, β_val, γ_val, δ_val
	global ζpUB_val, ζpLB_val, ζqUB_val, ζqLB_val
	global α_ctr, β_ctr, γ_ctr, δ_ctr
	global ζpUB_ctr, ζpLB_ctr, ζqUB_ctr, ζqLB_ctr
	global W_val, Wr_val, Wi_val
	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
        global nCuts, new_cut
        global sigK, zK

	new_cut = false
	H=spzeros(2*nbuses,2*nbuses)
        for i in N
          H[i,i] +=  α_val[i] * Y["shR"][i] - β_val[i] * Y["shI"][i]  + δ_val[i] - γ_val[i]
          H[nbuses+i,nbuses+i] += α_val[i] * Y["shR"][i] - β_val[i] * Y["shI"][i] + δ_val[i] - γ_val[i]
        end
        for l in L
          from = busIdx[lines[l].from];to = busIdx[lines[l].to] 
	  if x_val[l] == 0
            H[from,from] += α_val[from] * Y["ffR"][l] - β_val[from] * Y["ffI"][l]
            H[nbuses+from,nbuses+from] += α_val[from] * Y["ffR"][l] - β_val[from] * Y["ffI"][l]
            H[to,to] += α_val[to] * Y["ttR"][l] - β_val[to] * Y["ttI"][l]
            H[nbuses+to,nbuses+to] += α_val[to] * Y["ttR"][l] - β_val[to] * Y["ttI"][l]
            H[from,to] += 0.5*( α_val[from] * Y["ftR"][l] - β_val[from] * Y["ftI"][l] + α_val[to] * Y["tfR"][l] - β_val[to] * Y["tfI"][l] )
            H[to,from] += 0.5*( α_val[from] * Y["ftR"][l] - β_val[from] * Y["ftI"][l] + α_val[to] * Y["tfR"][l] - β_val[to] * Y["tfI"][l] )
            H[nbuses+from, nbuses+to] += 0.5*( α_val[from] * Y["ftR"][l] - β_val[from] * Y["ftI"][l] + α_val[to] * Y["tfR"][l] - β_val[to] * Y["tfI"][l] )
            H[nbuses+to, nbuses+from] += 0.5*( α_val[from] * Y["ftR"][l] - β_val[from] * Y["ftI"][l] + α_val[to] * Y["tfR"][l] - β_val[to] * Y["tfI"][l] )
            H[to, nbuses+from] += 0.5*( α_val[from] * Y["ftI"][l] - α_val[to] * Y["tfI"][l] + β_val[from] * Y["ftR"][l] - β_val[to] * Y["tfR"][l] )
            H[nbuses+from, to] += 0.5*( α_val[from] * Y["ftI"][l] - α_val[to] * Y["tfI"][l] + β_val[from] * Y["ftR"][l] - β_val[to] * Y["tfR"][l] )
            H[from,nbuses+to] -= 0.5*( α_val[from] * Y["ftI"][l] - α_val[to] * Y["tfI"][l] + β_val[from] * Y["ftR"][l] - β_val[to] * Y["tfR"][l] )
            H[nbuses+to,from] -= 0.5*( α_val[from] * Y["ftI"][l] - α_val[to] * Y["tfI"][l] + β_val[from] * Y["ftR"][l] - β_val[to] * Y["tfR"][l] )
	  end
        end
	# Reference bus angle is zero
	e_val = zeros(nbuses)
	f_val = zeros(nbuses)
	E=eigs(H,nev=maxNSG,which=:SR, maxiter=100000, tol=1e-8)

	η0Val = E[1][1]
	nSG = 0
	for cc=1:maxNSG
         
	 #if E[1][cc] <= -TOL
	  if nCuts == MAX_N_CUTS
	    dualSum = sum( cutDuals[s] for s in 1:nCuts)
	    if dualSum > 1e-6
	      for i in N
	       sg_α[i,1] = sum( cutDuals[s]*sg_α[i,s] for s in 1:nCuts )  
	       sg_β[i,1] = sum( cutDuals[s]*sg_β[i,s] for s in 1:nCuts )  
	       sg_δ[i,1] = sum( cutDuals[s]*sg_δ[i,s] for s in 1:nCuts )  
	       sg_γ[i,1] = sum( cutDuals[s]*sg_γ[i,s] for s in 1:nCuts )  
	       sg_α[i,1] /= dualSum
	       sg_β[i,1] /= dualSum
	       sg_δ[i,1] /= dualSum
	       sg_γ[i,1] /= dualSum
	       α_trl[i,1] = sum( cutDuals[s]*α_trl[i,s] for s in 1:nCuts )  
	       β_trl[i,1] = sum( cutDuals[s]*β_trl[i,s] for s in 1:nCuts )  
	       δ_trl[i,1] = sum( cutDuals[s]*δ_trl[i,s] for s in 1:nCuts )  
	       γ_trl[i,1] = sum( cutDuals[s]*γ_trl[i,s] for s in 1:nCuts )  
	       α_trl[i,1] /= dualSum
	       β_trl[i,1] /= dualSum
	       δ_trl[i,1] /= dualSum
	       γ_trl[i,1] /= dualSum
	       etaTrl[1] = sum( cutDuals[s]*etaTrl[s] for s in 1:nCuts)
	       etaTrl[1] /= dualSum
	      end
	      nCuts = 1
	    else
	      nCuts = 0
	    end
	  end
	  nCuts += 1
	  new_cut = true
	  maxMod = 0
	  for i in N
 	    e_val[i] = E[2][i,cc]
	    f_val[i] = E[2][nbuses+i,cc]
	    if sqrt(e_val[i]^2 + f_val[i]^2) > maxMod
	        maxMod = e_val[i]^2 + f_val[i]^2
	    end
	  end
#@show e_val
#@show f_val
	  etaTrl[nCuts]=E[1][cc]
	  for i in N
	    W_val[i] = e_val[i]^2 + f_val[i]^2
            sg_α[i,nCuts] = Y["shR"][i] * W_val[i] 
            sg_β[i,nCuts] = -Y["shI"][i] * W_val[i] 
	    sg_δ[i,nCuts] = W_val[i] 
	    sg_γ[i,nCuts] = -W_val[i]
	    α_trl[i,nCuts] = α_val[i]
	    β_trl[i,nCuts] = β_val[i]
	    γ_trl[i,nCuts] = γ_val[i]
	    δ_trl[i,nCuts] = δ_val[i]
	  end
	  for l in L
	    from = busIdx[lines[l].from]; to = busIdx[lines[l].to]
	    if x_val[l] == 0
	      Wr_val[l] = e_val[from]*e_val[to] + f_val[from]*f_val[to]
	      Wi_val[l] = e_val[to]*f_val[from] - e_val[from]*f_val[to]
	      sg_α[from,nCuts] += (Y["ffR"][l] * W_val[from] + Y["ftR"][l] * Wr_val[l] + Y["ftI"][l] * Wi_val[l])
	      sg_α[to,nCuts] += (Y["ttR"][l] * W_val[to] + Y["tfR"][l] * Wr_val[l] - Y["tfI"][l] * Wi_val[l])
	      sg_β[from,nCuts] += (-Y["ffI"][l] * W_val[from] - Y["ftI"][l] * Wr_val[l] + Y["ftR"][l]* Wi_val[l])
	      sg_β[to,nCuts] += (-Y["ttI"][l] * W_val[to] - Y["tfI"][l] * Wr_val[l] - Y["tfR"][l] * Wi_val[l])
	    end
	  end
	sgAngle = sum( sg_α_ctr[i]*sg_α[i,nCuts] + sg_β_ctr[i]*sg_β[i,nCuts] + sg_γ_ctr[i]*sg_γ[i,nCuts] + sg_δ_ctr[i]*sg_δ[i,nCuts] for i in N)
	sgAngleNorm = sqrt(norm(sg_α_ctr)^2+norm(sg_β_ctr)^2+norm(sg_γ_ctr)^2+norm(sg_δ_ctr)^2  )*sqrt( norm(sg_α[N,nCuts])^2 + norm(sg_β[N,nCuts])^2 + norm(sg_γ[N,nCuts])^2 + norm(sg_δ[N,nCuts])^2 )
	if sgAngleNorm > 1e-4
	  sgAngle /= sgAngleNorm
	else
	  sgAngle = 0
	end
	sgNorm = sqrt(norm(sg_α[N,nCuts])^2 + norm(sg_β[N,nCuts])^2 + norm(sg_γ[N,nCuts])^2 + norm(sg_δ[N,nCuts])^2 )

#=
        for i in N
            #sg_α[i,nCuts] /= sgNorm
	    if abs(sg_α[i,nCuts]) < 1e-10
              sg_α[i,nCuts] = 0
	    end
            #sg_β[i,nCuts] /= sgNorm 
	    if abs(sg_β[i,nCuts]) < 1e-10
              sg_β[i,nCuts] = 0
	    end
	    #sg_δ[i,nCuts] /= sgNorm 
	    if abs(sg_δ[i,nCuts]) < 1e-10
              sg_δ[i,nCuts] = 0
	    end
	    #sg_γ[i,nCuts] /= sgNorm
	    if abs(sg_γ[i,nCuts]) < 1e-10
              sg_γ[i,nCuts] = 0
	    end
	end
=#
#@show sgAngle
	  SSCVal= -(η0Val-etaCtr)/etaCtr 
	  tVal0 = tVal
	  if SSCVal > SSC || iter <= 1
	      gTd = sum( sg_α[i,nCuts]*α_val[i]  for i in N) + sum(sg_β[i,nCuts] * β_val[i] for i in N )  
		+ sum(sg_γ[i,nCuts] * γ_val[i] for i in N)  + sum(sg_δ[i,nCuts] * δ_val[i] for i in N) 
	        #println("Serious step")
   	        etaCtr = η0Val
	        for i in N
	         α_ctr[i] = α_val[i]
	         β_ctr[i] = β_val[i]
	         γ_ctr[i] = γ_val[i]
	         δ_ctr[i] = δ_val[i]
		 sg_α_ctr[i] = sg_α[i,nCuts] 
	 	 sg_β_ctr[i] = sg_β[i,nCuts] 
		 sg_γ_ctr[i] = sg_γ[i,nCuts] 
		 sg_δ_ctr[i] = sg_δ[i,nCuts] 
	        end
	        for g in G
	         ζpUB_ctr[g] = ζpUB_val[g]
	         ζpLB_ctr[g] = ζpLB_val[g]
	         ζqUB_ctr[g] = ζqUB_val[g]
	         ζqLB_ctr[g] = ζqLB_val[g]
	        end
########  Adjust penalty under serious step ###############
	        hk = 2*tVal*(1.0 - SSCVal)
		tProv = tVal0
		if SSCVal >= 0.5 && tCntr > 0
		  tProv = hk
		elseif tCntr > 3
		  tProv = 0.5*tVal0
		end
		tVal = max(max(tProv,0.1*tVal0),tMIN)
		if abs(tVal-tVal0)<1e-6
		  tCntr = max(tCntr+1,1)
		else
		  tCntr = 1
		end
######## End adjustment of penalty under serious step ################
@show objval,linobjval,etaCtr
	  else

########  Adjust penalty under null step ###############
	    hk = 2*tVal*(1.0 - SSCVal)
	    linErrorVal =  etaTrl[nCuts] - etaCtr + sum( sg_α[i,nCuts]*(α_trl[i,nCuts]-α_ctr[i])  for i in N) + sum(sg_β[i,nCuts] * (β_trl[i,nCuts] - β_ctr[i]) for i in N )  
			+ sum(sg_γ[i,nCuts] * (γ_trl[i,nCuts]-γ_ctr[i]) for i in N)  + sum(sg_δ[i,nCuts] * (δ_trl[i,nCuts]-δ_ctr[i]) for i in N)    
	    proxTermNorm = sqrt(norm(α_trl[N,nCuts]-α_ctr[N])^2  + norm(β_trl[N,nCuts] - β_ctr[N])^2 + norm(γ_trl[N,nCuts]-γ_ctr[N])^2 + norm(δ_trl[N,nCuts]-δ_ctr[N])^2  ) 
	    tProv = tVal0
	    if linErrorVal > max( -etaCtr + tVal0*proxTermNorm^2 + proxTermNorm, 10*etaCtr) && tCntr < -3
		tProv = hk
	    end
	    tVal = min(tProv,10*tVal0)
	    if abs(tVal-tVal0) < 1e-6
		tCntr = min(tCntr-1,-1)
	    else
		tCntr = -1
	    end
######## End adjustment of penalty under null step ################
	  end
       end # s=1:maxNSG
end

