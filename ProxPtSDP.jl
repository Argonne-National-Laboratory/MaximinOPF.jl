#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

MAX_TIME=24*3600

α_val = zeros(nbuses)
β_val = zeros(nbuses)
γ_val = zeros(nbuses)
δ_val = zeros(nbuses)
ζpUB_val = zeros(nbuses)
ζpLB_val = zeros(nbuses)
ζqUB_val = zeros(nbuses)
ζqLB_val = zeros(nbuses)
x_val = zeros(nlines)
W_val = zeros(nbuses)
Wr_val = zeros(nlines)
Wi_val = zeros(nlines)
e_val = zeros(nbuses)
f_val = zeros(nbuses)
objval = 0
bestBd = 0
η0Val = 0
nCuts = 0
new_cut = false
MAX_N_CUTS = 10000
maxNSG = 1
sg_α = zeros(nbuses,MAX_N_CUTS)
sg_β = zeros(nbuses,MAX_N_CUTS)
sg_γ = zeros(nbuses,MAX_N_CUTS) 
sg_δ = zeros(nbuses,MAX_N_CUTS) 

function solveMP(firstIt=false)
  global α_val, β_val, γ_val, δ_val
  global x_val
  global nCuts
  global sg_α, sg_β, sg_γ, sg_δ
  global objval

# The master problem MP
  mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=0,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_THREADS=1))
  @variable(mMP, -1 <= α[i=N] <= 1)
  @variable(mMP, -1 <= β[i=N] <= 1)
  @variable(mMP, γp[i=N] >= 0)
  @variable(mMP, γm[i=N] >= 0)
  @constraint(mMP, [i=N], γp[i]+γm[i] <= 1)
  @variable(mMP, ζpUB[g=G] >= 0)
  @variable(mMP, ζpLB[g=G] >= 0)
  @variable(mMP, ζqUB[g=G] >= 0)
  @variable(mMP, ζqLB[g=G] >= 0)
  @variable(mMP, 0 <= x[l=L] <= 1)
  @constraint(mMP, sum(x[l] for l in L) <= K)
  fixX=zeros(Int,nlines)
  #fixX[208]=1
  fixX[183]=1
  #fixX[41]=1
  #fixX[80]=1
  for l in L
    setlowerbound(x[l],fixX[l])
    setupperbound(x[l],fixX[l])
  end
  
  for i in N
   for g in BusGeners[i]
    @constraint(mMP, 
	-α[i] + ζpUB[g] - ζpLB[g] == 0 )
    @constraint(mMP, 
	-β[i] + ζqUB[g] - ζqLB[g] == 0 ) 
   end
  end

  if firstIt
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γm[i]*Wmin[i]-γp[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N)
    )
  else
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γm[i]*Wmin[i]-γp[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N)
	- 5*(sum( (α_val[i] -α[i])^2 + (β_val[i]- β[i])^2 + (γ_val[i] - γm[i])^2 + (δ_val[i] - γp[i])^2 for i in N) 
	+ sum( (ζpUB_val[g] - ζpUB[g])^2 + (ζpLB_val[g] - ζpLB[g])^2 for g in G) 
	+ sum( (ζqUB_val[g] - ζqUB[g])^2 + (ζqLB_val[g] - ζqLB[g])^2 for g in G) ) 
    )
  end

  if nCuts > 0
	@constraint(mMP, cutConstr[s=1:nCuts],  sum( sg_α[i,s] * α[i] for i in N) + sum(sg_β[i,s] * β[i] for i in N )  
		+ sum(sg_γ[i,s] * γm[i] for i in N)  + sum(sg_δ[i,s] * γp[i] for i in N) >= 0)
  end

  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.


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
	for l in L
    		x_val[l] = round(getvalue(x[l]))
	end
        objval=getobjectivevalue(mMP)
        getsolvetime(mMP)
        #println(" with optimal value ",solninfo[nlines+1])
        if status == :Stall
	  println("solveNodeAC: Return status $status")
        end	
  else
	println("solveNodeAC: Return status $status")
  end
end #end of function

function testECP()
   global objval, nCuts, new_cut
   start_time = time_ns()
   solveMP(true)
@show objval
   for kk=1:MAX_N_CUTS
     solveEta0Eigs()
     if new_cut
	solveMP()
@show nCuts
@show objval,η0Val
     else
	break
     end

   end
   end_time = time_ns()
   runtime = (end_time-start_time)/1e9
@show runtime
end

# Define callback function for generating and adding cuts
function solveEta0Eigs()

	global nSG
	global η0Val
	global α_val, β_val, γ_val, δ_val
	global W_val, Wr_val, Wi_val
	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
        global nCuts, new_cut

	new_cut = false
	H=spzeros(2*nbuses,2*nbuses)
        for i in N
          H[i,i] +=  α_val[i] * acYshR[i] - β_val[i] * acYshI[i]  + δ_val[i] - γ_val[i]
          H[nbuses+i,nbuses+i] += α_val[i] * acYshR[i] - β_val[i] * acYshI[i] + δ_val[i] - γ_val[i]
        end
        for l in L
          from = busIdx[lines[l].from];to = busIdx[lines[l].to] 
	  if x_val[l] == 0
            H[from,from] += α_val[from] * acYffR[l] - β_val[from] * acYffI[l]
            H[nbuses+from,nbuses+from] += α_val[from] * acYffR[l] - β_val[from] * acYffI[l]
            H[to,to] += α_val[to] * acYttR[l] - β_val[to] * acYttI[l]
            H[nbuses+to,nbuses+to] += α_val[to] * acYttR[l] - β_val[to] * acYttI[l]
            H[from,to] += 0.5*( α_val[from] * acYftR[l] - β_val[from] * acYftI[l] + α_val[to] * acYtfR[l] - β_val[to] * acYtfI[l] )
            H[to,from] += 0.5*( α_val[from] * acYftR[l] - β_val[from] * acYftI[l] + α_val[to] * acYtfR[l] - β_val[to] * acYtfI[l] )
            H[nbuses+from, nbuses+to] += 0.5*( α_val[from] * acYftR[l] - β_val[from] * acYftI[l] + α_val[to] * acYtfR[l] - β_val[to] * acYtfI[l] )
            H[nbuses+to, nbuses+from] += 0.5*( α_val[from] * acYftR[l] - β_val[from] * acYftI[l] + α_val[to] * acYtfR[l] - β_val[to] * acYtfI[l] )
            H[to, nbuses+from] += 0.5*( α_val[from] * acYftI[l] - α_val[to] * acYtfI[l] + β_val[from] * acYftR[l] - β_val[to] * acYtfR[l] )
            H[nbuses+from, to] += 0.5*( α_val[from] * acYftI[l] - α_val[to] * acYtfI[l] + β_val[from] * acYftR[l] - β_val[to] * acYtfR[l] )
            H[from,nbuses+to] -= 0.5*( α_val[from] * acYftI[l] - α_val[to] * acYtfI[l] + β_val[from] * acYftR[l] - β_val[to] * acYtfR[l] )
            H[nbuses+to,from] -= 0.5*( α_val[from] * acYftI[l] - α_val[to] * acYtfI[l] + β_val[from] * acYftR[l] - β_val[to] * acYtfR[l] )
	  end
        end
	# Reference bus angle is zero
	e_val = zeros(nbuses)
	f_val = zeros(nbuses)
	E=eigs(H,nev=maxNSG,which=:SR, maxiter=100000, tol=1e-8)

	η0Val = E[1][1]
	nSG = 0
	for cc=1:maxNSG
         
	 if E[1][cc] <= -TOL
	  nCuts += 1
	  new_cut = true
	  for i in N
 	    e_val[i] = E[2][i,cc]
	    f_val[i] = E[2][nbuses+i,cc]
	  end
	  
	  for i in N
	    W_val[i] = e_val[i]^2 + f_val[i]^2
            sg_α[i,nCuts] = acYshR[i] * W_val[i] 
            sg_β[i,nCuts] = -acYshI[i] * W_val[i] 
	    sg_δ[i,nCuts] = W_val[i] 
	    sg_γ[i,nCuts] = -W_val[i]
	  end
	  for l in L
	    from = busIdx[lines[l].from]; to = busIdx[lines[l].to]
	    if x_val[l] == 0
	      Wr_val[l] = e_val[from]*e_val[to] + f_val[from]*f_val[to]
	      Wi_val[l] = e_val[to]*f_val[from] - e_val[from]*f_val[to]
	      sg_α[from,nCuts] += (acYffR[l] * W_val[from] + acYftR[l] * Wr_val[l] + acYftI[l] * Wi_val[l])
	      sg_α[to,nCuts] += (acYttR[l] * W_val[to] + acYtfR[l] * Wr_val[l] - acYtfI[l] * Wi_val[l])
	      sg_β[from,nCuts] += (-acYffI[l] * W_val[from] - acYftI[l] * Wr_val[l] + acYftR[l]* Wi_val[l])
	      sg_β[to,nCuts] += (-acYttI[l] * W_val[to] - acYtfI[l] * Wr_val[l] - acYtfR[l] * Wi_val[l])
	    end
	  end
	else
	  break
	end
       end # s=1:maxNSG
end

