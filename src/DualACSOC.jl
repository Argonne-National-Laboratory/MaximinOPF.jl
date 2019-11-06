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
λF_val = zeros(nlines)
μF_val = zeros(nlines)
λT_val = zeros(nlines)
μT_val = zeros(nlines)
solninfo=zeros(nlines+2)
TOL=1e-4

# The dual problem 
  global α_val, β_val, γ_val, δ_val
  global λF_val, μF_val, λT_val, μT_val

  #mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=5,MSK_IPAR_NUM_THREADS=4))
  #mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_MIQCPSTRAT=1))
  mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50))

  @variable(mMP, -1 <= α[i=N] <= 1)
  @variable(mMP, -1 <= β[i=N] <= 1)
  @variable(mMP, γp[i=N] >= 0)
  @variable(mMP, γm[i=N] >= 0)
  @constraint(mMP, [i=N], γp[i]+γm[i] <= 1)
  @variable(mMP, ζpUB[g=G] >=0)
  @variable(mMP, ζpLB[g=G] >=0)
  @variable(mMP, ζqUB[g=G] >=0)
  @variable(mMP, ζqLB[g=G] >=0)
  @variable(mMP, λF[l=L])
  @variable(mMP, λT[l=L])
  @variable(mMP, μF[l=L])
  @variable(mMP, μT[l=L])
  @variable(mMP, x[l=L], Bin)
  @constraint(mMP, sum(x[l] for l in L) <= K)
  
  for i in N
   for g in BusGeners[i]
    @constraint(mMP, 
	-α[i] + ζpUB[g] - ζpLB[g] == 0 )
    @constraint(mMP, 
	-β[i] + ζqUB[g] - ζqLB[g] == 0 ) 
   end
  end

  @expression(mMP, C[i=1:(2*nbuses),j=1:(2*nbuses)], 0)
  for i in N
	C[i,i] += γp[i] - γm[i] + α[i]*acYshR[i] - β[i]*acYshI[i] 
        C[nbuses+i,nbuses+i] += γp[i] - γm[i] + α[i]*acYshR[i] - β[i]*acYshI[i] 
    for l in fromLines[i]
      C[i,i] += λF[l]*acYffR[l] - μF[l]*acYffI[l]  
      C[nbuses+i,nbuses+i] += λF[l]*acYffR[l] - μF[l]*acYffI[l]  
    end
    for l in toLines[i]
      C[i,i] += λT[l]*acYttR[l] - μT[l]*acYttI[l]  
      C[nbuses+i,nbuses+i] += λT[l]*acYttR[l] - μT[l]*acYttI[l]  
    end
  end
  for l in L
      from=busIdx[lines[l].from]; to=busIdx[lines[l].to]
      C[from,nbuses+to] -= 0.5*(  λF[l]*acYftI[l] + μF[l]*acYftR[l] - λT[l]*acYtfI[l] - μT[l]*acYtfR[l]  )   
      C[to,nbuses+from] += 0.5*(  λF[l]*acYftI[l] + μF[l]*acYftR[l] - λT[l]*acYtfI[l] - μT[l]*acYtfR[l]  )   
      C[from,to] += 0.5*(λF[l]*acYftR[l]  - μF[l]*acYftI[l] + λT[l]*acYtfR[l] - μT[l]*acYtfI[l])   
      C[nbuses+from,nbuses+to] += 0.5*(λF[l]*acYftR[l]  - μF[l]*acYftI[l] + λT[l]*acYtfR[l] - μT[l]*acYtfI[l])   
  end

  @constraint(mMP, SOCApprox[l=L], 
	norm([ 0.5*(C[busIdx[lines[l].from],busIdx[lines[l].from]] - C[busIdx[lines[l].to],busIdx[lines[l].to]]);C[busIdx[lines[l].from],busIdx[lines[l].to]];C[busIdx[lines[l].from],nbuses+busIdx[lines[l].to]]    ]) - 0.5*(C[busIdx[lines[l].from],busIdx[lines[l].from]] + C[busIdx[lines[l].to],busIdx[lines[l].to]])  <= 0)


  @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γm[i]*Wmin[i]-γp[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N))

  rel=1
  @constraint(mMP, [l in L], α[busIdx[lines[l].from]] - rel*x[l] - λF[l] <= 0)
  @constraint(mMP, [l in L], α[busIdx[lines[l].from]] + rel*x[l] - λF[l] >= 0)
  @constraint(mMP, [l in L], -rel*(1 - x[l]) - λF[l] <= 0)
  @constraint(mMP, [l in L], rel*(1 - x[l]) - λF[l] >= 0)
  @constraint(mMP, [l in L], α[busIdx[lines[l].to]] - rel*x[l] - λT[l] <= 0)
  @constraint(mMP, [l in L], α[busIdx[lines[l].to]] + rel*x[l] - λT[l] >= 0)
  @constraint(mMP, [l in L], -rel*(1 - x[l]) - λT[l] <= 0)
  @constraint(mMP, [l in L], rel*(1 - x[l]) - λT[l] >= 0)

  @constraint(mMP, [l in L], β[busIdx[lines[l].from]] - rel*x[l] - μF[l] <= 0)
  @constraint(mMP, [l in L], β[busIdx[lines[l].from]] + rel*x[l] - μF[l] >= 0)
  @constraint(mMP, [l in L], -rel*(1 - x[l]) - μF[l] <= 0)
  @constraint(mMP, [l in L], rel*(1 - x[l]) - μF[l] >= 0)
  @constraint(mMP, [l in L], β[busIdx[lines[l].to]] - rel*x[l] - μT[l] <= 0)
  @constraint(mMP, [l in L], β[busIdx[lines[l].to]] + rel*x[l] - μT[l] >= 0)
  @constraint(mMP, [l in L], -rel*(1 - x[l]) - μT[l] <= 0)
  @constraint(mMP, [l in L], rel*(1 - x[l]) - μT[l] >= 0)

  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.

  if HEUR == 1 
  	@constraint(mMP, LambdaMuConstr1[l in L], λF[l]*YftI[l] - λT[l]*YtfI[l] + μF[l]*YftR[l] - μT[l]*YtfR[l] == 0.0)
  elseif HEUR == 2
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l] == 0) 
	@constraint(mMP, muFequalsT[l in L], μF[l] - μT[l] == 0)
  elseif HEUR == 3
	@constraint(mMP, LambdaMuConstr2[l in L], λF[l]*YtfR[l] - λT[l]*YftR[l] - μF[l]*YtfI[l] + μT[l]*YftI[l] == 0.0)
  end
  
function generateCuts(cb)
  global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
  global α_val, β_val, γ_val, δ_val
  global λF_val, μF_val, λT_val, μT_val
  global η0Val, ncuts
      for i in N
	α_val[i]=getvalue(α[i])
	β_val[i]=getvalue(β[i]) 
	γ_val[i]=getvalue(γm[i]) 
	δ_val[i]=getvalue(γp[i])
      end
      for l in L
	λF_val[l]=getvalue(λF[l]) 
	μF_val[l]=getvalue(μF[l])
	λT_val[l]=getvalue(λT[l])
	μT_val[l]=getvalue(μT[l])
      end
  try
    solveEta0Eigs()
    #solveEta0SDP()
  catch exc
    println("Exception caught with eigs(), solving η0Val subproblem with Ipopt as recourse.") 
    println(exc)
    solveEta0SDP()
  end
#=
println(α_val)
println(β_val) 
println(γ_val) 
println(δ_val)
println(λF_val) 
println(μF_val)
println(λT_val)
println(μT_val)
=#

  if η0Val <= -TOL 
	 for s=1:nSG
	    @lazyconstraint(cb, 0.0 <= sum( sg_α[i,s]*α[i] for i in N) + sum(sg_β[i,s]*β[i] for i in N )  
		+ sum(sg_γ[i,s]*γm[i] for i in N)  + sum(sg_δ[i,s]*γp[i] for i in N)
		+ sum(sg_λF[l,s]*λF[l] + sg_λT[l,s]*λT[l] for l in L) + sum(sg_μF[l,s]*μF[l] + sg_μT[l,s]*μT[l] for l in L),
		localcut=useLocalCuts)
	    ncuts += 1
	 end
  else
    println("Tolerance met for not generating a new lazy cut.")
  end
end
addlazycallback(mMP, generateCuts,fractional=false)


function testDualACSOC(x_soln)
  global solninfo
  bestUBVal=1e20
  nNodes=0
  incVal=0
  runtime=0

  status=solve(mMP)
  if status == :Optimal || status == :Stall
      for l in L
        x_soln[l] = round(getvalue(x[l]))
      end
      bestUBVal=getobjectivebound(mMP)
      incVal = getobjectivevalue(mMP)
      runtime = getsolvetime(mMP)
      nNodes=getnodecount(mMP)
printX(x_soln)
      #println(" with optimal value ",solninfo[nlines+1])
      if status == :Stall
	println("solveNodeACSOC: Return status $status")
      end	
  else
	println("solveNodeACSOC: Return status $status")
  end

  return bestUBVal,nNodes,incVal,runtime
  #primalObjval = solveFullModelSDP(x_soln)
  #println("Primal value: ",primalObjval," and dual value: ",dualObjval)
end



sg_α = zeros(nbuses,2*nbuses)
sg_β = zeros(nbuses,2*nbuses)
sg_γ = zeros(nbuses,2*nbuses) 
sg_δ = zeros(nbuses,2*nbuses) 
sg_λF = zeros(nlines,2*nbuses)
sg_μF = zeros(nlines,2*nbuses)
sg_λT = zeros(nlines,2*nbuses)
sg_μT = zeros(nlines,2*nbuses)
W_val = zeros(nbuses)
Wr_val = zeros(nlines)
Wi_val = zeros(nlines)
e_val = zeros(nbuses)
f_val = zeros(nbuses)

maxNSG = 1


# Define callback function for generating and adding cuts
W_val=zeros(nbuses) 
Wr_val=zeros(nlines) 
Wi_val=zeros(nlines)
function solveEta0Eigs()
	global nSG
	global η0Val
	global α_val, β_val, γ_val, δ_val
	global λF_val, μF_val, λT_val, μT_val
	global W_val, Wr_val, Wi_val
	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
	H=spzeros(2*nbuses,2*nbuses)
	for i in N
	  H[i,i] =  α_val[i] * acYshR[i] - β_val[i] * acYshI[i]  + δ_val[i] - γ_val[i]
	end
	for l in L
	  from = busIdx[lines[l].from];to = busIdx[lines[l].to] 
	  H[from,from] += λF_val[l] * acYffR[l] - μF_val[l] * acYffI[l] 
	  H[to,to] += λT_val[l] * acYttR[l] - μT_val[l] * acYttI[l] 
	end
	for i in N
	  H[nbuses+i,nbuses+i] = H[i,i]
	end
	for l in L
	  from = busIdx[lines[l].from];to = busIdx[lines[l].to] 
	  H[from,to] = 0.5*( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] )  
	  H[to,from] = H[from,to] 
	  H[nbuses+from, nbuses+to] = H[from,to]
	  H[nbuses+to, nbuses+from] = H[from,to] 
	  H[from, nbuses+to] = -0.5*( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] ) 
	  H[to,nbuses+from] = -H[from,nbuses+to]
	  H[nbuses+from, to] = H[to,nbuses+from] 
	  H[nbuses+to,from] = H[from,nbuses+to]
	end
	# Reference bus angle is zero
	e_val = zeros(nbuses)
	f_val = zeros(nbuses)
	E=eigs(H,nev=maxNSG,which=:SR, maxiter=3000, tol=1e-6)

	η0Val = E[1][1]
	nSG = 0
	for s=1:maxNSG
	 if E[1][s] <= -TOL
	  nSG += 1
	  for i in N
 	    e_val[i] = E[2][i,s]
	    f_val[i] = E[2][nbuses+i,s]
	  end
	  
	  for i in N
	    W_val[i] = e_val[i]^2 + f_val[i]^2
            sg_α[i,s] = acYshR[i] * W_val[i] 
            sg_β[i,s] = -acYshI[i] * W_val[i] 
	    sg_δ[i,s] = W_val[i] 
	    sg_γ[i,s] = -W_val[i]
	  end
	  for l in L
	    from = busIdx[lines[l].from]; to = busIdx[lines[l].to]
	    Wr_val[l] = e_val[from]*e_val[to] + f_val[from]*f_val[to]
	    Wi_val[l] = e_val[to]*f_val[from] - e_val[from]*f_val[to]
	    sg_λF[l,s] = (acYffR[l] * W_val[from] + acYftR[l] * Wr_val[l] + acYftI[l] * Wi_val[l])
	    sg_λT[l,s] = (acYttR[l] * W_val[to] + acYtfR[l] * Wr_val[l] - acYtfI[l] * Wi_val[l])
	    sg_μF[l,s] = (-acYffI[l] * W_val[from] - acYftI[l] * Wr_val[l] + acYftR[l]* Wi_val[l])
	    sg_μT[l,s] = (-acYttI[l] * W_val[to] - acYtfI[l] * Wr_val[l] - acYtfR[l] * Wi_val[l])
	  end
#=
	    sgnorm = sqrt(sum( sg_α[i,s]^2 + sg_β[i,s]^2 + sg_δ[i,s]^2 + sg_γ[i,s]^2 for i in N) + 
	    	+ sum( sg_λF[l,s]^2 + sg_λT[l,s]^2 + sg_μF[l,s]^2 + sg_μT[l,s]^2 for l in L))
	    if sgnorm > 0
	     for i in N 
              sg_α[i,s] /= sgnorm;  sg_β[i,s] /= sgnorm;  sg_δ[i,s] /= sgnorm; sg_γ[i,s] /= sgnorm
	     end
	     for l in L
	      sg_λF[l,s] /= sgnorm;  sg_λT[l,s] /= sgnorm;  sg_μF[l,s] /= sgnorm;  sg_μT[l,s] /= sgnorm
	     end
	    end
=#	
	else
	  break
	end
       end # s=1:maxNSG
end


function solveEta0SDP()
	global nSG, η0Val
	global W_val, Wr_val, Wi_val, e_val, f_val
	global α_val, β_val, γ_val, δ_val
	global λF_val, μF_val, λT_val, μT_val
	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
#The QP subproblem
  mSDP = Model(solver=IpoptSolver())
  #@variable(mSDP, -1 <= e[i=N] <= 1, start=1) #Add bounds later
  #@variable(mSDP, -1 <= f[i=N] <= 1, start=0)
  @variable(mSDP, e[i=N], start=1) #Add bounds later
  @variable(mSDP, f[i=N], start=0)
  @NLexpression(mSDP, exprW[i=N], e[i]^2 + f[i]^2)
  @NLexpression(mSDP, exprWR[l=L], e[busIdx[lines[l].from]]*e[busIdx[lines[l].to]] + f[busIdx[lines[l].from]]*f[busIdx[lines[l].to]])
  @NLexpression(mSDP, exprWI[l=L], e[busIdx[lines[l].to]]*f[busIdx[lines[l].from]] - e[busIdx[lines[l].from]]*f[busIdx[lines[l].to]])

  @NLconstraint(mSDP, vMagSumUB, sum(exprW[i] for i in N) <= 1) ### Trust-region constraint
  ### Does having a reference angle make sense in the context of generating cuts? 
  ###   Since the vectors e and f that are computed are normalized by construction?
  ###   For now, there is no reference angle.

	η0Val = 0

	for i in N
		setvalue(e[i], 1)
		setvalue(f[i], 0)
		e_val[i]=1
		f_val[i]=0
	end

	# Adjust QP subproblem
	@NLobjective(mSDP, Min,
	  sum( ( α_val[i] * acYshR[i] - β_val[i] * acYshI[i] + δ_val[i] - γ_val[i]) * exprW[i] for i in N)
	  + sum(
		  ( λF_val[l] * acYffR[l] - μF_val[l] * acYffI[l]) 						    * exprW[busIdx[lines[l].from]]
	        + ( 						λT_val[l] * acYttR[l] - μT_val[l] * acYttI[l] ) * exprW[busIdx[lines[l].to]]
	  	+ ( λF_val[l] * acYftR[l] - μF_val[l] * acYftI[l] + λT_val[l] * acYtfR[l] - μT_val[l] * acYtfI[l] ) * exprWR[l]
	  	+ ( λF_val[l] * acYftI[l] - λT_val[l] * acYtfI[l] + μF_val[l] * acYftR[l] - μT_val[l] * acYtfR[l] ) * exprWI[l]
		for l in L))
	status = solve(mSDP)
	if status == :Optimal || status == :UserLimit
	    nSG = 1
	    η0Val = getobjectivevalue(mSDP)
	    for i in N
 	        e_val[i] = getvalue(e[i])
	        f_val[i] = getvalue(f[i])
	    end
	    if(status == :UserLimit)
	       println("solveEta0SDP solve status $status") 
	    end
	else
	   println("solveEta0SDP: solve status $status") 
	   η0Val = 0
	   nSG = 0
	   return
	end
	for i in N
	    W_val[i] = e_val[i]^2 + f_val[i]^2
            sg_α[i,1] = acYshR[i] * W_val[i] 
            sg_β[i,1] = -acYshI[i] * W_val[i] 
	    sg_δ[i,1] = W_val[i] 
	    sg_γ[i,1] = -W_val[i]
	end
	for l in L
	    from = busIdx[lines[l].from]; to = busIdx[lines[l].to]
	    Wr_val[l] = e_val[from]*e_val[to] + f_val[from]*f_val[to]
	    Wi_val[l] = e_val[to]*f_val[from] - e_val[from]*f_val[to]
	    sg_λF[l,1] = (acYffR[l] * W_val[from] + acYftR[l] * Wr_val[l] + acYftI[l] * Wi_val[l])
	    sg_λT[l,1] = (acYttR[l] * W_val[to] + acYtfR[l] * Wr_val[l] - acYtfI[l] * Wi_val[l])
	    sg_μF[l,1] = (-acYffI[l] * W_val[from] - acYftI[l] * Wr_val[l] + acYftR[l]* Wi_val[l])
	    sg_μT[l,1] = (-acYttI[l] * W_val[to] - acYtfI[l] * Wr_val[l] - acYtfR[l] * Wi_val[l])
	end
end

