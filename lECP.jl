#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

# The master problem MP
  mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50))
  #mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_LPMETHOD=4,CPX_PARAM_SOLUTIONTYPE=2,CPX_PARAM_STARTALG=4))
 # Define the model here
  @variable(mMP, x[l=L], Bin, start=0)
  @variable(mMP, -1 <= α[i=N] <= 1, start=0); @variable(mMP, -1 <= β[i=N] <= 1, start=0)
  @variable(mMP, δ[i=N] >= 0); @variable(mMP, γ[i=N] >= 0); @constraint(mMP, [i=N], δ[i]+γ[i] <= 1)

  @variable(mMP, λF[l=L], start=0); @variable(mMP, λT[l=L], start=0); @variable(mMP, μF[l=L], start=0); @variable(mMP, μT[l=L], start=0)

  @variable(mMP, ζpUB[g=G] >=0); @variable(mMP, ζpLB[g=G] >=0); @variable(mMP, ζqUB[g=G] >=0); @variable(mMP, ζqLB[g=G] >=0)

  for i in N
   for g in BusGeners[i]
    @constraint(mMP, -α[i] + ζpUB[g] - ζpLB[g] == 0 )
    @constraint(mMP, -β[i] + ζqUB[g] - ζqLB[g] == 0 ) 
   end
  end
#=
  if HEUR == 2
    @variable(mMP, xiPLB[l=L] >= 0); @variable(mMP, xiPUB[l=L] >= 0); @variable(mMP, xiQLB[l=L] >= 0); @variable(mMP, xiQUB[l=L] >= 0)
    hBd = 0.05
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γ[i]*Wmin[i]-δ[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N) - hBd*sum(xiPLB[l]+xiPUB[l]+xiQLB[l]+xiQUB[l]     for l in L)   )
  else
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γ[i]*Wmin[i]-δ[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N))
  end
=#
    @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γ[i]*Wmin[i]-δ[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N))

  @constraint(mMP, sum(x[l] for l in L) <= K)

 # McCormick inequalities enforcing bilinear equalities 
  @constraint(mMP, AMcf1[l in L], α[fromBus[l]] - x[l] <= λF[l]); @constraint(mMP, AMcf2[l in L], α[fromBus[l]] + x[l] >= λF[l])
  @constraint(mMP, AMcf3[l in L], -(1 - x[l]) <= λF[l]); @constraint(mMP, AMcf4[l in L], (1 - x[l]) >= λF[l])
  @constraint(mMP, AMct1[l in L], α[toBus[l]] - x[l] <= λT[l]); @constraint(mMP, AMct2[l in L], α[toBus[l]] + x[l] >= λT[l])
  @constraint(mMP, AMct3[l in L], -(1 - x[l]) <= λT[l]); @constraint(mMP, AMct4[l in L], (1 - x[l]) >= λT[l])

  @constraint(mMP, BMcf1[l in L], β[fromBus[l]] - x[l] <= μF[l]); @constraint(mMP, BMcf2[l in L], β[fromBus[l]] + x[l] >= μF[l])
  @constraint(mMP, BMcf3[l in L], -(1 - x[l]) <= μF[l]); @constraint(mMP, BMcf4[l in L], (1 - x[l]) >= μF[l])
  @constraint(mMP, BMct1[l in L], β[toBus[l]] - x[l] <= μT[l]); @constraint(mMP, BMct2[l in L], β[toBus[l]] + x[l] >= μT[l])
  @constraint(mMP, BMct3[l in L], -(1 - x[l]) <= μT[l]); @constraint(mMP, BMct4[l in L], (1 - x[l]) >= μT[l])

  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.

  if HEUR == 1 
  	@constraint(mMP, LambdaMuConstr1[l in L], λF[l]*YftI[l] - λT[l]*YtfI[l] + μF[l]*YftR[l] - μT[l]*YtfR[l] == 0.0)
  elseif HEUR == 2
#=
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l] + xiPUB[l] - xiPLB[l] == 0) 
	@constraint(mMP, muFequalsT[l in L], μF[l] - μT[l] + xiQUB[l] - xiQLB[l] == 0)
=#
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l]  == 0) 
	@constraint(mMP, muFequalsT[l in L], μF[l] - μT[l]  == 0)
  elseif HEUR == 3
	#@constraint(mMP, LambdaMuConstr2[l in L], λF[l]*YtfR[l] - λT[l]*YftR[l] - μF[l]*YtfI[l] + μT[l]*YftI[l] == 0.0)
	for l in L
	 if (l<34 || l>38) && l!=44 && l!=46 && l!=69 && l!=70
	  @constraint(mMP, λF[l] - λT[l]  == 0) 
	  @constraint(mMP, μF[l] - μT[l]  == 0)
	 end
	end
  end



sg_α = zeros(nbuses,3)
sg_β = zeros(nbuses,3)
sg_γ = zeros(nbuses,3) 
sg_δ = zeros(nbuses,3) 
sg_λF = zeros(nlines,3)
sg_μF = zeros(nlines,3)
sg_λT = zeros(nlines,3)
sg_μT = zeros(nlines,3)

nSG=0



maxNSG = 1
function generateCuts(cb)
  global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
  global η0Val, ncuts, nSG
  try
    #solveEta0Eigs()
    solveEta0SDP()
  catch exc
    println("Exception caught with eigs(), solving η0Val subproblem with Ipopt as recourse.") 
    println(exc)
    solveEta0SDP()
  end
  #if  η0Val <= -TOL 
  if nSG > 0 
	for s=1:nSG
	    @lazyconstraint(cb, 0.0 <= sum( (sg_α[i,s])* α[i] for i in N) + sum( (sg_β[i,s])* β[i] for i in N )  
		+ sum( (sg_γ[i,s])* γ[i] for i in N)  + sum( (sg_δ[i,s])* δ[i] for i in N)
		+ sum( (sg_λF[l,s])* λF[l] + (sg_λT[l,s])* λT[l] for l in L) + sum( (sg_μF[l,s])* μF[l] + (sg_μT[l,s])* μT[l] for l in L),
		localcut=useLocalCuts)
	    ncuts += 1
	end
  else
    println("Tolerance met for not generating a new lazy cut.")
  end
end

# Update Hessian
function updateHess(H)
        for i in N
	  α_val = getvalue(α[i]); β_val = getvalue(β[i]); γ_val = getvalue(γ[i]); δ_val = getvalue(δ[i]) 
          H[i,i] +=  α_val * acYshR[i] - β_val * acYshI[i]  + δ_val - γ_val
          H[nbuses+i,nbuses+i] += α_val * acYshR[i] - β_val * acYshI[i] + δ_val - γ_val
        end
        for l in L
          from = fromBus[l]; to = toBus[l] 
	  λF_val = getvalue(λF[l]); μF_val = getvalue(μF[l]); λT_val = getvalue(λT[l]); μT_val = getvalue(μT[l])
          H[from,from] += λF_val * acYffR[l] - μF_val * acYffI[l]
          H[nbuses+from,nbuses+from] += λF_val * acYffR[l] - μF_val * acYffI[l]
          H[to,to] += λT_val * acYttR[l] - μT_val * acYttI[l]
          H[nbuses+to,nbuses+to] += λT_val * acYttR[l] - μT_val * acYttI[l]
          H[from,to] += 0.5*( λF_val * acYftR[l] - μF_val * acYftI[l] + λT_val * acYtfR[l] - μT_val * acYtfI[l] )
          H[to,from] += 0.5*( λF_val * acYftR[l] - μF_val * acYftI[l] + λT_val * acYtfR[l] - μT_val * acYtfI[l] )
          H[nbuses+from, nbuses+to] += 0.5*( λF_val * acYftR[l] - μF_val * acYftI[l] + λT_val * acYtfR[l] - μT_val * acYtfI[l] )
          H[nbuses+to, nbuses+from] += 0.5*( λF_val * acYftR[l] - μF_val * acYftI[l] + λT_val * acYtfR[l] - μT_val * acYtfI[l] )
          H[to, nbuses+from] += 0.5*( λF_val * acYftI[l] - λT_val * acYtfI[l] + μF_val * acYftR[l] - μT_val * acYtfR[l] )
          H[nbuses+from, to] += 0.5*( λF_val * acYftI[l] - λT_val * acYtfI[l] + μF_val * acYftR[l] - μT_val * acYtfR[l] )
          H[from,nbuses+to] -= 0.5*( λF_val * acYftI[l] - λT_val * acYtfI[l] + μF_val * acYftR[l] - μT_val * acYtfR[l] )
          H[nbuses+to,from] -= 0.5*( λF_val * acYftI[l] - λT_val * acYtfI[l] + μF_val * acYftR[l] - μT_val * acYtfR[l] )
        end

end

# Define callback function for generating and adding cuts
function solveEta0Eigs()

	global nSG, maxNSG
	global η0Val
	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
	H=spzeros(2*nbuses,2*nbuses)
	updateHess(H)
	# Reference bus angle is zero
	E=eigs(H,nev=6,which=:SR, maxiter=100000, tol=1e-8)
#println("(",E[1][1],",",E[1][2],",",E[1][3],",",E[1][4],",",E[1][5],",",E[1][6],")")
	η0Val = E[1][1]
	nSG = 0
        for s=1:maxNSG
	 sInd = 2*(s-1)+1
	 if E[1][sInd] <= -TOL 
	  nSG += 1
	  for i in N
 	    e_val = E[2][i,sInd]; f_val = E[2][nbuses+i,sInd]; W_val = e_val^2 + f_val^2
            sg_α[i,s] = acYshR[i] * W_val; sg_β[i,s] = -acYshI[i] * W_val 
	    sg_δ[i,s] = W_val; sg_γ[i,s] = -W_val
	  end
	  for l in L
	    from = fromBus[l]; to = toBus[l]
 	    e_valF = E[2][from,sInd]; f_valF = E[2][nbuses+from,sInd]; W_valF = e_valF^2 + f_valF^2
 	    e_valT = E[2][to,sInd]; f_valT = E[2][nbuses+to,sInd]; W_valT = e_valT^2 + f_valT^2
	    Wr_val = e_valF*e_valT + f_valF*f_valT; Wi_val = e_valT*f_valF - e_valF*f_valT
	    sg_λF[l,s] = (acYffR[l] * W_valF + acYftR[l] * Wr_val + acYftI[l] * Wi_val)
	    sg_λT[l,s] = (acYttR[l] * W_valT + acYtfR[l] * Wr_val - acYtfI[l] * Wi_val)
	    sg_μF[l,s] = (-acYffI[l] * W_valF - acYftI[l] * Wr_val + acYftR[l] * Wi_val)
	    sg_μT[l,s] = (-acYttI[l] * W_valT - acYtfI[l] * Wr_val - acYtfR[l] * Wi_val)
	  end
	end
    end #s=1:3
#print("\n")
end

function solveEta0SDP()
	global nSG, η0Val
	global W_val, Wr_val, Wi_val, e_val, f_val
	global α_val, β_val, γ_val, δ_val
	global λF_val, μF_val, λT_val, μT_val
	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT
#The QP subproblem
  mSDP = Model(solver=IpoptSolver())
  #@variable(mSDP, -1 <= e[i=N] <= 1, start=0) #Add bounds later
  #@variable(mSDP, -1 <= f[i=N] <= 1, start=0)
  @variable(mSDP, e[i=N], start=0) #Add bounds later
  @variable(mSDP, f[i=N], start=0)
  @NLexpression(mSDP, exprW[i=N], e[i]^2 + f[i]^2)
  @NLexpression(mSDP, exprWR[l=L], e[fromBus[l]]*e[toBus[l]] + f[fromBus[l]]*f[toBus[l]])
  @NLexpression(mSDP, exprWI[l=L], e[toBus[l]]*f[fromBus[l]] - e[fromBus[l]]*f[toBus[l]])

  @NLconstraint(mSDP, vMagSumUB, sum(exprW[i] for i in N) <= 1) ### Trust-region constraint
  ### Does having a reference angle make sense in the context of generating cuts? 
  ###   Since the vectors e and f that are computed are normalized by construction?
  ###   For now, there is no reference angle.

	η0Val = 0

	for i in N
	  setvalue(e[i], 1); setvalue(f[i], 0)
	end

	H=spzeros(2*nbuses,2*nbuses)
	updateHess(H)

	# Adjust QP subproblem
	@NLobjective(mSDP, Min, sum( H[i,i]*(e[i]^2+f[i]^2) for i in N) 
			+ 2*sum( H[fromBus[l],toBus[l]]*(e[fromBus[l]]*e[toBus[l]]+f[fromBus[l]]*f[toBus[l]])   for l in L)
			- 2*sum( H[fromBus[l],nbuses+toBus[l]]*(f[fromBus[l]]*e[toBus[l]]-e[fromBus[l]]*f[toBus[l]])   for l in L)
	)
	status = solve(mSDP)
	if status == :Optimal || status == :UserLimit
	    nSG = 1
	    η0Val = getobjectivevalue(mSDP)
	    for i in N
	        W_val = getvalue(e[i])^2 + getvalue(f[i])^2
                sg_α[i,1] = acYshR[i] * W_val 
                sg_β[i,1] = -acYshI[i] * W_val 
	        sg_δ[i,1] = W_val 
	        sg_γ[i,1] = -W_val
	    end
	    for l in L
	       from = fromBus[l]; to = toBus[l]
	       W_valF = getvalue(e[from])^2 + getvalue(f[from])^2
	       W_valT = getvalue(e[to])^2 + getvalue(f[to])^2
	       Wr_val = getvalue(e[from])*getvalue(e[to]) + getvalue(f[from])*getvalue(f[to])
	       Wi_val = getvalue(e[to])*getvalue(f[from]) - getvalue(e[from])*getvalue(f[to])
	       sg_λF[l,1] = (acYffR[l] * W_valF + acYftR[l] * Wr_val + acYftI[l] * Wi_val)
	       sg_λT[l,1] = (acYttR[l] * W_valT + acYtfR[l] * Wr_val - acYtfI[l] * Wi_val)
	       sg_μF[l,1] = (-acYffI[l] * W_valF - acYftI[l] * Wr_val + acYftR[l]* Wi_val)
	       sg_μT[l,1] = (-acYttI[l] * W_valT - acYtfI[l] * Wr_val - acYtfR[l] * Wi_val)
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
end

function addMPCutsLazyCB(cb)
	global time_Eta0SP, time_MP, mpstart, avg_TimeMP, total_TimeMP, nMPUpdates, numcalls_Eta0SP, ncuts, attacks, noAttacks, bestAttack, noNoGoodCuts
  	global η0Val
  	global W_val, Wr_val, Wi_val
  	global α_val, β_val, γ_val, δ_val 
  	global λF_val, μF_val, λT_val, μT_val, x_val
  	global sg_α, sg_β, sg_γ, sg_δ, sg_λF, sg_μF, sg_λT, sg_μT

	time_MP= (time_ns() - mpstart)/1e9
	total_TimeMP += time_MP
	#nMPUpdates += 1
	#avg_TimeMP = total_TimeMP/nMPUpdates
	# get the variable values
	for i in N
		α_val[i] = getvalue(α[i])
		β_val[i] = getvalue(β[i])
		γ_val[i] = getvalue(γ[i])
		δ_val[i] = getvalue(δ[i]) 
	end
	for l in L
    	   x_val[l] = round(getvalue(x[l]))
		  λF_val[l] = getvalue(λF[l])
		  μF_val[l] = getvalue(μF[l]) 
		  λT_val[l] = getvalue(λT[l])
		  μT_val[l] = getvalue(μT[l]) 
	end

	# generate cuts
	time_Eta0Start = time_ns()
	#  generate cut(s)
	generateCuts(cb) ### the implementation depends on which eta0 SP code is included above
	time_Eta0SP += (time_ns()-time_Eta0Start)/1e9
	numcalls_Eta0SP += 1
	#mpstart = time_ns() #For updating MP
end

addlazycallback(mMP, addMPCutsLazyCB)
mpRootStartTime = time_ns() #For solving root node
mpstart = time_ns() #For updating MP
wtime = @elapsed status = solve(mMP)
@show wtime

@printf("\n********************FINAL RESULT FOR CASE %s WITH %d LINES CUT H%dR%d*******************\n",CASE_NUM,K, HEUR,FORM)

@printf("LAST ATTACK && ")
for l in L
	if round(getvalue(x[l])) == 1
		println("Cut line $l")
	end
end
@printf("\n")
#=
for l in L
    from = fromBus[l];to = toBus[l] 
    if getvalue(x[l]) < 0.5
	println("Line $l: (", getvalue(α[from]),",",getvalue(α[to]),")  (",getvalue(β[from]),",",getvalue(β[to]),")")
    end
end
=#
@printf("MPs took: avg: %f, total %f.\n",avg_TimeMP, total_TimeMP)

###Printing summary as one line to aid in latex writeup
@printf("Summary\t")
if(HEUR == 0 && FORM == AC)
    @printf("%d  &  ",K)
else
    @printf("    &   ")
end
if(FORM == AC )
    @printf("\$\\HSet^{%d}\$  &  ",HEUR)
else
    @printf("    &   ")
end

if FORM == AC
    @printf("AC &    ")
elseif FORM == SDP  
    @printf("SDP &    ")
elseif FORM == SOCP || FORM == SOCPBds 
    @printf("SOCP &    ")
elseif FORM == DC 
    @printf("DC &    ")
end

tShift = wtime - getsolvetime(mMP)

@printf("%d  & ", ncuts)
@printf("%.2f  & ", 100*time_Eta0SP/getsolvetime(mMP))
@printf("%d   &  ", getnodecount(mMP))
@printf("%.2f   &  ", getsolvetime(mMP)/getnodecount(mMP))
@printf("%d   &  ", getsolvetime(mMP))
#@printf("%.3f   &  ", getobjectivebound(mMP))

#@printf("%d    &   ", round(time_root-tShift))
#@printf("%d   &  ", round(time_Eta0SP))
#@printf("%d   &  ", round(total_TimeMP))
  for l in L
   if getvalue(x[l]) > 0.5
        @printf(" %d",l)
   end
  end
  #@printf("\t&\t%.3f \t&\t%.3f\t & \t%.3f", bestAttack[nlines+1+SDP], bestAttack[nlines+1+SOCP],bestAttack[nlines+1+DC])
@printf("  & %.3f ", getobjectivevalue(mMP))
if abs(getobjectivebound(mMP)-getobjectivevalue(mMP)) < 1e-3
  @printf("&  0.0\\%% ")
elseif incVal > 1e-3
  percentGap= 100*(getobjectivebound(mMP)-getobjectivevalue(mMP))/getobjectivevalue(mMP)
  @printf("&  %.1f\\%% ",percentGap)
else
  print("& ---     ")
end
@printf(" \\\\ \n")

