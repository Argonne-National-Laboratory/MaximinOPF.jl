#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

xiDVals = zeros(nlines)
xiRVals = zeros(nlines)
xiIVals = zeros(nlines)
xiNormVals = zeros(nlines)
nuVals = zeros(nlines)

# The master problem MP
  mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_THREADS=1))
  #mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=4))
  #mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=1))

 # Define the model here
  @variable(mMP, x[l=L], Bin)

  @variable(mMP, -rho_p <= α[i=N] <= rho_p)
  @variable(mMP, -rho_q <= β[i=N] <= rho_q)
  @variable(mMP, 0 <= γ[i=N] <= rhoU)
  @variable(mMP, 0 <= δ[i=N] <= rhoU)
  #@constraint(mMP, sum( γ[i] + δ[i] for i in N) <= rhoU)

  @variable(mMP, -rho_p <= λF[l=L] <= rho_p)
  @variable(mMP, -rho_p <= λT[l=L] <= rho_p)
  @variable(mMP, -rho_q <= μF[l=L] <= rho_q)
  @variable(mMP, -rho_q <= μT[l=L] <= rho_q)

  @variable(mMP, ζpUB[g=G] >= 0)
  @variable(mMP, ζpLB[g=G] >= 0)
  @variable(mMP, ζqUB[g=G] >= 0)
  @variable(mMP, ζqLB[g=G] >= 0)

  for i in N
   for g in BusGeners[i]
    @constraint(mMP, 
	-α[i] + ζpUB[g] - ζpLB[g] == 0 )
    @constraint(mMP, 
	-β[i] + ζqUB[g] - ζqLB[g] == 0 ) 
   end
  end

  #@variable(mMP, ζii[i=N] <= 1e6)

  #@objective(mMP, Max, sum(ζii[i] - Wmax[i]*δ[i] + Wmin[i]*γ[i] for i in N)) 
  @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γ[i]*Wmin[i]-δ[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N))

  @constraint(mMP, sum(x[l] for l in L) <= K)
  
  relFac = nlines #relaxation factor, needed to resolve a numerical problem noticeable in the lp relaxation
  @constraint(mMP, [l in L], α[busIdx[lines[l].from]] - relFac*rho_p*x[l] <= λF[l])
  @constraint(mMP, [l in L], α[busIdx[lines[l].from]] + relFac*rho_p*x[l] >= λF[l])
  @constraint(mMP, [l in L], -relFac*rho_p*(1 - x[l]) <= λF[l])
  @constraint(mMP, [l in L], relFac*rho_p*(1 - x[l]) >= λF[l])
  @constraint(mMP, [l in L], α[busIdx[lines[l].to]] - relFac*rho_p*x[l] <= λT[l])
  @constraint(mMP, [l in L], α[busIdx[lines[l].to]] + relFac*rho_p*x[l] >= λT[l])
  @constraint(mMP, [l in L], -relFac*rho_p*(1 - x[l]) <= λT[l])
  @constraint(mMP, [l in L], relFac*rho_p*(1 - x[l]) >= λT[l])

  @constraint(mMP, [l in L], β[busIdx[lines[l].from]] - relFac*rho_q*x[l] <= μF[l])
  @constraint(mMP, [l in L], β[busIdx[lines[l].from]] + relFac*rho_q*x[l] >= μF[l])
  @constraint(mMP, [l in L], -relFac*rho_q*(1 - x[l]) <= μF[l])
  @constraint(mMP, [l in L], relFac*rho_q*(1 - x[l]) >= μF[l])
  @constraint(mMP, [l in L], β[busIdx[lines[l].to]] - relFac*rho_q*x[l] <= μT[l])
  @constraint(mMP, [l in L], β[busIdx[lines[l].to]] + relFac*rho_q*x[l] >= μT[l])
  @constraint(mMP, [l in L], -relFac*rho_q*(1 - x[l]) <= μT[l])
  @constraint(mMP, [l in L], relFac*rho_q*(1 - x[l]) >= μT[l])

#=
  #The following four constraints are implied by solving the D1 subproblem with alphas, betas set to zero.
    #These constraints actually provide an exact characterization of the D1 dual function
    @constraint(mMP, EtaIInit1[i in N,(PD[i]-PminI[i]!=0 || QD[i]-QminI[i]!=0)], 
	ζii[i] <= (PD[i] - sum(Pmin[g] for g in BusGeners[i]))*α[i] + (QD[i] - sum(Qmin[g] for g in BusGeners[i]))*β[i])
    @constraint(mMP, EtaIInit2[i in N,(PD[i]-PminI[i]!=0 || QD[i]-QmaxI[i]!=0)], 
	ζii[i] <= (PD[i] - sum(Pmin[g] for g in BusGeners[i]))*α[i] + (QD[i] - sum(Qmax[g] for g in BusGeners[i]))*β[i])
    @constraint(mMP, EtaIInit3[i in N,(PD[i]-PmaxI[i]!=0 || QD[i]-QminI[i]!=0)], 
	ζii[i] <= (PD[i] - sum(Pmax[g] for g in BusGeners[i]))*α[i] + (QD[i] - sum(Qmin[g] for g in BusGeners[i]))*β[i])
    @constraint(mMP, EtaIInit4[i in N,(PD[i]-PmaxI[i]!=0 || QD[i]-QmaxI[i]!=0)], 
	ζii[i] <= (PD[i] - sum(Pmax[g] for g in BusGeners[i]))*α[i] + (QD[i] - sum(Qmax[g] for g in BusGeners[i]))*β[i])
 
    for i in N
      if( (PD[i]-PminI[i])==0 && (QD[i]-QminI[i])==0 && (PD[i]-PmaxI[i])==0 && (QD[i]-QmaxI[i])==0)
	setlowerbound(ζii[i],0)
	setupperbound(ζii[i],0)
      end
    end
=#


  #These constraints are not valid in general, but their inclusion often results in much faster time to near optimal solution.

  if HEUR == 1 
  	@constraint(mMP, LambdaMuConstr1[l in L], λF[l]*YftI[l] - λT[l]*YtfI[l] + μF[l]*YftR[l] - μT[l]*YtfR[l] == 0.0)
  elseif HEUR == 2
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l] == 0) 
	@constraint(mMP, muFequalsT[l in L], μF[l] - μT[l] == 0)
  elseif HEUR == 3
	@constraint(mMP, LambdaMuConstr2[l in L], λF[l]*YtfR[l] - λT[l]*YftR[l] - μF[l]*YtfI[l] + μT[l]*YftI[l] == 0.0)
  end

  println("Done with initial setup.")

# The dual SOCP constraints 
@expression(mMP, wCoeff[i=N], 
    α[i]*acYshR[i] - β[i]*acYshI[i] + δ[i] - γ[i]
    	+ sum( λF[l]*acYffR[l] - μF[l]*acYffI[l] for l in fromLines[i])   
	+ sum( λT[l]*acYttR[l] - μT[l]*acYttI[l] for l in toLines[i])  
	)
@expression(mMP, wRCoeff[l=L], λF[l]*acYftR[l] - μF[l]*acYftI[l] + λT[l]*acYtfR[l] - μT[l]*acYtfI[l])
@expression(mMP, wICoeff[l=L], λF[l]*acYftI[l] - λT[l]*acYtfI[l] + μF[l]*acYftR[l] - μT[l]*acYtfR[l])

@variable(mMP, xiD[l=L])
@variable(mMP, xiR[l=L])
@variable(mMP, xiI[l=L])
@variable(mMP, nu[l=L] >= 0)
@variable(mMP, WiiNonNeg[i=N] >= 0)

@constraint(mMP, dSOCP[i=N], 0.5*(sum(xiD[l]-nu[l] for l in fromLines[i]) - sum(xiD[l]+nu[l] for l in toLines[i])) - WiiNonNeg[i] + wCoeff[i]==0)
@constraint(mMP, dRSOCP[l=L], xiR[l] + wRCoeff[l] == 0)
@constraint(mMP, dISOCP[l=L], xiI[l] + wICoeff[l] == 0)

#@constraint(mMP, cSOCP[l=L], norm([xiD[l];xiR[l];xiI[l]]) <= nu[l])
function generateCuts(cb)
  global xiDVals, xiRVals, xiIVals, xiNormVals, nuVals, ncuts
  for l in L
    xiDVals[l]=getvalue(xiD[l])
    xiRVals[l]=getvalue(xiR[l])
    xiIVals[l]=getvalue(xiI[l])
    nuVals[l]=getvalue(nu[l])
    xiNormVals[l]=sqrt(xiDVals[l]^2+xiRVals[l]^2+xiIVals[l]^2)
    if xiNormVals[l] - nuVals[l] > 1e-6 && xiNormVals[l] > 1e-6
	    @lazyconstraint(cb, 
		(xiDVals[l]*xiD[l] + xiRVals[l]*xiR[l] + xiIVals[l]*xiI[l])/xiNormVals[l] - nu[l] <= 0)
      ncuts += 1

    end
  end
end
addlazycallback(mMP, generateCuts,fractional=false)


 wtime = @elapsed status = solve(mMP)
 @show wtime
 for l in L
    finalXSoln[l] = getvalue(x[l])
 end
 bestUBVal = getobjectivebound(mMP)
 nNodes = getnodecount(mMP)
 incVal = getobjectivevalue(mMP)
 finalXSoln[nlines+1] = getobjectivevalue(mMP)
 runtime = getsolvetime(mMP)


