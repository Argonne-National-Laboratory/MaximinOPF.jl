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
x_vals = zeros(nlines)
bestBd=-1

# The master problem MP
  mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=1,CPX_PARAM_THREADS=1))
  #mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=4))
  #mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=1))

 # Define the model here
  @variable(mMP, x[l=L], Bin)

  @variable(mMP, -rho_p <= α[i=N] <= rho_p)
  @variable(mMP, -rho_q <= β[i=N] <= rho_q)
  @variable(mMP, 0 <= γ[i=N] <= rhoU)
  @variable(mMP, 0 <= δ[i=N] <= rhoU)
  @constraint(mMP, [i=N], γ[i] + δ[i] <= rhoU)

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
  @variable(mMP, eta)

  @constraint(mMP, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] + ζqLB[g]*Qmin[g] - ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( γ[i]*Wmin[i]-δ[i]*Wmax[i] + α[i]*PD[i] + β[i]*QD[i] for i in N) - eta >= 0)
  @objective(mMP, Max, eta)

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

nCuts=0
function generateCuts(cb)
  global xiDVals, xiRVals, xiIVals, xiNormVals, nuVals, x_vals, nCuts, bestBd
  println("Starting generateCuts()")
  for l in L
    x_vals[l]=getvalue(x[l])
    if x_vals[l] < 0.5
      x_vals[l] = 0
    else
      x_vals[l] = 1
    end
    xiDVals[l]=getvalue(xiD[l])
    xiRVals[l]=getvalue(xiR[l])
    xiIVals[l]=getvalue(xiI[l])
    nuVals[l]=getvalue(nu[l])
    xiNormVals[l]=sqrt(xiDVals[l]^2+xiRVals[l]^2+xiIVals[l]^2)
  end
#=
  if bestBd < MathProgBase.cbgetobj(cb)
	bestBd = max(MathProgBase.cbgetobj(cb),0)
@show bestBd
	@lazyconstraint(cb, eta -1e-4 >= bestBd)
  end
=#
#@show x_vals
#@show xiNormVals
#@show nuVals

    #solve relaxed NLP with x fixed to x_vals
     nlpM = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=1))
     #nlpM = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_THREADS=1))
     @variable(nlpM, -rho_p <= r_α[i=N] <= rho_p)
     @variable(nlpM, -rho_q <= r_β[i=N] <= rho_q)
     @variable(nlpM, 0 <= r_γ[i=N] <= rhoU)
     @variable(nlpM, 0 <= r_δ[i=N] <= rhoU)
     @constraint(nlpM, [i=N], r_γ[i] + r_δ[i] <= rhoU)
     @variable(nlpM, r_ζpUB[g=G] >= 0)
     @variable(nlpM, r_ζpLB[g=G] >= 0)
     @variable(nlpM, r_ζqUB[g=G] >= 0)
     @variable(nlpM, r_ζqLB[g=G] >= 0)

     for i in N
       for g in BusGeners[i]
        @constraint(nlpM, 
	 -r_α[i] + r_ζpUB[g] - r_ζpLB[g] == 0 )
        @constraint(nlpM, 
	 -r_β[i] + r_ζqUB[g] - r_ζqLB[g] == 0 ) 
       end
     end
     @objective(nlpM, Max, sum(r_ζpLB[g]*Pmin[g] - r_ζpUB[g]*Pmax[g] + r_ζqLB[g]*Qmin[g] - r_ζqUB[g]*Qmax[g]  for g in G) 
	+ sum( r_γ[i]*Wmin[i]-r_δ[i]*Wmax[i] + r_α[i]*PD[i] + r_β[i]*QD[i] for i in N))


     fromIdx=zeros(Int,nlines);toIdx=zeros(Int,nlines)
     for l in L
      fromIdx[l]=busIdx[lines[l].from]; toIdx[l]=busIdx[lines[l].to]
     end
     for l in L 
      if x_vals[l] < 0.5
       if HEUR == 1 
  	@constraint(nlpM, r_α[fromIdx[l]]*YftI[l] - r_α[toIdx[l]]*YtfI[l] + r_β[fromIdx[l]]*YftR[l] - r_β[toIdx[l]]*YtfR[l] == 0.0)
       elseif HEUR == 2
	@constraint(nlpM, r_α[fromIdx[l]] - r_α[toIdx[l]] == 0) 
	@constraint(nlpM, r_β[fromIdx[l]] - r_β[toIdx[l]] == 0)
       elseif HEUR == 3
	@constraint(nlpM, r_α[fromIdx[l]]*YtfR[l] - r_α[toIdx[l]]*YftR[l] - r_β[fromIdx[l]]*YtfI[l] + r_β[toIdx[l]]*YftI[l] == 0.0)
       end
      end
     end

     # The dual SOCP constraints 
      @expression(nlpM, r_wCoeff[i=N], 
        r_α[i]*acYshR[i] - r_β[i]*acYshI[i] + r_δ[i] - r_γ[i]
    	+ sum( (1-x_vals[l])*(r_α[fromIdx[l]]*acYffR[l] - r_β[fromIdx[l]]*acYffI[l]) for l in fromLines[i])   
	+ sum( (1-x_vals[l])*(r_α[toIdx[l]]*acYttR[l] - r_β[toIdx[l]]*acYttI[l]) for l in toLines[i])  
	)
      @expression(nlpM, r_wRCoeff[l=L], (1-x_vals[l])*(r_α[fromIdx[l]]*acYftR[l] - r_β[fromIdx[l]]*acYftI[l] + r_α[toIdx[l]]*acYtfR[l] - r_β[toIdx[l]]*acYtfI[l]) )
      @expression(nlpM, r_wICoeff[l=L], (1-x_vals[l])*(r_α[fromIdx[l]]*acYftI[l] - r_α[toIdx[l]]*acYtfI[l] + r_β[fromIdx[l]]*acYftR[l] - r_β[toIdx[l]]*acYtfR[l])  )
      @variable(nlpM, r_xiD[l=L])
      @variable(nlpM, r_xiR[l=L])
      @variable(nlpM, r_xiI[l=L])
      @variable(nlpM, r_nu[l=L] >= 0)
      @variable(nlpM, r_WiiNonNeg[i=N] >= 0)

      @constraint(nlpM, dSOCP[i=N], 0.5*(sum(r_xiD[l]-r_nu[l] for l in fromLines[i]) - sum(r_xiD[l]+r_nu[l] for l in toLines[i])) - r_WiiNonNeg[i] + r_wCoeff[i]==0)
      @constraint(nlpM, dRSOCP[l=L], r_xiR[l] + r_wRCoeff[l] == 0)
      @constraint(nlpM, dISOCP[l=L], r_xiI[l] + r_wICoeff[l] == 0)

      @constraint(nlpM, cSOCP[l=L], norm([r_xiD[l];r_xiR[l];r_xiI[l]]) - r_nu[l] <= 0)

      status = solve(nlpM)
#@show status
#@show getobjectivevalue(nlpM)
#println("Testing to generate cuts")
   noCuts=true
   maxInfeas = 0
   for l in L
     maxInfeas = max(maxInfeas, xiNormVals[l] - nuVals[l])  
   end
@show maxInfeas
   if maxInfeas > 1e-6
    noCuts=false
    @lazyconstraint(cb,sum( x_vals[l]*(x_vals[l]-x[l]) + (1-x_vals[l])*(x[l]-x_vals[l]) for l in L) >= 1 )
   end
    for l in L
     if xiNormVals[l] - nuVals[l] > 1e-6 
      nCuts += 1
      xiDVals[l]=getvalue(r_xiD[l])
      xiRVals[l]=getvalue(r_xiR[l])
      xiIVals[l]=getvalue(r_xiI[l])
      nuVals[l]=getvalue(r_nu[l])
      xiNormVals[l]=sqrt(xiDVals[l]^2+xiRVals[l]^2+xiIVals[l]^2)
      @lazyconstraint(cb, 
		(xiDVals[l]*xiD[l] + xiRVals[l]*xiR[l] + xiIVals[l]*xiI[l])/xiNormVals[l] - nu[l] <= 0)
     end
    end
println("Ending generateCuts(): noCuts=",noCuts)
end

addlazycallback(mMP, generateCuts,fractional=false)


 wtime = @elapsed status = solve(mMP)
 @show status
 @show wtime
 for l in L
    finalXSoln[l] = getvalue(x[l])
 end
 bestUBVal = getobjectivebound(mMP)
 nNodes = getnodecount(mMP)
 incVal = getobjectivevalue(mMP)
 finalXSoln[nlines+1] = getobjectivevalue(mMP)
 runtime = getsolvetime(mMP)


#=
println("Testing for final feasibility.")
   for l in L
     xiDVals[l]=getvalue(xiD[l])
     xiRVals[l]=getvalue(xiR[l])
     xiIVals[l]=getvalue(xiI[l])
     xiNormVals[l]=sqrt(xiDVals[l]^2+xiRVals[l]^2+xiIVals[l]^2)
     nuVals[l]=getvalue(nu[l])
     @show xiNormVals[l] - nuVals[l]  
   end
=#
@show nNodes
@show nCuts
