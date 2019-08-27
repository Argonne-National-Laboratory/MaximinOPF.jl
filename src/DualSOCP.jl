#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

function solveDSOCP(xSoln)
# The master problem MP
  #mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_THREADS=1))
  #mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=4))
  mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=10,MSK_IPAR_NUM_THREADS=1))

 # Define the model here
  @variable(mMP, x[l=L], Bin)
  rho_p=1;rho_q=1;rhoU=1

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
@constraint(mMP, cSOCP[l=L], norm([xiD[l];xiR[l];xiI[l]]) <= nu[l])

 wtime = @elapsed status = solve(mMP)
 @show wtime
 for l in L
  xSoln[l] = getvalue(x[l])
 end
 bestUBVal = getobjectivebound(mMP)
 nNodes = 0 #getnodecount(mMP)
 incVal = getobjectivevalue(mMP)
 runtime = getsolvetime(mMP)
 return bestUBVal,nNodes,incVal,runtime
end #end of function


finalXSoln = zeros(Int,nlines)
bestUBVal,nNodes,incVal,runtime = solveDSOCP(finalXSoln)

@printf("\n********************FINAL RESULT FOR CASE %s WITH %d LINES CUT H%dR%d*******************\n",CASE_NUM,K, HEUR,FORM)

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

printResults=true
if printResults
if FORM == AC
    @printf("AC &    ")
elseif FORM == SDP  
    @printf("SDP &    ")
elseif FORM == SOCP || FORM == SOCPBds 
    @printf("SOCP &    ")
elseif FORM == DC 
    @printf("DC &    ")
elseif FORM == ACSOC 
    @printf("ACSOC &    ")
else
    println("Not implemented")
    printResults=false
end

#@printf("%d & ", ncuts)
#@printf("%.2f  & ", 100*time_Eta0SP/runtime)
@printf("%d &  ", nNodes)
@printf("%.2f &  ", runtime/nNodes)
#@printf("%d & ", ncuts)
# @printf("%d &   ", round(time_root-tShift))
#@printf("%d &  ", round(time_Eta0SP))
#@printf("%d &  ", round(total_TimeMP))
@printf("%d &  ", runtime)
printX(finalXSoln)
@printf(" & %.3f &  ", incVal)
if abs(bestUBVal-incVal) < 1e-3
  @printf("0.0\\%% &")
elseif incVal > 1e-3
  percentGap= 100*(bestUBVal-incVal)/incVal
  @printf("%.1f\\%% &",percentGap)
else
  print("---    & ")
end

  #@printf("\t&\t%.3f \t&\t%.3f\t & \t%.3f", optSDP, optSOCP, optDC)
@printf(" \\\\ \n")

println("No. Nodes: ", nNodes)
println("Best bound:  ", bestUBVal)
@printf("Objective value: %.3f\n", incVal)
@show runtime
@show finalXSoln 
end

