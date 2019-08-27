#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

MAX_TIME=24*3600

function solveDualDC(solninfo)
# The dual problem 

  mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=0))
  @variable(mMP, -rho_p <= α[i=N] <= rho_p)
  @variable(mMP, ζpUB[g=G] >=0)
  @variable(mMP, ζpLB[g=G] >=0)
  @variable(mMP, λF[l=L])
  @variable(mMP, λT[l=L])
  @variable(mMP, x[l=L], Bin)
  @constraint(mMP, sum(x[l] for l in L) <= K)
#=
  setlowerbound(x[5],1)
  setlowerbound(x[9],1)
  setlowerbound(x[10],1)
  setlowerbound(x[40],1)
=#

  
  for i in N
   for g in BusGeners[i]
    @constraint(mMP, 
	-α[i] + ζpUB[g] - ζpLB[g] == 0 )
   end
  end

   @variable(mMP, psiF[l=L]>=0)
   @variable(mMP, psiT[l=L]>=0)
   @expression(mMP, A[i=N], 0)
   for l in L
      from=busIdx[lines[l].from]; to=busIdx[lines[l].to]
      A[from] +=   psiF[l] - psiT[l] + λF[l]*dcYftI[l] - λT[l]*dcYtfI[l]
      A[to]   += (-psiF[l] + psiT[l] - λF[l]*dcYftI[l] + λT[l]*dcYtfI[l])
   end

   @constraint( mMP, dfTheta[i=N], A[i] == 0)


  @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] for g in G) 
     + sum( α[i]*(PD[i] + dcYshR[i])  for i in N) 
     + sum(λF[l]*(dcYffR[l] +  dcYftR[l]  ) + λT[l]*(dcYttR[l] +  dcYtfR[l] ) - (psiF[l] + psiT[l])*pi for l in L))

  @constraint(mMP, [l in L], α[busIdx[lines[l].from]] - rho_p*x[l] - λF[l] <= 0)
  @constraint(mMP, [l in L], α[busIdx[lines[l].from]] + rho_p*x[l] - λF[l] >= 0)
  @constraint(mMP, [l in L], -rho_p*(1 - x[l]) - λF[l] <= 0)
  @constraint(mMP, [l in L], rho_p*(1 - x[l]) - λF[l] >= 0)
  @constraint(mMP, [l in L], α[busIdx[lines[l].to]] - rho_p*x[l] - λT[l] <= 0)
  @constraint(mMP, [l in L], α[busIdx[lines[l].to]] + rho_p*x[l] - λT[l] >= 0)
  @constraint(mMP, [l in L], -rho_p*(1 - x[l]) - λT[l] <= 0)
  @constraint(mMP, [l in L], rho_p*(1 - x[l]) - λT[l] >= 0)

  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.

  if HEUR == 1 
  	@constraint(mMP, LambdaMuConstr1[l in L], λF[l]*YftI[l] - λT[l]*YtfI[l] == 0.0)
  elseif HEUR == 2
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l] == 0) 
  elseif HEUR == 3
	@constraint(mMP, LambdaMuConstr2[l in L], λF[l]*YtfR[l] - λT[l]*YftR[l]  == 0.0)
  end

  status=solve(mMP)
  if status == :Optimal || status == :Stall
      for l in L
        solninfo[l] = getvalue(x[l])
      end
      solninfo[nlines+1] = getobjectivevalue(mMP)
      solninfo[nlines+2] = getsolvetime(mMP)
      if status == :Stall
	println("solveNodeDC: Return status $status")
      end	
  else
	println("solveNodeDC: Return status $status")
  end
  return
end #end of function

function testDualDC()
  solninfo=zeros(nlines+2)
  solveDualDC(solninfo)
  dualObjval = solninfo[nlines+1]
  x_val = solninfo[1:nlines]
  primalObjval = solveFullModelDC(x_val)
  println("Primal value: ",primalObjval," and dual value: ",dualObjval)
end

