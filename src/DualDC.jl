#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

MAX_TIME=24*3600
include("utils.jl")

function solveDualDC(opfdata,xsoln)
  # OBTAIN PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = 1:nbuses, 1:nlines, 1:ngens 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_DC
    Pmin,Pmax,PD=opfdata.Pmin,opfdata.Pmax,opfdata.PD
  rho_p=1;rho_q=1;rhoU=1
# The dual problem 
  mMP = Model(with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=4))
  @variable(mMP, -rho_p <= α[i=N] <= rho_p)
  @variable(mMP, ζpUB[g=G] >=0)
  @variable(mMP, ζpLB[g=G] >=0)
  @variable(mMP, λF[l=L])
  @variable(mMP, λT[l=L])
  @variable(mMP, x[l=L], Bin)
  @constraint(mMP, sum(x[l] for l in L) <= K)

  
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
      from=fromBus[l]; to=toBus[l]
      A[from] +=   psiF[l] - psiT[l] + λF[l]*Y["ftI"][l] - λT[l]*Y["tfI"][l]
      A[to]   += (-psiF[l] + psiT[l] - λF[l]*Y["ftI"][l] + λT[l]*Y["tfI"][l])
   end

   @constraint( mMP, dfTheta[i=N], A[i] == 0)


  @objective(mMP, Max, sum(ζpLB[g]*Pmin[g] - ζpUB[g]*Pmax[g] for g in G) 
     + sum( α[i]*(PD[i] + Y["shR"][i])  for i in N) 
     + sum(λF[l]*(Y["ffR"][l] +  Y["ftR"][l]  ) + λT[l]*(Y["ttR"][l] +  Y["tfR"][l] ) - (psiF[l] + psiT[l])*pi for l in L))

  @constraint(mMP, [l in L], α[fromBus[l]] - rho_p*x[l] - λF[l] <= 0)
  @constraint(mMP, [l in L], α[fromBus[l]] + rho_p*x[l] - λF[l] >= 0)
  @constraint(mMP, [l in L], -rho_p*(1 - x[l]) - λF[l] <= 0)
  @constraint(mMP, [l in L], rho_p*(1 - x[l]) - λF[l] >= 0)
  @constraint(mMP, [l in L], α[toBus[l]] - rho_p*x[l] - λT[l] <= 0)
  @constraint(mMP, [l in L], α[toBus[l]] + rho_p*x[l] - λT[l] >= 0)
  @constraint(mMP, [l in L], -rho_p*(1 - x[l]) - λT[l] <= 0)
  @constraint(mMP, [l in L], rho_p*(1 - x[l]) - λT[l] >= 0)

  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.

  if HEUR == 1 
  	@constraint(mMP, LambdaMuConstr1[l in L], λF[l]*Y["ftI"][l] - λT[l]*Y["tfI"][l] == 0.0)
  elseif HEUR == 2
	@constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l] == 0) 
  elseif HEUR == 3
	@constraint(mMP, LambdaMuConstr2[l in L], λF[l]*Y["tfR"][l] - λT[l]*Y["ftR"][l]  == 0.0)
  end

   wtime = @elapsed JuMP.optimize!(mMP)
   status=JuMP.termination_status(mMP)
   @show wtime
   for l in L
     xsoln[l] = JuMP.value(x[l])
   end
   bestUBVal = JuMP.objective_bound(mMP)
   nNodes = 0 #getnodecount(mMP)
   incVal = JuMP.objective_value(mMP)
   runtime = wtime #JuMP.solve_time(mMP)
   println("solveNodeDC: Return status $status")
   return bestUBVal,nNodes,incVal,runtime

end #end of function

function testDualDC(opfdata)
  xsoln=zeros(opfdata.nlines)
  optval,nNodes,incVal,runtime=solveDualDC(opfdata,xsoln)
  println("Optimal value: ", optval)
end

