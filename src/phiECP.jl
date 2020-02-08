#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

#include("utils.jl")
include("../../MaximinOPF/src/MaximinOPF.jl")
using JuMP
using CPLEX 
using GLPK 
using Mosek, MosekTools
TOL=1e-4
"After the feasibility model is solved, we can query a subgradient of its value functions"
function computeFPSG(pm_orig, x_vals=Dict{Int64,Float64}() )
    optBD = JuMP.objective_value(pm_orig.model)
    sg=Dict{Int64,Float64}()
        pm = instantiate_model(pm_orig.data, SparseSDPWRMPowerModel, MaximinOPF.WRConicPost_PF)
        for l in ids(pm,pm.cnw,:branch)
          if !haskey(x_vals,l)
	    x_vals[l]=0
	    println("WARNING: x_vals not set in computeFPSG().")
          end
        end
        up_br1 = var(pm,pm.cnw,pm.ccnd,:up_br1)
        uq_br1 = var(pm,pm.cnw,pm.ccnd,:uq_br1)
        up_br0 = var(pm,pm.cnw,pm.ccnd,:up_br0)
        uq_br0 = var(pm,pm.cnw,pm.ccnd,:uq_br0)
        protected_arcs = filter(a->(a[1] in pm.data["protected_branches"]),ref(pm,pm.cnw,:arcs))
        inactive_arcs = filter(a->(a[1] in pm.data["inactive_branches"]),ref(pm,pm.cnw,:arcs))
        undecided_arcs = filter(a->!(a in protected_arcs || a in inactive_arcs),ref(pm,pm.cnw,:arcs))
        JuMP.@constraint(pm.model,  
	   sum( (1-x_vals[a[1]])*(up_br1[a,0] + uq_br1[a,0] + up_br1[a,1] + uq_br1[a,1]) 
	       + x_vals[a[1]]*(up_br0[a,0] + uq_br0[a,0] + up_br0[a,1] + uq_br0[a,1]) for a in undecided_arcs)
	 + sum( (up_br1[a,0] + uq_br1[a,0] + up_br1[a,1] + uq_br1[a,1])  for a in protected_arcs)
	 + sum( (up_br0[a,0] + uq_br0[a,0] + up_br0[a,1] + uq_br0[a,1])  for a in inactive_arcs) <= optBD+1e-6
        )
        JuMP.@variable(pm.model, slk[l in ids(pm,pm.cnw,:branch)] >= 0)
        JuMP.@variable(pm.model, xiUB[l in ids(pm,pm.cnw,:branch)] >= 0)
        JuMP.@variable(pm.model, xiLB[l in ids(pm,pm.cnw,:branch)] >= 0)
	for l in ids(pm,pm.cnw,:branch)
	  if x_vals[l] > 0.01
	    JuMP.set_upper_bound(xiLB[l],0) #FIX TO ZERO
	  end
	  if x_vals[l] < 0.99
	    JuMP.set_upper_bound(xiUB[l],0) #FIX TO ZERO
	  end
	end
	alpha = 0
	attack_budget_binding = false
	if abs(sum(x_vals[l] for l in ids(pm,:branch)) - pm.data["attacker_budget"]) < 1e-6
	  attack_budget_binding = true
	end
	if attack_budget_binding
	  alpha = JuMP.@variable(pm.model, alpha >= 0) ###For projecting the subgradient onto the norm of the attacker budget constraint when it is binding
	end
	JuMP.@expression(pm.model, phiSGExpr[l in ids(pm,pm.cnw,:branch)], 
	  sum( -(up_br1[a,0] + up_br1[a,1] + uq_br1[a,0] + uq_br1[a,1]) + (up_br0[a,0] + up_br0[a,1] + uq_br0[a,0] + uq_br0[a,1]) for a in filter(a->(a[1]==l),ref(pm,pm.cnw,:arcs))) )
	#JuMP.@expression(pm.model, sgExpr[l in ids(pm,pm.cnw,:branch)], phiSGExpr[l] + alpha + xiUB[l] - xiLB[l] )
	JuMP.@expression(pm.model, sgExpr[l in ids(pm,pm.cnw,:branch)], phiSGExpr[l] )
        JuMP.@constraint(pm.model, slkUBind[l in ids(pm,pm.cnw,:branch)], sgExpr[l] - slk[l] <= 0 )
        JuMP.@constraint(pm.model, slkLBind[l in ids(pm,pm.cnw,:branch)], -sgExpr[l] - slk[l] <= 0)
	JuMP.@variable(pm.model, psi)

	socvec = []
	push!(socvec,psi)
	for l in ids(pm, :branch)
	  push!(socvec,slk[l])
	end
	JuMP.@constraint(pm.model, socvec in SecondOrderCone())

	#JuMP.@constraint(pm.model, sum(slk[l] for l in ids(pm, :branch)) - psi <= 0)
	#JuMP.@constraint(pm.model, infNormSlk[l in ids(pm,pm.cnw,:branch)], slk[l] - psi <= 0)
	JuMP.@objective(pm.model,Min, psi)
        #objective_feasibility_problem(pm,x_vals)
	set_optimizer(pm.model,Mosek.Optimizer)
        set_parameters(pm.model,"MSK_IPAR_LOG"=>0)
	JuMP.optimize!(pm.model)
	if attack_budget_binding
	  alpha_val = JuMP.value(alpha)
	else
	  alpha_val = 0
	end

        for l in ids(pm, :branch)
          branch = ref(pm, pm.cnw, :branch, l)
          f_bus = branch["f_bus"]
          t_bus = branch["t_bus"]
          f_idx = (l, f_bus, t_bus)
          t_idx = (l, t_bus, f_bus)
#=
          upf0 = abs(  JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br0)[f_idx,0]) + JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br0)[f_idx,1]) )
          upt0 = abs(  JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br0)[t_idx,0]) + JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br0)[t_idx,1]) )
          uqf0 = abs(  JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br0)[f_idx,0]) + JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br0)[f_idx,1]) )
          uqt0 = abs(  JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br0)[t_idx,0]) + JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br0)[t_idx,1]) )
          upf1 = abs(  JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br1)[f_idx,0]) + JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br1)[f_idx,1]) )
          upt1 = abs(  JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br1)[t_idx,0]) + JuMP.value(var(pm, pm.cnw, pm.ccnd,:up_br1)[t_idx,1]) )
          uqf1 = abs(  JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br1)[f_idx,0]) + JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br1)[f_idx,1]) )
          uqt1 = abs(  JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br1)[t_idx,0]) + JuMP.value(var(pm, pm.cnw, pm.ccnd,:uq_br1)[t_idx,1]) )
=#
          sg[l] = JuMP.value(phiSGExpr[l])       
        end
    #println("SG: ",sg)
	sgNorm = 0.0
	for l in keys(sg)
	  sgNorm += sg[l]^2
	end
	sgNorm = sqrt(sgNorm)
	println("||sg||=",JuMP.value(psi)," it should be ", sgNorm)
    return sg, JuMP.value(psi)
end

function solvePhiECP(pm_data)

    #mMP = Model(with_optimizer(CPLEX.Optimizer))
    #mMP = direct_model(CPLEX.Optimizer())
    #MOI.set(mMP, MOI.RawParameter("CPX_PARAM_THREADS"), 1)
#=,
      CPX_PARAM_SCRIND=1,
      CPX_PARAM_TILIM=MAX_TIME,
      CPX_PARAM_MIPDISPLAY=4,
      CPX_PARAM_MIPINTERVAL=1,
      # CPX_PARAM_NODELIM=1,
      # CPX_PARAM_HEURFREQ=-1,
      CPX_PARAM_THREADS=1,
      CPX_PARAM_ADVIND=0))
=#

    #mMP = Model(with_optimizer(Cplex.Optimizer,CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_LPMETHOD=4,CPX_PARAM_SOLUTIONTYPE=2,CPX_PARAM_STARTALG=4))

    sg=Dict{Int64,Dict{Int64,Float64}}()
    x0=Dict{Int64,Dict{Int64,Float64}}()
    phi0=Dict{Int64,Float64}()
    x0[0] = Dict{Int64,Float64}()
    
    bestVal = -1e20
    bestSoln = Dict{Int64,Float64}()

    time_Start = time_ns()
    for ncuts=0:1000

      fp_pm = MaximinOPF.PF_FeasModel(pm_data, SparseSDPWRMPowerModel,x0[ncuts])
      pf_optmodel = fp_pm.model
      set_optimizer(pf_optmodel,Mosek.Optimizer)
      set_parameters(pf_optmodel,"MSK_IPAR_LOG"=>0)
      JuMP.optimize!(pf_optmodel)
      status=JuMP.termination_status(pf_optmodel)
      if status != OPTIMAL
        println("FLAGGING: Solve status=",status)
      end
      phi0[ncuts] = JuMP.objective_value(fp_pm.model)
      if phi0[ncuts] > bestVal
	bestVal=phi0[ncuts]
	bestSoln=copy(x0[ncuts])
      end
      sg[ncuts],sgNorm=computeFPSG(fp_pm, x0[ncuts])


    # The master problem MP
      mMP = direct_model(GLPK.Optimizer())
     # each x[l] is either 0 (line l active) or 1 (line l is cut)
      @variable(mMP, x[l in ids(fp_pm, :branch)], Bin)
      @variable(mMP, phi)
      @constraint(mMP, sum(x[l] for l in ids(fp_pm,:branch)) <= pm_data["attacker_budget"])
      @constraint(mMP, CP[kk=0:ncuts], phi0[kk] + sum( sg[kk][l]*(x[l]-x0[kk][l]) for l in ids(fp_pm, :branch))  - phi >= 0)
      @objective(mMP, Max, phi)
      JuMP.optimize!(mMP)
      status=JuMP.termination_status(mMP)
      if status != OPTIMAL
        println("FLAGGING: Solve status=",status)
      end
      mp_phival = JuMP.value(phi)
      x0[ncuts+1] = Dict{Int64,Float64}()
      for l in ids(fp_pm, :branch)
        x0[ncuts+1][l]=JuMP.value(x[l])
      end
      println("MP Optimal value of: ",mp_phival, " versus LB: ",phi0[ncuts], " versus best ",bestVal)
      println("MP soln has the following lines cut: ")
      for l in ids(fp_pm, :branch)
	if x0[ncuts+1][l] > 0.99
	  print(" ",l)
	end
      end
      print("\n")
      if mp_phival - bestVal < 1e-4  
	println("Terminating due to closure of UB-LB gap after ",ncuts+1," cuts.")
	break
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Finishing after ", time_End," seconds.")
    println("Best optimal value of: ",bestVal,".")
    println("Best soln has the following lines cut: ")
    for l in keys(bestSoln)
	if bestSoln[l] > 0.99
	  print(" ",l)
	end
    end
    

#=
  # CALLBACK FUNCTION FOR ADDING CUTS
    function addMPCutsLazyCB(cb)
println("#####################")
      x0 = Dict{Int64,Float64}()
      for l in ids(fp_pm,:branch)
        x0[l] = round(callback_value(cb,x[l]))
      end
      mp_phi_val = callback_value(cb,phi)
      MaximinOPF.reset_objective_feasibility_problem(fp_pm,x0)
      JuMP.optimize!(fp_pm.model)
      status=JuMP.termination_status(fp_pm.model)
      if status != OPTIMAL
        println("FLAGGING: Solve status=",status)
      end
      phi0 = JuMP.objective_value(fp_pm.model)
      sg=MaximinOPF.computeFPSG(fp_pm)
println("x0: ",x0)
println("sg: ",sg)
      # generate cut(s)
        if mp_phi_val - phi0 > 1e-4  
          con = @build_constraint( phi0 + sum( sg[l]*(x[l]-x0[l]) for l in ids(fp_pm, :branch)) - phi >= 0 )
	  MOI.submit(mMP, MOI.LazyConstraint(cb), con)
          ncuts += 1
        else
          println("Tolerance met for not generating a new lazy cut.")
        end
println("phi versus phi0: ",mp_phi_val, "   ", phi0, " and ncuts: ",ncuts)
println("*******************")
  end

  MOI.set(mMP, MOI.LazyConstraintCallback(), addMPCutsLazyCB)
  ncuts = 0   # collect number of cuts applied
  wtime = @elapsed JuMP.optimize!(mMP)
  status=JuMP.termination_status(mMP)
  @show wtime
  @printf("LAST ATTACK COMPUTED&& ")
  for l in L
    if round(JuMP.value(x[l])) == 1
        println("Cut line $l")
    end
  end
  @printf("\n")
  println("Optimal value of: ",JuMP.objective_value(mMP))
=#
end
