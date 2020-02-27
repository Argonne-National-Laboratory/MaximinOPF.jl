include("../../MaximinOPF/src/MaximinOPF.jl")
using JuMP, MathOptInterface
using PowerModels
using Ipopt
using Mosek
using MosekTools
using CPLEX
#using SCIP
using LinearAlgebra
using Arpack
using Printf
PowerModels.silence()

function convertSOCtoPSD(model::JuMP.Model)
    #BRIDGE SOC CONSTRAINTS
    model_rsoc_moi = MOI.Utilities.Model{Float64}()
    rsoc_bridged_model = MOI.Bridges.Constraint.SOCtoPSD{Float64}(model_rsoc_moi)
    MOI.copy_to(rsoc_bridged_model,backend(model))

    model_psd_moi = MOI.Utilities.Model{Float64}()
    psd_bridged_model = MOI.Bridges.Constraint.RSOCtoPSD{Float64}(model_psd_moi)
    MOI.copy_to(psd_bridged_model,model_rsoc_moi)
    model_psd = JuMP.Model()
    MOI.copy_to(backend(model_psd),model_psd_moi)
    return model_psd
end

function solveMaxminECP(pm_data,pm_form,pm_optimizer,io=Base.stdout)
    #println(io,"Formulating and solving the form ",pm_form, " for problem ",pm_data["name"], " with attack budget K=",pm_data["attacker_budget"],".")
    MAX_N_COLS=100
    time_Start = time_ns()
    feas_pm = MaximinOPF.PF_FeasModel(pm_data, pm_form)
    #println(io,feas_pm.model)

    maxmin=Dict{String,Any}()
    maxmin["branch_ids"] = pm_data["undecided_branches"]
    println(maxmin["branch_ids"])
    maxmin["x_soln"]=Dict{Int64,Float64}()
    for l in maxmin["branch_ids"]
        maxmin["x_soln"][l] = 0
    end
    best_x_soln=copy(maxmin["x_soln"])
    bestLB=-1e20

    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; rm_rsoc=true, rm_therm_line_lim=true)
    maxmin["model"] = convertSOCtoPSD(base_maxmin)
    JuMP.set_optimizer(maxmin["model"],with_optimizer(CPLEX.Optimizer))
    JuMP.set_parameter(maxmin["model"],"CPXPARAM_ScreenOutput",0)


    feas=Dict{String,Any}()
    feas["branch_ids"] = pm_data["undecided_branches"]
    feas["x_soln"]=maxmin["x_soln"] #Shallow copy is OK here
    feas["opt_val"]=0

    ###Prepare to modify maxmin model, create feas model
    gatherPSDConInfo(maxmin)

    feas["model"],feas["ref_map"]=copy_model(maxmin["model"])
    JuMP.set_optimizer(feas["model"],pm_optimizer)
    gatherPSDConInfo(feas)
    removePSD_Constraints(maxmin) 
	add_initial_cuts(maxmin;io=io)
    println(io,maxmin["model"])
    println(io,maxmin["psd_con"])
#=
    relax_integrality(maxmin)
    JuMP.optimize!(maxmin["model"])
    println(io,"Solve status: ",JuMP.termination_status(maxmin["model"]))
    println(io,"Int-relaxed maxmin obj val:: ",JuMP.objective_value(maxmin["model"]))
    eval_cut_dual(maxmin; io=io)
    enforce_integrality(maxmin)
=#

    relax_integrality(feas)
	add_initial_cuts(feas;io=io)

    cut_via_dual_feas(feas,maxmin; add_cut=true, io=io)
    resolveMP(maxmin,feas; io=io)
    println(io,"Iteration 0: ")
    println("Iteration 0: ")
    print("\tUB: ",maxmin["UB_Val"]," with status: ",JuMP.termination_status(maxmin["model"])," with x solution: ")
    printX(maxmin["x_soln"])
    println("\tLB: ",feas["opt_val"], " with status: ",JuMP.termination_status(feas["model"]),", best value is: ",bestLB)
    #eval_cut_con(maxmin; io=io)

    #computeNewConstraintEig(maxmin;io=io)


  for ii=1:1
    cut_via_dual_feas(feas,maxmin; add_cut=true, io=io)
    if feas["opt_val"] > bestLB
      bestLB=feas["opt_val"]
      best_x_soln=copy(maxmin["x_soln"])
      println(io,"Updating incumbent solution: ")
      printX(best_x_soln,io)
      println(io,"The incumbent value is now: ",bestLB,".")
      println("Updating incumbent solution: ")
      printX(best_x_soln)
      println("The incumbent value is now: ",bestLB,".")
    end
    resolveMP(maxmin,feas; io=io)
    println(io,"Iteration $ii: ")
    println("Iteration $ii: ")
    print("\tUB: ",maxmin["UB_Val"]," with status: ",JuMP.termination_status(maxmin["model"])," with x solution: ")
      printX(maxmin["x_soln"])
    println("\tLB: ",feas["opt_val"], " with status: ",JuMP.termination_status(feas["model"]),", best value is: ",bestLB)
    #cut_refs[ii]=computeNewConstraintEig(maxmin;io=io)
    if maxmin["UB_Val"] - bestLB < 1e-3
	println("Terminating at iteration $ii since the relaxed solution is feasible.")
	print("Final attack solution: ")
	printX(best_x_soln)
	println("With value: ", bestLB)
	break
    end
    #delete_inactive_cuts(maxmin; io=io)

  end
end

function gatherPSDConInfo(model_dict)
    con_types=list_of_constraint_types(model_dict["model"])
    n_con_types=length(con_types)
    model_dict["psd_con"] = Dict{Int64,Dict{String,Any}}()
    nn = 0
    for cc=1:n_con_types
        if con_types[cc][2]==MathOptInterface.PositiveSemidefiniteConeTriangle
            psd_con = all_constraints(model_dict["model"], con_types[cc][1], con_types[cc][2]) 
            n_psd_con = length(psd_con)
            for kk=1:n_psd_con
                nn += 1
                model_dict["psd_con"][nn] = Dict{String,Any}()
                model_dict["psd_con"][nn]["ref"] = psd_con[kk]
                model_dict["psd_con"][nn]["expr"] = @expression(model_dict["model"], constraint_object(psd_con[kk]).func)
                model_dict["psd_con"][nn]["set"] = constraint_object(psd_con[kk]).set
                n_els = length(model_dict["psd_con"][nn]["expr"])
                model_dict["psd_con"][nn]["expr_val"] = zeros(n_els)
                model_dict["psd_con"][nn]["dual_val"] = zeros(n_els)
                model_dict["psd_con"][nn]["sg"] = Array{Float64,1}(undef, n_els)
            end
	    end
    end
end

function add_initial_cuts(model_dict; io=Base.stdout)
    model_dict["cuts"]=Dict{Int64,Any}()
    n_expr = length(model_dict["psd_con"])
    for nn=1:n_expr
	    diag_psd_vars = Dict{Int64,Any}()
        psd_expr=model_dict["psd_con"][nn]["expr"]
        n_els = length(psd_expr)
	    jj,ii=1,1
        for mm in 1:n_els
            if ii==jj
                diag_psd_vars[jj] = psd_expr[mm]
                model_dict["cuts"][ length(model_dict["cuts"]) ]=JuMP.@constraint(model_dict["model"],psd_expr[mm] >= 0)
                model_dict["cuts"][ length(model_dict["cuts"]) ]=JuMP.@constraint(model_dict["model"],psd_expr[mm] <= 1e3)
                jj += 1
                ii = 1
            else
                model_dict["cuts"][ length(model_dict["cuts"]) ]=JuMP.@constraint(model_dict["model"], diag_psd_vars[ii] + 2*psd_expr[mm] + psd_expr[mm + jj-ii] >= 0)
                model_dict["cuts"][ length(model_dict["cuts"]) ]=JuMP.@constraint(model_dict["model"], diag_psd_vars[ii] - 2*psd_expr[mm] + psd_expr[mm + jj-ii] >= 0)
                ii += 1
            end
        end
    end
    for kk in keys(model_dict["cuts"])
        println(io, model_dict["cuts"][kk])
    end
end

function removePSD_Constraints(model_dict) 
    n_expr = length(model_dict["psd_con"])
    for nn=1:n_expr
        delete(model_dict["model"],model_dict["psd_con"][nn]["ref"])
        delete!(model_dict["psd_con"][nn],"ref")
    end
end

function relax_integrality(model_dict)
    for l in model_dict["branch_ids"]
	    if JuMP.is_integer(variable_by_name(model_dict["model"],"x[$l]_1"))
            JuMP.unset_integer(variable_by_name(model_dict["model"],"x[$l]_1"))
            JuMP.set_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"),0)
            JuMP.set_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"),1)
	    else
            JuMP.set_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"),0)
            JuMP.set_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"),1)
	    end
    end
end

function enforce_integrality(model_dict)
    for l in model_dict["branch_ids"]
        if has_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
            delete_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
        end
        if has_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
            delete_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
        end
        JuMP.set_integer(variable_by_name(model_dict["model"],"x[$l]_1"))
    end
end

function convertVecToSG(Vec::Array{Float64,1},SG::Array{Float64,1})
	n_el=length(Vec)
	@assert n_el == length(SG)
	jj=1
        for nn=1:n_el
          if nn == (Int)(((jj+1)*jj)/2)
	    SG[nn] = Vec[nn]
	    jj += 1
          else
	    SG[nn] = 2*Vec[nn]
	  end
	end
end
function convertSGToVec(SG::Array{Float64,1},Vec::Array{Float64,1})
	n_el=length(Vec)
	@assert n_el == length(SG)
	jj=1
        for nn=1:n_el
          if nn == (Int)(((jj+1)*jj)/2)
	    Vec[nn] = SG[nn]
	    jj += 1
          else
	    Vec[nn] = 0.5*SG[nn]
	  end
	end
end

function cut_via_dual_feas(feas,maxmin; add_cut=true, fix_x=true, io=Base.stdout)
    n_con=length(feas["psd_con"])
    if fix_x
      for l in feas["branch_ids"]
        JuMP.fix(variable_by_name(feas["model"],"x[$l]_1"),maxmin["x_soln"][l];force=true)
      end
    end
    JuMP.optimize!(feas["model"])
    for nn=1:n_con
	feas["psd_con"][nn]["dual_val"][:] = JuMP.dual.(feas["psd_con"][nn]["ref"])[:]
    end
    
    feas["opt_val"] = JuMP.objective_value(feas["model"])
    println(io,"Optimal value of the dual feas problem is: ",feas["opt_val"], " with status: ",JuMP.termination_status(feas["model"]))

    if add_cut
      for nn=1:n_con
	    n_el=length(feas["psd_con"][nn]["dual_val"])
        convertVecToSG(feas["psd_con"][nn]["dual_val"],maxmin["psd_con"][nn]["sg"])
	    feas["cuts"][length(feas["cuts"])] = @constraint(feas["model"],   sum(maxmin["psd_con"][nn]["sg"][ii]*feas["psd_con"][nn]["expr"][ii] for ii=1:n_el)  >= 0)
        maxmin["cuts"][length(maxmin["cuts"])] = @constraint(maxmin["model"], sum(maxmin["psd_con"][nn]["sg"][ii]*maxmin["psd_con"][nn]["expr"][ii] for ii=1:n_el)  >= 0)
      end
    end

    if fix_x
      for l in feas["branch_ids"]
        JuMP.unfix(variable_by_name(feas["model"],"x[$l]_1"))
      end
    end
end

function eval_cut_con(model_dict; io=Base.stdout)
  println(io,"Evaluating cut values: ")
  for cc in keys(model_dict["cuts"])
      println(io,"$cc => ",trunc(JuMP.value(model_dict["cuts"][cc]);digits=5))
  end
end
function eval_cut_dual(model_dict; io=Base.stdout)
  println(io,"Evaluating cut dual values: ")
  for cc in keys(model_dict["cuts"])
      println(io,"$cc => ",trunc(JuMP.dual(model_dict["cuts"][cc]);digits=5))
  end
end

#=
function delete_inactive_cuts(maxmin; io=Base.stdout)
  n_del=0
  for ii in keys(maxmin["cut_str_ref"])
    for kk in keys(maxmin["cut_str_ref"][ii])
	if JuMP.value( maxmin["cut_str_ref"][ii][kk] ) > 1e-2
	  delete(maxmin["model"],maxmin["cut_str_ref"][ii][kk])
	  delete!(maxmin["model"],maxmin["cut_str_ref"][ii],kk)
	  n_del += 1
	end
    end
  end
  println(io,"Number of cuts deleted: ",n_del)
end
=#

###Inputs are symmetric matrices in triangular form
function triFrobIP(tri_mat1::Array{Float64,1},tri_mat2::Array{Float64,1})
   n_el=length(tri_mat1)
   @assert(n_el==length(tri_mat2))
   result = 0.0
   jj=1
   
   for nn=1:n_el
      if nn == (Int)(((jj+1)*jj)/2)
	result += tri_mat1[nn]*tri_mat2[nn]
	jj += 1
      else
	result += 2*tri_mat1[nn]*tri_mat2[nn]
      end
    end
    return result
end
###Input is symmetric matrix in triangular form
function computeEigenDecomp(tri_mat::Array{Float64,1})
   n_el=length(tri_mat)
   ncols=(Int)((-1+sqrt(1+8*n_el))/2)
   psd_mat = zeros(ncols,ncols)
   nn=1
   sq_tri_map=Dict{Tuple{Int,Int},Int}()
   tri_sq_map=Dict{Int,Tuple{Int,Int}}()
   for jj=1:ncols
      for ii=1:jj
	 tri_sq_map[nn]=(ii,jj)
	 sq_tri_map[(ii,jj)]=nn
	 psd_mat[ii,jj] = tri_mat[nn] 
	 psd_mat[jj,ii] = psd_mat[ii,jj]
	 nn += 1
      end
   end
   E = eigs(psd_mat, which=:SR, nev=(ncols-1))
   return E, sq_tri_map, tri_sq_map
end

function computeNewConstraintEig(maxmin;io=Base.stdout)
  n_cuts = 0
  cut_refs = Dict{Int64,Any}()
  for kk in keys(maxmin["psd_con"])
println(io,"maxmin status $kk: ",JuMP.termination_status(maxmin["model"]))
      PSD_Vec = maxmin["psd_con"][kk]["expr"]
      PSD0 = maxmin["psd_con"][kk]["expr_val"]
#println(io,"PSD0[$kk}: ",PSD0)
      E, sq_tri_map, tri_sq_map = computeEigenDecomp(PSD0)
      n_eigs = length(E[1])
println(io,"eigenval_$kk: ", E[1])
      n_el = length(maxmin["psd_con"][kk]["expr"])
   if E[1][1] < -1e-6
     SG = zeros(n_el)
     eig_avg = sum(-E[1][mm] for mm in filter(mmm->(E[1][mmm]) < 0, 1:n_eigs))
     for mm in filter(mmm->(E[1][mmm] <= 0), 1:n_eigs)
#=
        for nn=1:n_el
	  (ii,jj)=tri_sq_map[nn]
	  if ii==jj
	     SG[nn] -= (E[1][mm]/eig_avg)*E[2][ii,mm]*E[2][jj,mm]  
	  else
	     SG[nn] -= (E[1][mm]/eig_avg)*(E[2][ii,mm]*E[2][jj,mm] + E[2][jj,mm]*E[2][ii,mm]) 
	  end
        end
=#
        for nn=1:n_el
	  (ii,jj)=tri_sq_map[nn]
	  if ii==jj
	     SG[nn] = E[2][ii,mm]*E[2][jj,mm]  
	  else
	     SG[nn] = (E[2][ii,mm]*E[2][jj,mm] + E[2][jj,mm]*E[2][ii,mm]) 
	  end
        end
        cref=@constraint(maxmin["model"], sum(trunc(SG[nn];digits=10)*(PSD_Vec[nn]-PSD0[nn]) for nn=1:n_el) + E[1][mm]   >= 0)
        n_cuts += 1
        cut_refs[n_cuts]=cref
     end
     #cref=@constraint(maxmin["model"], sum(trunc(SG[nn];digits=10)*PSD_Vec[nn] for nn=1:n_el)  >= 0)
     #n_cuts += 1
     #cut_refs[n_cuts]=cref
   end
  end
  println("Eigs: Adding $n_cuts cuts.")
  return cut_refs
end


function resolveMP(maxmin,feas;io=Base.stdout)
  JuMP.optimize!(maxmin["model"])
  for kk in keys(maxmin["psd_con"])
    maxmin["psd_con"][kk]["expr_val"][:] = JuMP.value.(maxmin["psd_con"][kk]["expr"])[:]
  end
  maxmin["UB_Val"]=JuMP.objective_value(maxmin["model"])
  for l in keys(maxmin["x_soln"])
    x_var = variable_by_name(maxmin["model"],"x[$l]_1")
    x_val = JuMP.value(x_var)
    if x_val > 1.0-1.0e-8
	x_val = 1
    elseif x_val < 1.0e-8
	x_val = 0
    end
    maxmin["x_soln"][l]= x_val
  end
end

function printX(x_soln::Dict{Int64,Float64},io=Base.stdout)
  for l in keys(x_soln)
      if x_soln[l] > 0.01
	if x_soln[l] > 0.99
	 print(io," ",l)
	else
	 print(io," ($l,",trunc(x_soln[l];digits=3),")")
	end
      end
  end
  print(io,"\n")
end




testcase = Dict(
	"file" => "data/case30.m", 
 	"name" => "case30K4",  	
 	"attack_budget" => 4,
 	"inactive_indices" => [],
 	"protected_indices" => []
	)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

pm_optimizer=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=8)

io = open(string(testcase["name"],".out"), "w")
#io = Base.stdout
#sdp_relax=[SDPWRMPowerModel, SparseSDPWRMPowerModel]
#pm_form=SparseSDPWRMPowerModel
pm_form=SDPWRMPowerModel
solveMaxminECP(pm_data,pm_form,pm_optimizer,io)
close(io)
