include("../src/MaximinOPF.jl")
include("ProxSDP.jl")
using JuMP, MathOptInterface
using PowerModels
using Ipopt
using Mosek
using MosekTools
using CPLEX
#using SCIP
using LinearAlgebra
using SparseArrays
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
    MAX_N_ITER=50
    # "Mode definitions"
    PBM,ADMM=1,2
    #MODE=ADMM
    MODE=PBM

    time_Start = time_ns()

    maxmin=Dict{String,Any}()

    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false, rm_rsoc=true, rm_therm_line_lim=false)
    psd_base_maxmin = convertSOCtoPSD(base_maxmin)


    JuMP.set_optimizer(psd_base_maxmin,with_optimizer(Ipopt.Optimizer))
    #JuMP.set_optimizer(psd_base_maxmin,with_optimizer(CPLEX.Optimizer))
    #JuMP.set_parameter(psd_base_maxmin,"CPXPARAM_ScreenOutput",0)

    maxmin = prepare_to_solve_PSD_via_ProxPt( psd_base_maxmin; io=io )
    psd_expr = maxmin["model"][:psd_expr]

    #cplex_backend = CPLEX.Optimizer()
    #MOI.copy_to(cplex_backend,backend(psd_base_maxmin)) 
    #println(typeof(maxmin["model"]))
    #maxmin["model"] = direct_model(cplex_backend)

    maxmin["model"] = psd_base_maxmin
    maxmin["branch_ids"] = pm_data["undecided_branches"]
    maxmin["x_soln"]=Dict{Int64,Float64}()
    maxmin["x_soln_01"]=spzeros(Int64,maximum(maxmin["branch_ids"]))
    for l in maxmin["branch_ids"]
        maxmin["x_soln"][l] = 0
    end
    best_x_soln=copy(maxmin["x_soln"])
    bestLB=-1e20


    maxmin["cuts"]=Dict{SparseVector{Int64,Int64},Dict{Int64,Any}}()

    #enforce_integrality(maxmin["model"],maxmin["branch_ids"])
    fix_integer_vals(maxmin)
    if MODE==PBM
        solve_PSD_via_ProxPt(maxmin; max_n_iter=10, prox_t=0, io=io)
    elseif MODE==ADMM
        solve_PSD_via_ADMM(maxmin; max_n_iter=20, prox_t=1, io=io)
        total_sum_neg_eigs = 0
        for kk in keys(maxmin["psd_con"])
            PSD = maxmin["psd_con"][kk]
            total_sum_neg_eigs += PSD["neg_eigs_sum"]
            if PSD["neg_eigs_sum"] < 0
                JuMP.@constraint(maxmin["model"], sum( PSD["ip"][nn]*PSD["C"][nn]*psd_expr[kk,nn] for nn in 1:PSD["vec_len"]) >= 0 )
                #JuMP.@constraint(maxmin["model"], sum( PSD["ip"][nn]*PSD["sg"][nn]*psd_expr[kk,nn] for nn in 1:PSD["vec_len"]) >= 0 )
            end
        end
    end
    
    for ii=1:MAX_N_ITER
        resolveMP(maxmin; io=io)
        println(io,"Iteration $ii: ")
        println("Iteration $ii: ")
        print("\tUB: ",maxmin["opt_val"]," with status: ", maxmin["solve_status"]," with x solution: ")
        printX(maxmin["x_soln"])
        if PSDProjections(maxmin;io=io)
            total_sum_neg_eigs = 0
            fix_integer_vals(maxmin)
            if MODE==PBM
                solve_PSD_via_ProxPt(maxmin; max_n_iter=10, prox_t=1, io=io)
                total_sum_neg_eigs = 0
                for kk in keys(maxmin["psd_con"])
                    PSD = maxmin["psd_con"][kk]
                    if PSD["neg_eigs_sum"] < 0
                        total_sum_neg_eigs += PSD["neg_eigs_sum"]
                    end
                end
            elseif MODE==ADMM
                solve_PSD_via_ADMM(maxmin; max_n_iter=20,prox_t=1, io=io)
                total_sum_neg_eigs = 0
                for kk in keys(maxmin["psd_con"])
                    PSD = maxmin["psd_con"][kk]
                    total_sum_neg_eigs += PSD["neg_eigs_sum"]
                    if PSD["neg_eigs_sum"] < 0
                        JuMP.@constraint(maxmin["model"], sum( PSD["ip"][nn]*PSD["C"][nn]*psd_expr[kk,nn] for nn in 1:PSD["vec_len"]) >= 0 )
                        JuMP.@constraint(maxmin["model"], sum( PSD["ip"][nn]*PSD["sg"][nn]*psd_expr[kk,nn] for nn in 1:PSD["vec_len"]) >= 0 )
                    end
                end
            end
        end
        println("\tTotal sum of negative eigenvalues: ",total_sum_neg_eigs)
        if total_sum_neg_eigs >= -1e-6
            println("Terminating: feasible within tolerance, optimality verified.") 
            break
        end
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

function unfix_enforce_integrality(model_dict)
    for l in model_dict["branch_ids"]
        if has_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
            delete_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
        end
        if has_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
            delete_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"))
        end
        if is_fixed(variable_by_name(model_dict["model"],"x[$l]_1"))
            unfix(variable_by_name(model_dict["model"],"x[$l]_1"))
        end
        JuMP.set_integer(variable_by_name(model_dict["model"],"x[$l]_1"))
        set_lower_bound(variable_by_name(model_dict["model"],"x[$l]_1"),0)
        set_upper_bound(variable_by_name(model_dict["model"],"x[$l]_1"),1)
    end
end
function enforce_integrality(model,branch_ids)
    for l in branch_ids
        if has_lower_bound(variable_by_name(model,"x[$l]_1"))
            delete_lower_bound(variable_by_name(model,"x[$l]_1"))
        end
        if has_upper_bound(variable_by_name(model,"x[$l]_1"))
            delete_upper_bound(variable_by_name(model,"x[$l]_1"))
        end
        if is_fixed(variable_by_name(model,"x[$l]_1"))
            unfix(variable_by_name(model,"x[$l]_1"))
        end
        JuMP.set_integer(variable_by_name(model,"x[$l]_1"))
        set_lower_bound(variable_by_name(model,"x[$l]_1"),0)
        set_upper_bound(variable_by_name(model,"x[$l]_1"),1)
    end
end
function unfix_vars(model,branch_ids)
    for l in branch_ids
        if is_fixed(variable_by_name(model,"x[$l]_1"))
            unfix(variable_by_name(model,"x[$l]_1"))
        end
        set_lower_bound(variable_by_name(model,"x[$l]_1"),0)
        set_upper_bound(variable_by_name(model,"x[$l]_1"),1)
    end
end

function resolveMP(maxmin;io=Base.stdout)
    psd_expr = maxmin["model"][:psd_expr]
    JuMP.optimize!(maxmin["model"])
    for kk in keys(maxmin["psd_con"])
        n_el = maxmin["psd_con"][kk]["vec_len"]
        for nn=1:n_el
            maxmin["psd_con"][kk]["expr_val"][nn] = JuMP.value.(psd_expr[kk,nn])
        end
    end
    maxmin["opt_val"]=JuMP.objective_value(maxmin["model"])
    maxmin["solve_status"]=JuMP.termination_status(maxmin["model"])
    maxmin["x_soln"]=Dict{Int64,Float64}()
    for l in maxmin["branch_ids"]
        x_var = variable_by_name(maxmin["model"],"x[$l]_1")
        x_val = JuMP.value(x_var)
        if x_val > 1.0-1.0e-8
	        x_val = 1
            maxmin["x_soln_01"][l] = 1
        elseif x_val < 1.0e-8
	        x_val = 0
        end
        maxmin["x_soln"][l] = x_val
    end
end
#=
function resolveMP(maxmin;io=Base.stdout)
    #psd_expr_bd = maxmin["model"][:psd_expr_lbs]
    @objective(maxmin["model"], maxmin["obj_sense"], maxmin["model"][:linobj_expr])
    unfix_vars(maxmin["model"],maxmin["branch_ids"])
    mp_mip = JuMP.Model()
    JuMP.set_optimizer(mp_mip,with_optimizer(CPLEX.Optimizer))
    JuMP.set_parameter(mp_mip,"CPXPARAM_ScreenOutput",0)
    ref_index=MOI.copy_to(backend(mp_mip),backend(maxmin["model"]))
    enforce_integrality(mp_mip,maxmin["branch_ids"])
    #println("ref_index: ",ref_index[maxmin["model"][:psd_expr])
    #psd_expr = mp_mip[ref_index[:psd_expr]]
    # "1   CPXPROB_MILP"
    #CPLEX.set_prob_type!(maxmin["model"], 1) ### "Change to unfixed MIP"
    #JuMP.optimize!(maxmin["model"])
    JuMP.optimize!(mp_mip)
#=
    for kk in keys(maxmin["psd_con"])
        n_el = maxmin["psd_con"][kk]["vec_len"]
        for nn=1:n_el
            maxmin["psd_con"][kk]["expr_val"][nn] = JuMP.value(constraint_ref_with_index(mp_mip, ref_index[psd_expr_bd[kk,nn].index]))
        end
    end
=#

    maxmin["opt_val"]=JuMP.objective_value(mp_mip)
    maxmin["solve_status"]=JuMP.termination_status(mp_mip)
    maxmin["x_soln_01"]=spzeros(Int64,maximum(maxmin["branch_ids"]))
    maxmin["x_soln"]=Dict{Int64,Float64}()
    for l in maxmin["branch_ids"]
        x_var = variable_by_name(mp_mip,"x[$l]_1")
        x_val = JuMP.value(x_var)
        if x_val > 1.0-1.0e-8
	        x_val = 1
            maxmin["x_soln_01"][l] = 1
        elseif x_val < 1.0e-8
	        x_val = 0
        end
        maxmin["x_soln"][l] = x_val
    end

    # "8   CPXPROB_FIXEDMIQP"
    #CPLEX.set_prob_type!(maxmin["model"], 8) ### "Change to fixed quadratic MIP"
end
=#

function fix_integer_vals(model_info)
    relax_integrality(model_info)
    for l in model_info["branch_ids"]
        if is_fixed(variable_by_name(model_info["model"],"x[$l]_1"))
            unfix(variable_by_name(model_info["model"],"x[$l]_1"))
        end
        fix(variable_by_name(model_info["model"],"x[$l]_1"), model_info["x_soln"][l]; force=true)
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
pm_form=SparseSDPWRMPowerModel
#pm_form=SDPWRMPowerModel
solveMaxminECP(pm_data,pm_form,pm_optimizer,io)
close(io)
