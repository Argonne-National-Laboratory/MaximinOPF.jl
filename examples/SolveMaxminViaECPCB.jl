#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("../../MaximinOPF/src/MaximinOPF.jl")
include("../../MaximinOPF/src/utils.jl")
include("psd_utils.jl")
using JuMP, MathOptInterface
using Mosek, MosekTools
using CPLEX
using GLPK
using JuMP
using PowerModels
PowerModels.silence()




function solveMaxminViaCB(pm_data,pm_form; use_dual_minmax=true)
    #global MAX_TIME
    K=pm_data["attacker_budget"]
    art_bd=ceil(getMaxShedding(pm_data);digits=6)

    start_time = time_ns()

    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false, rm_rsoc=true, rm_therm_line_lim=true)
    global branch_ids=sort(collect(pm_data["undecided_branches"]))
    println("Branch_ids: ",branch_ids)

    cp_base_maxmin = convertSOCtoPSD(base_maxmin)
    JuMP.set_optimizer(cp_base_maxmin,with_optimizer(CPLEX.Optimizer))
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_ScreenOutput",0)
    #JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_Simplex_Tolerances_Optimality",1e-8)
	JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_Parallel", 1)
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_Threads", 1)
    global MIP_INT_TOL=1e-5
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_MIP_Tolerances_Integrality", MIP_INT_TOL)
    #JuMP.set_optimizer(cp_base_maxmin,with_optimizer(GLPK.Optimizer))

    global cp_model_info=Dict{String,Any}("branch_ids"=>branch_ids,"x_soln"=>Dict{Int64,Float64}(),
            "x_soln_str"=>"","heur_x_soln"=>Dict{Int64,Float64}(),"attacker_budget"=>K)
    cp_model_info["model"] = cp_base_maxmin
    ### "Sets `model_info[\"cp_model\"] = cp_base_maxmin` and also the 'psd_info' entry"
    gatherPSDConInfo(cp_model_info) ### "Sets the 'psd_info' entry of model_info"
    removePSD_Constraints(cp_model_info["psd_info"])
    #println(cp_model_info["model"])
    #bound_obj(cp_model_info; bd_mag=art_bd)
    println("Artifical bound on the objective function: ",art_bd)
	add_psd_initial_cuts(cp_model_info; bdmag=1e3,io=devnull)
    
    global feasXs = Dict{String,Dict{String,Any}}()
    feasXs[""]=Dict("x_soln"=>Dict{Int64,Float64}(),"x_soln_str"=>"", "bound_value"=>0, "value"=>0, "n_visits"=>1) ### CREATE INITIAL INCUMBENT SOLN

    for l in branch_ids
        cp_model_info["x_soln"][l] = 0
    end

#=
    n_init_cuts=100
    psd_expr = cp_model_info["model"][:psd_expr]
    PSD=cp_model_info["psd_info"]
    for cc=1:n_init_cuts
        solveSP(cp_model_info; fix_x=false, compute_projection=true, compute_psd_dual=false)
        if cp_model_info["orth_norm"] >= 1e-4
            for kk in keys(PSD)
                JuMP.@constraint(cp_model_info["model"],
                    #(1/PSD[kk]["orth_norm"])*sum( PSD[kk]["ip"][nn]*PSD[kk]["orth_expr_val"][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) <= 0 
                    sum( PSD[kk]["ip"][nn]*PSD[kk]["sg"][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) >= 0 
                )
            end
        end
        #add_cuts(cp_model_info, cp_model_info["psd_info"]; cut_type="sg")
        println("Init cut $cc root node bound value: ",cp_model_info["opt_val"], " infeas: ",cp_model_info["orth_norm"])
    end
=#

    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_ScreenOutput",1)
    global n_sdp_solves = 1
    feasXs[""]["value"] = 0.0
    global IncX = feasXs[""]
    global nXs = 1

    global x_vars = Dict{Int64,Any}()
    for l in branch_ids
        x_vars[l] = JuMP.variable_by_name(cp_model_info["model"],"x[$l]_1")
        if has_lower_bound(x_vars[l])
            delete_lower_bound(x_vars[l])
        end
        if has_upper_bound(x_vars[l])
            delete_upper_bound(x_vars[l])
        end
        JuMP.set_integer(x_vars[l])
        set_lower_bound(x_vars[l],0)
        set_upper_bound(x_vars[l],1)
    end

    function my_callback_function(cb_data)
        global branch_ids, feasXs,nXs,IncX,n_sdp_solves
        global MIP_INT_TOL
        global x_vars
        global cp_model_info
        psd_expr = cp_model_info["model"][:psd_expr]
        PSD = cp_model_info["psd_info"]
        x_is_int=true
        tol = 1e-5
        x_soln_str=""
        #println("callback at x_soln:")
        obj_expr=JuMP.objective_function(cp_model_info["model"], AffExpr)
        bd_val = 0.0
        for tt in keys(obj_expr.terms)
            bd_val += obj_expr.terms[tt]*JuMP.callback_value(cb_data, tt)
        end
        for l in branch_ids
            cp_model_info["x_soln"][l] = JuMP.callback_value(cb_data, x_vars[l])
            if cp_model_info["x_soln"][l] <= MIP_INT_TOL
                cp_model_info["x_soln"][l] = 0
            elseif cp_model_info["x_soln"][l] >= 1-MIP_INT_TOL
                cp_model_info["x_soln"][l] = 1
                x_soln_str = string(x_soln_str," $l")
            else
                x_is_int=false
            end
        end
        #println("\n")
        if x_is_int
            for kk in keys(PSD)
                for nn=1:PSD[kk]["vec_len"]
                    PSD[kk]["expr_val"][nn] = 0
                    if typeof(psd_expr[kk,nn])==VariableRef
                        PSD[kk]["expr_val"][nn] = JuMP.callback_value(cb_data,psd_expr[kk,nn])
                    else
                        for tt in keys(psd_expr[kk,nn].terms)
                            PSD[kk]["expr_val"][nn] += psd_expr[kk,nn].terms[tt]*JuMP.callback_value(cb_data,tt)
                        end
                    end
                end
            end
            PSDProjections(cp_model_info)
            added_cut=false
            if cp_model_info["orth_norm"] >= 1e-3
                added_cut=true
                for kk in keys(PSD)
                    con=JuMP.@build_constraint(
                    #(1/PSD[kk]["orth_norm"])*sum( PSD[kk]["ip"][nn]*PSD[kk]["orth_expr_val"][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) <= 0 
                        sum( PSD[kk]["ip"][nn]*PSD[kk]["sg"][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) >= 0 
                    )
                    MOI.submit(cp_model_info["model"], MOI.LazyConstraint(cb_data), con)
                end
            end
            if !haskey(feasXs,x_soln_str)
	            nXs += 1
                feasXs[x_soln_str]=Dict{String,Any}("value"=>-1e20,"bound_value"=>bd_val,"n_visits"=>1)
                feasXs[x_soln_str]["x_soln_str"] = x_soln_str
                println("Encountering solution: ", x_soln_str, " having bound value ",bd_val, " and infeas: ",cp_model_info["orth_norm"])
            else
                feasXs[x_soln_str]["n_visits"]+=1
                if mod(feasXs[x_soln_str]["n_visits"],100)==0
                    println("Revisiting solution: ",x_soln_str," having value: ", 
                        feasXs[x_soln_str]["value"], " and having bound value: ", bd_val, " and infeas: ",cp_model_info["orth_norm"])
                end
            end
            if !added_cut
                feasXs[x_soln_str]["value"] = bd_val
                println("Solution: ",x_soln_str," has value: ", feasXs[x_soln_str]["value"], " within tolerance, no cut added.")
	            if IncX["value"] < feasXs[x_soln_str]["value"]
	                    IncX = feasXs[x_soln_str]
	                    println("\t***New incumbent***",IncX["x_soln_str"]," with value ",round(IncX["value"];digits=5))
	            end
            end
        end
    end
    MOI.set(cp_model_info["model"], MOI.LazyConstraintCallback(), my_callback_function)
    MOI.set(cp_model_info["model"], MOI.NumberOfThreads(), 1)
    #MOI.set(cp_model_info["model"], CPLEX.CallbackFunction(), my_callback_function)
    JuMP.optimize!(cp_model_info["model"])
    cp_model_info["opt_val"]=JuMP.objective_value(cp_model_info["model"])
    cp_model_info["solve_status"]=JuMP.termination_status(cp_model_info["model"])

    end_time = time_ns()
    runtime = (end_time-start_time)/1e9
    println("Final best solution: ",IncX["x_soln_str"]," with value ",round(IncX["value"];digits=5))
    println("Runtime: ",runtime)
    #println("Optimal value: ",cp_model_info["opt_val"]," with status ", cp_model_info["solve_status"])
    #println("Number of nodes processed: ",nNodes)
    println("Number of feasibility problem solves: ",n_sdp_solves)
    #return bestUBVal,nNodes,runtime
end

testcase = Dict(
	"file" => "data/case9.m", 
 	"name" => "case9K3",  	
 	"attack_budget" => 3,
 	"inactive_indices" => [],
 	"protected_indices" => []
	)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

solveMaxminViaCB(pm_data,SparseSDPWRMPowerModel,use_dual_minmax=true)
#solveMaxminViaECPBnB(pm_data,SDPWRMPowerModel,use_dual_minmax=true)
