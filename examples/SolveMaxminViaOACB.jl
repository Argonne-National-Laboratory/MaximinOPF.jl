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

function solveMaxminViaOACB(pm_data,pm_form)
    #global MAX_TIME
    K=pm_data["attacker_budget"]

    start_time = time_ns()

    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false, rm_rsoc=true, rm_therm_line_lim=true)
    global branch_ids=sort(collect(pm_data["undecided_branches"]))
    println("Branch_ids: ",branch_ids)

    psd_base_maxmin = convertSOCtoPSD(base_maxmin)
    psd_optimizer=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=8)
    JuMP.set_optimizer(psd_base_maxmin,psd_optimizer)

    global psd_model_info=Dict{String,Any}("branch_ids"=>branch_ids,"x_soln"=>Dict{Int64,Float64}(),
            "x_soln_str"=>"","heur_x_soln"=>Dict{Int64,Float64}(),"attacker_budget"=>K)
    for l in branch_ids
        psd_model_info["x_soln"][l]=0
    end
    psd_model_info["model"] = psd_base_maxmin
    gatherPSDConInfo(psd_model_info) ### "Sets the 'psd_info' key"


    cp_base_maxmin = convertSOCtoPSD(base_maxmin)
    JuMP.set_optimizer(cp_base_maxmin,with_optimizer(CPLEX.Optimizer))
    #JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_ScreenOutput",0)
    #JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_Simplex_Tolerances_Optimality",1e-8)
	JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_Parallel", 1)
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_Threads", 1)
    global MIP_INT_TOL=1e-5
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_MIP_Tolerances_Integrality", MIP_INT_TOL)
    #JuMP.set_optimizer(cp_base_maxmin,with_optimizer(GLPK.Optimizer))

    global cp_model_info=Dict{String,Any}("branch_ids"=>branch_ids,"x_soln"=>psd_model_info["x_soln"],
            "x_soln_str"=>"","heur_x_soln"=>Dict{Int64,Float64}(),"attacker_budget"=>K)
    cp_model_info["model"] = cp_base_maxmin
    ### "Sets `model_info[\"cp_model\"] = cp_base_maxmin` and also the 'psd_info' entry"
    gatherPSDConInfo(cp_model_info) ### "Sets the 'psd_info' entry of model_info"
    removePSD_Constraints(cp_model_info["psd_info"])
    
    global feasXs = Dict{String,Dict{String,Any}}()
    feasXs[""]=Dict("x_soln"=>Dict{Int64,Float64}(),"x_soln_str"=>"", "bound_value"=>0, "value"=>0) ### CREATE INITIAL INCUMBENT SOLN

    for l in branch_ids
        psd_model_info["x_soln"][l] = 0
    end
    global n_sdp_solves = 1
    global IncX = feasXs[""]
    global nXs = 1
    feasXs[""]["value"] = 0.0

    solveSP(psd_model_info; fix_x=true, compute_projection=false, compute_psd_dual=true)
    PSD=psd_model_info["psd_info"]
    psd_expr = cp_model_info["model"][:psd_expr]
    for kk in keys(PSD)
        vec_norm = norm(PSD[kk]["dual_val"][:])
        if vec_norm > 1e-10
            JuMP.@constraint(cp_model_info["model"], (1/vec_norm)*sum( PSD[kk]["ip"][nn]*PSD[kk]["dual_val"][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) >= 0 )
        end
    end

    solveSP(psd_model_info; fix_x=false,compute_projection=false, compute_psd_dual=true)
    PSD=psd_model_info["psd_info"]
    psd_expr = cp_model_info["model"][:psd_expr]
    for kk in keys(PSD)
        vec_norm = norm(PSD[kk]["dual_val"][:])
        if vec_norm > 1e-10
            JuMP.@constraint(cp_model_info["model"], (1/vec_norm)*sum( PSD[kk]["ip"][nn]*PSD[kk]["dual_val"][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) >= 0 )
        end
    end
    n_sdp_solves += 1
    
	add_psd_initial_cuts(cp_model_info; bdmag=1e3,io=devnull)
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
        global psd_model_info
        psd_expr = cp_model_info["model"][:psd_expr]
        x_is_int=true
        tol = 1e-5
        x_soln_str=""
        #println("callback at x_soln:")
        bd_val = 0.0
        obj_expr=JuMP.objective_function(cp_model_info["model"], AffExpr)
        for tt in keys(obj_expr.terms)
            bd_val += obj_expr.terms[tt]*JuMP.callback_value(cb_data, tt)
        end
        for l in branch_ids
            psd_model_info["x_soln"][l] = JuMP.callback_value(cb_data, x_vars[l])
            #print(" ",psd_model_info["x_soln"][l])
            if psd_model_info["x_soln"][l] <= MIP_INT_TOL
                psd_model_info["x_soln"][l] = 0
            elseif psd_model_info["x_soln"][l] >= 1-MIP_INT_TOL
                psd_model_info["x_soln"][l] = 1
                x_soln_str = string(x_soln_str," $l")
            else
                x_is_int=false
            end
        end
        #println("\n")
        if x_is_int
            if !haskey(feasXs,x_soln_str)
	            nXs += 1
                solveSP(psd_model_info; fix_x=true, compute_projection=false, compute_psd_dual=true)
                n_sdp_solves += 1
                feasXs[x_soln_str]=Dict("value"=>psd_model_info["opt_val"])
                feasXs[x_soln_str]["x_soln_str"] = x_soln_str
                incumbent_update=""
	            if IncX["value"] < feasXs[x_soln_str]["value"]
	                    IncX = feasXs[x_soln_str]
	                    incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"]," with value ",round(IncX["value"];digits=5))
	            end
                println("Encountering solution: ", x_soln_str, " having value: ",round(psd_model_info["opt_val"];digits=4),incumbent_update)
                PSD=psd_model_info["psd_info"]
                for kk in keys(PSD)
                    vec_norm = norm(PSD[kk]["dual_val"][:])
                    if vec_norm > 1e-10
                        con=JuMP.@build_constraint((1/vec_norm)*sum( PSD[kk]["ip"][nn]*PSD[kk]["dual_val"][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) >= 0 )
                        MOI.submit(cp_model_info["model"], MOI.LazyConstraint(cb_data), con)
                    end
                end
            elseif bd_val - feasXs[x_soln_str]["value"] > 1e-2
                println("Revisiting solution: ",x_soln_str," having value: ", 
                    round(feasXs[x_soln_str]["value"];digits=4), 
                    " and having bound value: ", round(bd_val;digits=4))
                solveSP(psd_model_info; fix_x=true, compute_projection=false, compute_psd_dual=true)
                n_sdp_solves += 1
                PSD=psd_model_info["psd_info"]
                for kk in keys(PSD)
                    vec_norm = norm(PSD[kk]["dual_val"][:])
                    if vec_norm > 1e-10
                        con=JuMP.@build_constraint((1/vec_norm)*sum( PSD[kk]["ip"][nn]*PSD[kk]["dual_val"][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) >= 0 )
                        MOI.submit(cp_model_info["model"], MOI.LazyConstraint(cb_data), con)
                    end
                end
            end
        end
    end
    MOI.set(cp_model_info["model"], MOI.LazyConstraintCallback(), my_callback_function)
    MOI.set(cp_model_info["model"], MOI.NumberOfThreads(), 1)
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
	"file" => "data/case57.m", 
 	"name" => "case57K2",  	
 	"attack_budget" => 2,
 	"inactive_indices" => [],
 	"protected_indices" => []
	)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

#solveMaxminViaOACB(pm_data,SDPWRMPowerModel)
solveMaxminViaOACB(pm_data,SparseSDPWRMPowerModel)
