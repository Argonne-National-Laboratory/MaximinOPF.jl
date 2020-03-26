#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("../../MaximinOPF/src/MaximinOPF.jl")
include("../../MaximinOPF/src/utils.jl")
include("ProxSDP.jl")
using JuMP, MathOptInterface
using Mosek, MosekTools
using CPLEX
using JuMP
using PowerModels
PowerModels.silence()


function solveNodeMinmaxSP(model_info, ndata)
    model = model_info["model"]
    branch_ids = model_info["branch_ids"]
    psd_expr = model[:psd_expr]
    PSD=model_info["psd_info"] 

    unfix_vars(model,branch_ids)
    for l in ndata["inactive_branches"]
        fix(variable_by_name(model,"x[$l]_1"), 1; force=true)
    end
    for l in ndata["protected_branches"]
        fix(variable_by_name(model,"x[$l]_1"), 0; force=true)
    end
    JuMP.optimize!(model)
    model_info["opt_val"]=JuMP.objective_value(model)
    model_info["solve_status"]=JuMP.termination_status(model)
    model_info["x_soln_str"]=""
    for l in branch_ids
        x_var = variable_by_name(model,"x[$l]_1")
        x_val = JuMP.value(x_var)
        if x_val > 1.0-1.0e-8
	        x_val = 1
            model_info["x_soln_str"] = string(model_info["x_soln_str"]," $l")
        elseif x_val < 1.0e-8
	        x_val = 0
        end
        model_info["x_soln"][l] = x_val
    end
    ndata["x_soln"] = copy(model_info["x_soln"])
    ndata["bound_value"] = model_info["opt_val"]
    for kk in keys(PSD)
        PSD[kk]["dual_val"][:] = JuMP.dual.(PSD[kk]["cref"])[:]
        for nn=1:PSD[kk]["vec_len"]
            PSD[kk]["expr_val"][nn] = JuMP.value(psd_expr[kk,nn])
        end
    end
end #end of function

function solveNodeSPFixed(model_info; compute_psd_dual=true)
    model = model_info["model"]
    branch_ids = model_info["branch_ids"]
    PSD=model_info["psd_info"] 
    psd_expr = model[:psd_expr]
    unfix_vars(model,branch_ids)
    model_info["x_soln_str"]=""
    for l in branch_ids
        fix(variable_by_name(model,"x[$l]_1"), model_info["x_soln"][l]; force=true)
        if model_info["x_soln"][l] > 1.0-1.0e-8
            model_info["x_soln_str"] = string(model_info["x_soln_str"]," $l")
        end
    end
    JuMP.optimize!(model)
    model_info["opt_val"]=JuMP.objective_value(model)
    model_info["solve_status"]=JuMP.termination_status(model)
    model_info["heur_x_soln"] = copy(model_info["x_soln"])
    model_info["heur_opt_val"] = model_info["opt_val"]
    sg_info = Dict{Tuple{Int64,Int64},Array{Float64,1}}()
    if compute_psd_dual
        for kk in keys(PSD)
            PSD[kk]["dual_val"][:] = JuMP.dual.(PSD[kk]["cref"])[:]
            sg_info[kk] = zeros(PSD[kk]["vec_len"])
            for nn=1:PSD[kk]["vec_len"]
                PSD[kk]["expr_val"][nn] = JuMP.value(psd_expr[kk,nn])
                sg_info[kk][nn] = PSD[kk]["dual_val"][nn]
            end
        end
    end
    unfix_vars(model,branch_ids)
    return sg_info
end #end of function

function add_cuts(model_info::Dict{String,Any}, PSD::Dict{Tuple{Int64,Int64},Dict{String,Any}}; cut_type::String="sg", tol=1e-5)
    model = model_info["model"]
    psd_expr = model[:psd_expr]
    for kk in keys(PSD)
        if cut_type == "sg"  
            JuMP.@constraint(model, sum( PSD[kk]["ip"][nn]*PSD[kk][cut_type][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) >= 0 )
        elseif cut_type == "dual_val"
            vec_norm = norm(PSD[kk][cut_type][:])
            if vec_norm > tol
                JuMP.@constraint(model, (1/vec_norm)*sum( PSD[kk]["ip"][nn]*PSD[kk][cut_type][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) >= 0 )
            end
        elseif cut_type == "orth_expr_val"
            if PSD[kk]["orth_norm"] > tol
                JuMP.@constraint(model, 
                    (1/PSD[kk]["orth_norm"])*sum( PSD[kk]["ip"][nn]*PSD[kk][cut_type][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) <= 0 
                )
            end
        else
            println(Base.stderr,"add_cuts(): cut_type=",cut_type," is unrecognized.")
        end
    end
end

function solveNodeSP_ECP(model_info, ndata; compute_psd_dual=false, io=devnull)
    model = model_info["model"]
    branch_ids = model_info["branch_ids"]
    PSD = model_info["psd_info"]

    unfix_vars(model,branch_ids)
    for l in ndata["inactive_branches"]
        fix(variable_by_name(model,"x[$l]_1"), 1; force=true)
    end
    for l in ndata["protected_branches"]
        fix(variable_by_name(model,"x[$l]_1"), 0; force=true)
    end

    psd_expr = model[:psd_expr]
    JuMP.optimize!(model)
    sg_info = Dict{Tuple{Int64,Int64},Array{Float64,1}}()
    for kk in keys(PSD)
        for nn=1:PSD[kk]["vec_len"]
            PSD[kk]["expr_val"][nn] = JuMP.value(psd_expr[kk,nn])
        end
        if compute_psd_dual
            PSD[kk]["dual_val"][:] = JuMP.dual.(PSD[kk]["cref"])[:]
            sg_info[kk] = zeros(PSD[kk]["vec_len"])
            for nn=1:PSD[kk]["vec_len"]
                sg_info[kk][nn] = PSD[kk]["dual_val"][nn]
            end
        end
    end
    model_info["opt_val"]=JuMP.objective_value(model)
    model_info["solve_status"]=JuMP.termination_status(model)

    model_info["x_soln_str"]=""
    for l in branch_ids
        x_var = variable_by_name(model,"x[$l]_1")
        x_val = JuMP.value(x_var)
        if x_val > 1.0-1.0e-8
	        x_val = 1
            model_info["x_soln_str"] = string(model_info["x_soln_str"]," $l")
        elseif x_val < 1.0e-8
	        x_val = 0
        end
        model_info["x_soln"][l] = x_val
    end
    ndata["x_soln"] = copy(model_info["x_soln"])
    ndata["bound_value"] = model_info["opt_val"]
    return sg_info
end

function applyPrimalHeuristic(model_info)
    model = model_info["model"]
    branch_ids = model_info["branch_ids"]
    unfix_vars(model,branch_ids)
    idx=sortperm(sort(collect(keys(model_info["x_soln"]))), by=kk->model_info["x_soln"][kk]; rev=true)
    n_idx = length(idx)
    n_attack = min(model_info["attacker_budget"],n_idx)
    model_info["x_soln_str"]=""
    for kk in idx
        l=idx[kk]
        if kk <= n_attack && model_info["x_soln"][l] > 0
            model_info["heur_x_soln"][l] = 1
            fix(variable_by_name(model,"x[$l]_1"), 1; force=true)
            model_info["x_soln_str"] = string(model_info["x_soln_str"]," $l")
        else
            model_info["heur_x_soln"][l] = 0
            fix(variable_by_name(model,"x[$l]_1"), 0; force=true)
        end
    end
    JuMP.optimize!(model)
    model_info["heur_opt_val"]=JuMP.objective_value(model)
    model_info["solve_status"]=JuMP.termination_status(model)
    return model_info["heur_opt_val"]
end

function findNextNode(BnBTree::Dict{Int64,Dict{String,Any}})
  nodekey=-1
  weakestUBVal = -1e20
  for kk in keys(BnBTree) 
    if BnBTree[kk]["bound_value"] > weakestUBVal
        nodekey = kk
        weakestUBVal = BnBTree[kk]["bound_value"]
    end
  end
  return nodekey,weakestUBVal
end

function findNextIndex(undecided_branches,x_vals)
  return reduce((k1, k2) -> abs(x_vals[k1]-0.5) <= abs(x_vals[k2]-0.5) ? k1 : k2, undecided_branches)
end

function printNode(node::Dict{String,Any}; pretext="",posttext="")
    println(pretext,"Node ",node["node_key"]," has bound ",node["bound_value"],": inactive=",node["inactive_branches"]," protected: ",node["protected_branches"],posttext)
end

function solveMaxminViaBnB(pm_data,pm_form; use_dual_minmax=true)
    #global MAX_TIME
    K=pm_data["attacker_budget"]
    art_bd=10

    start_time = time_ns()
    init_node=Dict("inactive_branches"=>[], "protected_branches"=>[], "bound_value"=>1e20, "attacker_budget"=>K, "node_key"=>0) ### CREATE ROOT NODE

    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false, rm_rsoc=true, rm_therm_line_lim=false)
    branch_ids=sort(collect(pm_data["undecided_branches"]))
    println("Branch_ids: ",branch_ids)

    psd_base_maxmin = convertSOCtoPSD(base_maxmin)
    psd_optimizer=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=8)
    JuMP.set_optimizer(psd_base_maxmin,psd_optimizer)

    psd_model_info=Dict{String,Any}("branch_ids"=>branch_ids,"x_soln"=>Dict{Int64,Float64}(),
            "x_soln_str"=>"","heur_x_soln"=>Dict{Int64,Float64}(),"attacker_budget"=>K)
    for l in branch_ids
        psd_model_info["x_soln"][l]=0
    end
    psd_model_info["model"] = psd_base_maxmin
    gatherPSDConInfo(psd_model_info) ### "Sets the 'psd_info' key"
    #add_artificial_var_bds(psd_model_info; bd_mag=art_bd, io=devnull)

    prox_base_maxmin = convertSOCtoPSD(base_maxmin)
    prox_optimizer=with_optimizer(Ipopt.Optimizer)
    JuMP.set_optimizer(prox_base_maxmin,prox_optimizer)
    prox_model_info=Dict{String,Any}("branch_ids"=>branch_ids,"x_soln"=>psd_model_info["x_soln"],
            "x_soln_str"=>"","heur_x_soln"=>Dict{Int64,Float64}(),"attacker_budget"=>K)
    prox_model_info["model"] = prox_base_maxmin
    prepare_to_solve_PSD_via_ProxPt(prox_base_maxmin, prox_model_info; io=devnull)

    cp_base_maxmin = convertSOCtoPSD(base_maxmin)
    JuMP.set_optimizer(cp_base_maxmin,with_optimizer(CPLEX.Optimizer))
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_ScreenOutput",0)
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_Simplex_Tolerances_Optimality",1e-8)

    cp_model_info=Dict{String,Any}("branch_ids"=>branch_ids,"x_soln"=>psd_model_info["x_soln"],
            "x_soln_str"=>"","heur_x_soln"=>Dict{Int64,Float64}(),"attacker_budget"=>K)
    cp_model_info["model"] = cp_base_maxmin
    ### "Sets `model_info[\"cp_model\"] = cp_base_maxmin` and also the 'psd_info' entry"
    gatherPSDConInfo(cp_model_info) ### "Sets the 'psd_info' entry of model_info"
    #add_artificial_var_bds(cp_model_info; io=devnull)
    removePSD_Constraints(cp_model_info["psd_info"])
    #println(cp_model_info["model"])
    #print_var_bds(cp_model_info)
    #add_artificial_var_bds(cp_model_info; bd_mag=art_bd)

    BnBTree = Dict{Int64,Dict{String,Any}}()
    BnBTree[0] = init_node ### ADD ROOT NODE TO TREE
    nodekey=0

    feasXs = Dict{String,Dict{String,Any}}()
    feasXs[""]=Dict("x_soln"=>Dict{Int64,Float64}(),"x_soln_str"=>"", "bound_value"=>0, "value"=>0, "cheat_value"=>0) ### CREATE INITIAL INCUMBENT SOLN
    IncX = feasXs[""]
    nXs = 1
    solveNodeSP_ECP(psd_model_info, BnBTree[0]; compute_psd_dual=true)
    add_cuts(psd_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
    add_cuts(cp_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
    bound_obj(psd_model_info; bd_mag=4*psd_model_info["opt_val"])
    bound_obj(cp_model_info; bd_mag=4*psd_model_info["opt_val"])
	add_psd_initial_cuts(psd_model_info;io=devnull)
	add_psd_initial_cuts(cp_model_info;io=devnull)

    bestUBVal=1e20
    println("Initial processing of root node yields an upper bound of : ",bestUBVal)
    #cp_model_info=psd_model_info

    nNodes=1
    maxidx=1
    maxval=-1
    n_switch_direct = 0
    recycle_threshold = 10
    cheat_threshold = 0
    while true
        tree_report=string(" Nodes left: ",length(BnBTree)," out of ", nNodes,)
        incumbent_update=""
        currNode = pop!(BnBTree,nodekey)
        if currNode["bound_value"] <= IncX["value"]
            println("\nFathoming node $nodekey due to bound ",currNode["bound_value"]," <= ",IncX["value"])
        else
            printNode(currNode;pretext="\n",posttext=tree_report)
            inn_while_cntr = 0
            direct_mode=false
            while true
#=
                if inn_while_cntr >= recycle_threshold
                    println("\tRecycling current node")
	                BnBTree[nNodes]=currNode ### "Recycle current node"
                    nNodes += 1
                    break
=#
                model_info = cp_model_info
                if direct_mode
                    model_info = psd_model_info
                end
                solveNodeSP_ECP(model_info, currNode; compute_psd_dual=direct_mode)
#=
                test_passes=test_artificial_bds(model_info)
                if !test_passes
                    continue
                end
=#
                #PSDProjections(model_info)
                println("\tECP opt val: ",model_info["opt_val"], " with status: ", model_info["solve_status"])
                #println("\tECP opt val: ",model_info["opt_val"]," infeas: ", model_info["orth_norm"], " with status: ", model_info["solve_status"])
                            

                currNode["bound_value"] = min(model_info["opt_val"],currNode["bound_value"])
                if currNode["bound_value"] <= IncX["value"]
                    println("\tFathoming due to bound ",currNode["bound_value"]," <= ",IncX["value"])
                    break
                end
	            idx=findNextIndex(model_info["branch_ids"],model_info["x_soln"])
                if model_info["x_soln"][idx] == 0 || model_info["x_soln"][idx] == 1
                    x_soln_str = model_info["x_soln_str"]
                    if !haskey(feasXs,x_soln_str)
	                    nXs += 1
                        feasXs[x_soln_str]=Dict("bound_value"=>model_info["opt_val"], "value"=>0, "cheat_value"=>0)
                        feasXs[x_soln_str]["x_soln"] = copy(model_info["x_soln"])
                        feasXs[x_soln_str]["x_soln_str"] = x_soln_str
                        solveNodeSPFixed(psd_model_info)
                        #test_artificial_bds(psd_model_info)
                        feasXs[x_soln_str]["cheat_value"] = psd_model_info["heur_opt_val"]
                        println("New solution a priori known optimal value: ",feasXs[x_soln_str]["cheat_value"])
                    end
                    if !direct_mode
                        PSDProjections(cp_model_info)
                        # PSDSubgradient(cp_model_info)  ### "Set 'total_sum_neg_eigs' key of cp_model_info"
                        add_cuts(cp_model_info, cp_model_info["psd_info"]; cut_type="orth_expr_val")
                        add_cuts(psd_model_info, cp_model_info["psd_info"]; cut_type="orth_expr_val")
                        if cp_model_info["total_sum_neg_eigs"] < -1e-4 
                            println("\tSoln ",x_soln_str,": PSD infeas: ", cp_model_info["orth_norm"],", cp bound: ",
                                cp_model_info["opt_val"],", versus known value: ",feasXs[x_soln_str]["cheat_value"])
                            if inn_while_cntr >= cheat_threshold  
                                direct_mode=true
                                println("\tSwitching to direct psd node subproblem...")
                                n_switch_direct += 1
                            else    
                                inn_while_cntr += 1
                            end
                            continue
                        else
	                        println("\tFathoming due to optimality. Adding new attack solution: ",currNode["inactive_branches"])		
                            feasXs[x_soln_str]["value"] = cp_model_info["opt_val"]
	                        if IncX["value"] < feasXs[x_soln_str]["value"]
	                            IncX = feasXs[x_soln_str]
	                            incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"])
	                        end
                            break
                        end
                    else
	                    println("\tFathoming due to optimality. Adding new attack solution: ",currNode["inactive_branches"])		
                        add_cuts(cp_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
                        add_cuts(psd_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
                        feasXs[x_soln_str]["value"] = psd_model_info["opt_val"]
	                    if IncX["value"] < feasXs[x_soln_str]["value"]
	                        IncX = feasXs[x_soln_str]
	                        incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"])
	                    end
                        break
                    end
                else
	                # "BRANCH"
                    bd_val=currNode["bound_value"]
	                println("\tAdding two new nodes, the branching index is: ",idx)
                    node0=Dict{String,Any}("node_key"=>nNodes, "bound_value"=>bd_val,"attacker_budget"=>currNode["attacker_budget"],
                        "inactive_branches"=>copy(currNode["inactive_branches"]), "protected_branches"=>copy(currNode["protected_branches"])
                    )
                    push!(node0["protected_branches"],idx)
	                BnBTree[nNodes]=node0
                    nNodes += 1

                    node1=Dict{String,Any}("node_key"=>nNodes, "bound_value"=>bd_val,"attacker_budget"=>currNode["attacker_budget"],
                        "inactive_branches"=>copy(currNode["inactive_branches"]), "protected_branches"=>copy(currNode["protected_branches"])
                    )
                    push!(node1["inactive_branches"],idx)
	                node1["attacker_budget"] -= 1 ### This node will have another branch inactive, debiting from the inherited budget
	                BnBTree[nNodes]=node1
                    nNodes += 1
                    break
                end
            end ### "while loop"
            println("\tNew node bound value: ",currNode["bound_value"])
        end ### "if initial fathoming due to bound"

        if length(BnBTree) > 0
            nodekey,bestUBVal=findNextNode(BnBTree)
        else
            bestUBVal = IncX["value"]
            break
        end
#=
        if (time_ns()-start_time)/1e9 > MAX_TIME
            break
        end
=#
        println("\tBound gap: [",IncX["value"],", ",bestUBVal,"]",incumbent_update)
        if bestUBVal < IncX["value"]
            println("Terminating since all other nodes will be fathomed by bound")
            break
        end
    end # while
    end_time = time_ns()

    runtime = (end_time-start_time)/1e9
    println("Final best solution: ",IncX["x_soln_str"]," with value ",IncX["value"])
    println("Runtime: ",runtime)
    println("Number of nodes processed: ",nNodes)
    println("Number of times cheating: ",n_switch_direct)
    return bestUBVal,nNodes,runtime
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

solveMaxminViaBnB(pm_data,SparseSDPWRMPowerModel,use_dual_minmax=true)
#solveMaxminViaBnB(pm_data,SDPWRMPowerModel,use_dual_minmax=true)
