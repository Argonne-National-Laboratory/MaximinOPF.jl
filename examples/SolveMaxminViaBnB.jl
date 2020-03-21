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

function solveNodeSPFixed(model_info)
    model = model_info["model"]
    branch_ids = model_info["branch_ids"]
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
    PSD=model_info["psd_info"] 
    sg_info = Dict{Tuple{Int64,Int64},Array{Float64,1}}()
    for kk in keys(PSD)
        PSD[kk]["dual_val"][:] = JuMP.dual.(PSD[kk]["cref"])[:]
        sg_info[kk] = zeros(PSD[kk]["vec_len"])
        for nn=1:PSD[kk]["vec_len"]
            PSD[kk]["expr_val"][nn] = JuMP.value(psd_expr[kk,nn])
            sg_info[kk][nn] = PSD[kk]["dual_val"][nn]
        end
    end
    return sg_info
end #end of function

function add_cut_with_sg(model_info::Dict{String,Any}, sg_info)
    model = model_info["model"]
    psd_expr = model[:psd_expr]
    PSD=model_info["psd_info"] 
    for kk in keys(PSD)
        #if PSD[kk]["neg_eigs_sum"] < -1e-4 
        JuMP.@constraint(model, sum( PSD[kk]["ip"][nn]*sg_info[kk][nn]*psd_expr[kk,nn] for nn in 1:PSD[kk]["vec_len"]) >= 0 )
        #end
    end
end

function solveNodeSP_ECP(model_info, ndata; io=devnull)
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
    for kk in keys(PSD)
        for nn=1:PSD[kk]["vec_len"]
            PSD[kk]["expr_val"][nn] = JuMP.value(psd_expr[kk,nn])
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
    add_artificial_var_bds(psd_model_info; io=devnull)
	#add_psd_initial_cuts(psd_model_info;io=devnull)


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
	#add_psd_initial_cuts(cp_model_info;io=devnull)
    removePSD_Constraints(cp_model_info["psd_info"])
    add_artificial_var_bds(cp_model_info; io=devnull)

    BnBTree = Dict{Int64,Dict{String,Any}}()
    BnBTree[0] = init_node ### ADD ROOT NODE TO TREE
    nodekey=0

    feasXs = Dict{String,Dict{String,Any}}()
    feasXs[""]=Dict("x_soln"=>Dict{Int64,Float64}(),"x_soln_str"=>"", "value"=>0) ### CREATE INITIAL INCUMBENT SOLN
    IncX = feasXs[""]
    nXs = 1
    bestUBVal=1e20
    #cp_model_info=psd_model_info

    nNodes=1
    maxidx=1
    maxval=-1
    while true
        tree_report=string(" Nodes left: ",length(BnBTree)," out of ", nNodes,)
        incumbent_update=""
        currNode = pop!(BnBTree,nodekey)
        if currNode["bound_value"] <= IncX["value"]
            println("\nFathoming node $nodekey due to bound ",currNode["bound_value"]," <= ",IncX["value"])
        else
            printNode(currNode;pretext="\n",posttext=tree_report)
            min_sum=1e20
            inn_while_cntr = 0
            while true
                println("\tSolving node $nodekey subproblem...")
                model_info = cp_model_info
                if inn_while_cntr >= 20
                    model_info = psd_model_info
                    println("\tSwitching to direct psd node subproblem...")
                else
                    inn_while_cntr += 1
                end
                solveNodeSP_ECP(model_info, currNode)
                currNode["bound_value"] = model_info["opt_val"]
                if currNode["bound_value"] <= IncX["value"]
                    println("\tFathoming due to bound ",currNode["bound_value"]," <= ",IncX["value"])
                    break
                end
	            idx=findNextIndex(model_info["branch_ids"],model_info["x_soln"])
                if model_info["x_soln"][idx] == 0 || model_info["x_soln"][idx] == 1
                    x_soln_str = model_info["x_soln_str"]
                    is_new_soln = !haskey(feasXs,x_soln_str)
                    if is_new_soln
	                    nXs += 1
                        feasXs[x_soln_str]=Dict("bound_value"=>model_info["opt_val"])
                        feasXs[x_soln_str]["x_soln"] = copy(model_info["x_soln"])
                        feasXs[x_soln_str]["x_soln_str"] = x_soln_str
                    else
                        feasXs[x_soln_str]["bound_value"] = model_info["opt_val"]
                        if feasXs[x_soln_str]["bound_value"] - feasXs[x_soln_str]["value"] < 1e-4
                            break
                        end
                    end
                    sg_info=PSDSubgradient(model_info)  ### "Set 'total_sum_neg_eigs' key of model_info"
                    if model_info["total_sum_neg_eigs"] < -1e-4 
                        if abs(model_info["total_sum_neg_eigs"]) < min_sum
                            min_sum = abs(model_info["total_sum_neg_eigs"])
                        end
                        add_cut_with_sg(model_info, sg_info)
                        add_cut_with_sg(psd_model_info, sg_info)
	                    # APPLY A PRIMAL HEURISTIC
                        if is_new_soln
                            println("\tEvaluating new solution: ",x_soln_str)
                            sg_info = solveNodeSPFixed(psd_model_info)
	                        println("\t\t cp_optval=",model_info["opt_val"]," psd_optval=",psd_model_info["opt_val"])		
                            add_cut_with_sg(model_info, sg_info)
                            add_cut_with_sg(psd_model_info, sg_info)
                            PSDSubgradient(psd_model_info)  ### "Set 'total_sum_neg_eigs' key of model_info"
                            println("\t\tVerifying PSD feas: ",psd_model_info["total_sum_neg_eigs"])
                            feasXs[x_soln_str]["value"] = psd_model_info["heur_opt_val"]
	                        if IncX["value"] < feasXs[x_soln_str]["value"]
	                            IncX = feasXs[x_soln_str]
	                            incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"])
	                        end
                        end
                        println("\tPSD infeas: ", model_info["total_sum_neg_eigs"],", cp bound: ",
                        model_info["opt_val"],", versus known value: ",feasXs[x_soln_str]["value"])
                        continue
                    else
	                    println("\tFathoming due to optimality. Adding new attack solution: ",currNode["inactive_branches"])		
                        if is_new_soln
                            feasXs[x_soln_str]["value"] = model_info["opt_val"]
	                        if IncX["value"] < feasXs[x_soln_str]["value"]
	                            IncX = feasXs[x_soln_str]
	                            incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"])
	                        end
                        end
                        break
                    end
                else
	                # "BRANCH"
	                println("\tAdding two new nodes, the branching index is: ",idx)
                    node0=Dict{String,Any}("node_key"=>nNodes, "bound_value"=>currNode["bound_value"],"attacker_budget"=>currNode["attacker_budget"],
                        "inactive_branches"=>copy(currNode["inactive_branches"]), "protected_branches"=>copy(currNode["protected_branches"])
                    )
                    push!(node0["protected_branches"],idx)
	                BnBTree[nNodes]=node0
                    nNodes += 1

                    node1=Dict{String,Any}("node_key"=>nNodes, "bound_value"=>currNode["bound_value"],"attacker_budget"=>currNode["attacker_budget"],
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
