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
using JuMP
using PowerModels
PowerModels.silence()

function solveNodeSP_ECP(model_info, ndata; compute_projection=false, compute_psd_dual=false, io=devnull)
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
        if compute_psd_dual
            PSD[kk]["dual_val"][:] = JuMP.dual.(PSD[kk]["cref"])[:]
        end
    end
    if compute_projection
        PSDProjections(model_info)
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

function solveNodeSPFixed(model_info; compute_projection=true, compute_psd_dual=false, io=devnull)
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
    for kk in keys(PSD)
        for nn=1:PSD[kk]["vec_len"]
            PSD[kk]["expr_val"][nn] = JuMP.value(psd_expr[kk,nn])
        end
        if compute_psd_dual
            PSD[kk]["dual_val"][:] = JuMP.dual.(PSD[kk]["cref"])[:]
        end
    end
    if compute_projection
        PSDProjections(model_info)
    end
    unfix_vars(model,branch_ids)
end #end of function

### "Precondition: Make sure that the relevant data structure, be it 'sg', 'dual_val', or 'orth_expr_val' has been computed"
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
    model_info["opt_val"]=JuMP.objective_value(model)
    model_info["solve_status"]=JuMP.termination_status(model)
    return model_info["opt_val"]
end

function deleteNodesByBound(BnBTree,bdval)
    for bb in keys(BnBTree)
        if BnBTree[bb]["bound_value"] <= bdval
            delete!(BnBTree,bb)
        end
    end
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
    println(pretext,"Node ",node["node_key"]," has bound ",round(node["bound_value"];digits=5),
        ": inactive=",node["inactive_branches"]," protected: ",node["protected_branches"],posttext)
end

function getMaxShedding(pm_data)
    maxShed=0.0
    for ll in keys(pm_data["load"])
        maxShed += abs(pm_data["load"][ll]["pd"]) + abs(pm_data["load"][ll]["qd"])
    end
    for gg in keys(pm_data["gen"])
        maxShed += max(pm_data["gen"][gg]["pmin"],0) + max(pm_data["gen"][gg]["qmin"],0)
    end
    return maxShed 
end

function solveMaxminViaECPBnB(pm_data,pm_form; use_dual_minmax=true)
    #global MAX_TIME
    K=pm_data["attacker_budget"]
    art_bd=ceil(getMaxShedding(pm_data);digits=6)

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


    cp_base_maxmin = convertSOCtoPSD(base_maxmin)
    JuMP.set_optimizer(cp_base_maxmin,with_optimizer(CPLEX.Optimizer))
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_ScreenOutput",0)
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_Simplex_Tolerances_Optimality",1e-8)

    cp_model_info=Dict{String,Any}("branch_ids"=>branch_ids,"x_soln"=>psd_model_info["x_soln"],
            "x_soln_str"=>"","heur_x_soln"=>Dict{Int64,Float64}(),"attacker_budget"=>K)
    cp_model_info["model"] = cp_base_maxmin
    ### "Sets `model_info[\"cp_model\"] = cp_base_maxmin` and also the 'psd_info' entry"
    gatherPSDConInfo(cp_model_info) ### "Sets the 'psd_info' entry of model_info"
    removePSD_Constraints(cp_model_info["psd_info"])
    #println(cp_model_info["model"])

    BnBTree = Dict{Int64,Dict{String,Any}}()
    BnBTree[0] = init_node ### ADD ROOT NODE TO TREE
    nodekey=0

    feasXs = Dict{String,Dict{String,Any}}()
    feasXs[""]=Dict("x_soln"=>Dict{Int64,Float64}(),"x_soln_str"=>"", "bound_value"=>0, "value"=>0) ### CREATE INITIAL INCUMBENT SOLN
    IncX = feasXs[""]
    nXs = 1
    solveNodeSPFixed(psd_model_info; compute_projection=false, compute_psd_dual=true)
    n_sdp_solves = 1
    add_cuts(cp_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
    feasXs[""]["value"] = psd_model_info["opt_val"]

    
    bound_obj(cp_model_info; bd_mag=art_bd)
    println("Artifical bound on the objective function: ",art_bd)
	add_psd_initial_cuts(cp_model_info;io=devnull)

    bestUBVal=1e20
    println("Initial processing of root node yields an upper bound of : ",round(bestUBVal;digits=5))

    nNodes=1
    maxidx=1
    maxval=-1
    n_extra_cuts = 0
    while length(BnBTree) > 0
        nodekey,bestUBVal=findNextNode(BnBTree)
        incumbent_update=""
        currNode = pop!(BnBTree,nodekey)
        if currNode["bound_value"] <= IncX["value"]
            println("\nFathoming node $nodekey due to bound ",round(currNode["bound_value"];digits=5)," <= ",round(IncX["value"];digits=5))
            tree_report=string(" Nodes left: ",length(BnBTree)," out of ", nNodes,)
            println("\tBound gap: [",round(IncX["value"];digits=5),", ",round(bestUBVal;round=5),"]",incumbent_update,",\t",tree_report)
            continue
        end
        printNode(currNode;pretext="\n",posttext="")
        while true
            solveNodeSP_ECP(cp_model_info, currNode; compute_projection=false, compute_psd_dual=false)
            println("\tNew node bound: ",round(cp_model_info["opt_val"];digits=5), 
                " with status: ", cp_model_info["solve_status"])
            for nn=1:n_extra_cuts
                PSDProjections(cp_model_info)
                add_cuts(cp_model_info, cp_model_info["psd_info"];cut_type="orth_expr_val")
                #add_cuts(psd_model_info, cp_model_info["psd_info"];cut_type="orth_expr_val")
                solveNodeSP_ECP(cp_model_info, currNode; compute_psd_dual=false)
                println("\tNew node bound: ",round(cp_model_info["opt_val"];digits=5), " with status: ", cp_model_info["solve_status"])
            end

            currNode["bound_value"] = min(cp_model_info["opt_val"],currNode["bound_value"])
            if currNode["bound_value"] <= IncX["value"]
                println("\tFathoming due to bound ",round(currNode["bound_value"];digits=5)," <= ",round(IncX["value"];digits=5))
                break
            end
            idx=findNextIndex(cp_model_info["branch_ids"],cp_model_info["x_soln"])
            if cp_model_info["x_soln"][idx] == 0 || cp_model_info["x_soln"][idx] == 1
                x_soln_str = cp_model_info["x_soln_str"]
                if !haskey(feasXs,x_soln_str)
	                nXs += 1
                    feasXs[x_soln_str]=Dict("bound_value"=>cp_model_info["opt_val"], "value"=>0)
                    feasXs[x_soln_str]["x_soln"] = copy(cp_model_info["x_soln"])
                    feasXs[x_soln_str]["x_soln_str"] = x_soln_str
                    solveNodeSPFixed(psd_model_info; compute_projection=false, compute_psd_dual=true)
                    n_sdp_solves += 1
                    add_cuts(cp_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
                    #add_cuts(psd_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
                    feasXs[x_soln_str]["value"] = psd_model_info["opt_val"]
                    println("\tNew solution has known optimal value: ",round(feasXs[x_soln_str]["value"];digits=5))
                    continue
                else
                    println("\tSoln: ",x_soln_str, ": cp_val ",cp_model_info["opt_val"], " vs psd_val: ",feasXs[x_soln_str]["value"])
	                println("\tFathoming due to optimality. Adding new attack solution: ",currNode["inactive_branches"])		
	                if IncX["value"] < feasXs[x_soln_str]["value"]
	                    IncX = feasXs[x_soln_str]
	                    incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"]," with value ",round(IncX["value"];digits=5))
                        deleteNodesByBound(BnBTree,IncX["value"])
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

#=
        if (time_ns()-start_time)/1e9 > MAX_TIME
            break
        end
=#
        tree_report=string(" Nodes left: ",length(BnBTree)," out of ", nNodes,)
        println("\tBound gap: [",round(IncX["value"];digits=5),", ",round(bestUBVal;digits=5),
            "]",incumbent_update,",\t",tree_report)
    end # "while BnB tree not empty"
    bestUBVal = IncX["value"]
    end_time = time_ns()

    runtime = (end_time-start_time)/1e9
    println("Final best solution: ",IncX["x_soln_str"]," with value ",round(IncX["value"];digits=5))
    println("Runtime: ",runtime)
    println("Number of nodes processed: ",nNodes)
    println("Number of feasibility problem solves: ",n_sdp_solves)
    return bestUBVal,nNodes,runtime
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

solveMaxminViaECPBnB(pm_data,SparseSDPWRMPowerModel,use_dual_minmax=true)
#solveMaxminViaECPBnB(pm_data,SDPWRMPowerModel,use_dual_minmax=true)
