#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("../../MaximinOPF/src/MaximinOPF.jl")
include("../../MaximinOPF/src/utils.jl")
include("../../MaximinOPF/src/utils.jl")
include("psd_utils.jl")
using JuMP, MathOptInterface
using Mosek, MosekTools
using CPLEX
using JuMP
using PowerModels
PowerModels.silence()


function solveNodeSP(model_info, ndata; fix_all_x=false, compute_psd_dual=false, io=devnull)
    model = model_info["model"]
    branch_ids = model_info["branch_ids"]

    unfix_vars(model,branch_ids)
    for l in branch_ids
        if l in ndata["inactive_branches"]
            fix(variable_by_name(model,"x[$l]_1"), 1; force=true)
	elseif l in ndata["protected_branches"] || fix_all_x
            fix(variable_by_name(model,"x[$l]_1"), 0; force=true)
	end
    end

    JuMP.optimize!(model)
    sg_info = Dict{Tuple{Int64,Int64},Array{Float64,1}}()
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

function deleteNodesByBound(BnBTree,bdval)
    for bb in keys(BnBTree)
        if BnBTree[bb]["bound_value"] <= bdval
            delete!(BnBTree,bb)
        end
    end
end

function solveMaxminViaBnB(pm_data,pm_form; use_dual_minmax=true)
    #global MAX_TIME
    K=pm_data["attacker_budget"]
    art_bd=10

    start_time = time_ns()
    init_node=Dict("inactive_branches"=>[], "protected_branches"=>[], "bound_value"=>1e20, "attacker_budget"=>K, "node_key"=>0) ### CREATE ROOT NODE

    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false, rm_rsoc=true, rm_therm_line_lim=true)
    branch_ids=sort(collect(pm_data["undecided_branches"]))
    println("Branch_ids: ",branch_ids)
    if pm_data["use_rf"]
	println("Using reversible flow relaxation.")
	for l in branch_ids
    	    ii,jj = pm_data["branch"]["$l"]["f_bus"],pm_data["branch"]["$l"]["t_bus"]

            pi_f=JuMP.variable_by_name(base_maxmin,"pi_f[$l]_1")
            alpha_ii=JuMP.variable_by_name(base_maxmin,"p_bal[$ii]_1")
            pi_t=JuMP.variable_by_name(base_maxmin,"pi_t[$l]_1")
            alpha_jj=JuMP.variable_by_name(base_maxmin,"p_bal[$jj]_1")
	    JuMP.@constraint(base_maxmin, pi_f + alpha_ii - pi_t - alpha_jj == 0)

            phi_f=JuMP.variable_by_name(base_maxmin,"phi_f[$l]_1")
            beta_ii=JuMP.variable_by_name(base_maxmin,"q_bal[$ii]_1")
            phi_t=JuMP.variable_by_name(base_maxmin,"phi_t[$l]_1")
            beta_jj=JuMP.variable_by_name(base_maxmin,"q_bal[$jj]_1")
	    JuMP.@constraint(base_maxmin, phi_f + beta_ii - phi_t - beta_jj == 0)
	end
    end

    psd_optimizer=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=1)
    JuMP.set_optimizer(base_maxmin,psd_optimizer)

    psd_model_info=Dict{String,Any}("branch_ids"=>branch_ids,"x_soln"=>Dict{Int64,Float64}(),
            "x_soln_str"=>"","heur_x_soln"=>Dict{Int64,Float64}(),"attacker_budget"=>K)
    for l in branch_ids
        psd_model_info["x_soln"][l]=0
    end
    psd_model_info["model"] = base_maxmin

    BnBTree = Dict{Int64,Dict{String,Any}}()
    BnBTree[0] = init_node ### ADD ROOT NODE TO TREE
    nodekey=0

    feasXs = Dict{String,Dict{String,Any}}()
    feasXs[""]=Dict("x_soln"=>Dict{Int64,Float64}(),"x_soln_str"=>"", "bound_value"=>0, "value"=>0) ### CREATE INITIAL INCUMBENT SOLN
    IncX = feasXs[""]
    nXs = 1

    bestUBVal=1e20
    println("Initial processing of root node yields an upper bound of : ",bestUBVal)

    nNodes=1
    maxidx=1
    maxval=-1

    while true
        incumbent_update=""
        currNode = pop!(BnBTree,nodekey)
        if currNode["bound_value"] <= IncX["value"]
            println("\nFathoming node $nodekey due to bound ",currNode["bound_value"]," <= ",IncX["value"])
        else
            printNode(currNode;pretext="\n",posttext="")
            inn_while_cntr = 0
            while true
	x_soln_str=""
	for l in branch_ids
            psd_model_info["x_soln"][l]=0
	end
    	cut_lines=sort(collect(currNode["inactive_branches"]))
    	for l in cut_lines
            x_soln_str=string(x_soln_str," $l")
            psd_model_info["x_soln"][l]=1
    	end
        if !haskey(feasXs,x_soln_str)
	    nXs += 1
            feasXs[x_soln_str]=Dict{String,Any}("x_soln_str"=>x_soln_str)
            solveNodeSP(psd_model_info, currNode; fix_all_x=true, compute_psd_dual=false)
            feasXs[x_soln_str]["value"] = psd_model_info["opt_val"]
            println("\tNew solution ",x_soln_str," has known optimal value: ",round(feasXs[x_soln_str]["value"];digits=5))
	    if IncX["value"] < feasXs[x_soln_str]["value"]
	    	IncX = feasXs[x_soln_str]
	        incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"]," with value ",round(IncX["value"];digits=5))
                deleteNodesByBound(BnBTree,IncX["value"])
	    end
	end
                solveNodeSP(psd_model_info, currNode; compute_psd_dual=false)

                currNode["bound_value"] = min(psd_model_info["opt_val"],currNode["bound_value"])
                if currNode["bound_value"] <= IncX["value"]
                    println("\tFathoming due to bound ",currNode["bound_value"]," <= ",IncX["value"])
                    break
                end
	            idx=findNextIndex(psd_model_info["branch_ids"],psd_model_info["x_soln"])
                if psd_model_info["x_soln"][idx] == 0 || psd_model_info["x_soln"][idx] == 1
                    x_soln_str = psd_model_info["x_soln_str"]
                    if !haskey(feasXs,x_soln_str)
	                    nXs += 1
                        feasXs[x_soln_str]=Dict("bound_value"=>psd_model_info["opt_val"], "value"=>psd_model_info["opt_val"])
                        feasXs[x_soln_str]["x_soln"] = copy(psd_model_info["x_soln"])
                        feasXs[x_soln_str]["x_soln_str"] = x_soln_str
                        println("New solution has known optimal value: ",feasXs[x_soln_str]["value"])
                    end
	                println("\tFathoming due to optimality. Adding new attack solution: ",currNode["inactive_branches"])		
                    feasXs[x_soln_str]["value"] = psd_model_info["opt_val"]
	                if IncX["value"] < feasXs[x_soln_str]["value"]
	                    IncX = feasXs[x_soln_str]
	                    incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"])
                        deleteNodesByBound(BnBTree,IncX["value"])
	                end
                    break
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
        tree_report=string(" Nodes left: ",length(BnBTree)," out of ", nNodes,)
        println("\tBound gap: [",IncX["value"],", ",bestUBVal,"]",incumbent_update,",\t",tree_report)
        if bestUBVal < IncX["value"]
            println("Terminating since all other nodes will be fathomed by bound")
            break
        end
    end # while

    end_time = time_ns()
    runtime = (end_time-start_time)/1e9

    println("Final best solution: ",IncX["x_soln_str"]," with value ",IncX["value"])
    println("Final bound gap: [",IncX["value"],", ",bestUBVal,"], percent gap: ",round(100*(bestUBVal-IncX["value"])/IncX["value"];digits=1))
    println("Runtime: ",runtime)
    println("Number of nodes processed: ",nNodes)
    println("case \t budget \t nNodes \t sec_per_node \t secs ")
    println(" & $K & $nNodes & ",round(runtime/nNodes;digits=2),
	" & ",round(runtime),"\\\\ \n")
    return bestUBVal,nNodes,runtime
end

global case_instance="30"
global attack_budget=4
global form_str="SDP"
global use_rf_relax=false
global case_spec=false
global budget_spec=false
for aa in 1:length(ARGS)
    global case_instance
    global attack_budget
    global use_rf_relax
    if occursin("--case=",ARGS[aa])
            case_instance=ARGS[aa][(length("--case=")+1):length(ARGS[aa])]
            println("case being set to ",case_instance)
            case_spec=true
    elseif occursin("--K=",ARGS[aa])
            attack_budget=parse(Int64,ARGS[aa][(length("--K=")+1):length(ARGS[aa])])
            println("attack budget being set to ",attack_budget)
            budget_spec=true
    elseif occursin("--use_rf",ARGS[aa])
    	if occursin("=no",ARGS[aa]) || occursin("=false",ARGS[aa])
            use_rf_relax=false
	else
            use_rf_relax=true
	end
    elseif occursin("--form_str=",ARGS[aa])
            form_str=ARGS[aa][(length("--form_str=")+1):length(ARGS[aa])]
    else
            println(Base.stderr,"Argument ",ARGS[aa]," not recognized.")
    end
end

testcase = Dict(
    "file" => string("data/case",case_instance,".m"), 
    "PMOption" => SparseSDPWRMPowerModel,
    "name" => string("case",case_instance,"K",attack_budget,form_str),  	
    "attack_budget" => attack_budget,
    "inactive_indices" => [],
    "protected_indices" => [],
    "use_rf"=>use_rf_relax
)
println("Testing ADMM for root node relaxation of Maxmin problem: ")
println(testcase)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry
pm_data["use_rf"]=testcase["use_rf"]

solveMaxminViaBnB(pm_data,SparseSDPWRMPowerModel,use_dual_minmax=true)
#solveMaxminViaBnB(pm_data,SDPWRMPowerModel,use_dual_minmax=true)
