#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("../../MaximinOPF/src/MaximinOPF.jl")
include("../../MaximinOPF/src/utils.jl")
using Mosek, MosekTools
using JuMP
using PowerModels


function solveNodeMinmaxSP(pm_data,pm_form,ndata; use_dual=true)
  pm_data["inactive_branches"] = ndata["inactive_branches"] ###Adding another key and entry
  pm_data["protected_branches"] = ndata["protected_branches"] ###Adding another key and entry
  pm_data["attacker_budget"]=ndata["attacker_budget"]
  pm_data["x_vals"]=Dict{Int64,Float64}()
  
  model=nothing
  if pm_data["attacker_budget"] > 0
    if use_dual
      model,pm = MaximinOPF.SolveMinmaxDual(pm_data, pm_form,with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0))
    else
      model,pm = MaximinOPF.SolveMinmax(pm_data, pm_form,with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0))
    end
  else
    pm = MaximinOPF.SolveFP(pm_data,pm_form,with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0))
    pm_data["undecided_branches"]= []
    model=pm.model
  end

  ndata["bound_value"]=JuMP.objective_value(model)
  return pm
end #end of function

function solveNodeMinmaxSP(model_info, ndata)
    unfix_vars(model_info["model"],model_info["branch_ids"])
    for l in ndata["inactive_branches"]
        fix(variable_by_name(model_info["model"],"x[$l]_1"), 1; force=true)
    end
    for l in ndata["protected_branches"]
        fix(variable_by_name(model_info["model"],"x[$l]_1"), 0; force=true)
    end
    JuMP.optimize!(model_info["model"])
    model_info["opt_val"]=JuMP.objective_value(model_info["model"])
    model_info["solve_status"]=JuMP.termination_status(model_info["model"])
    for l in model_info["branch_ids"]
        x_var = variable_by_name(model_info["model"],"x[$l]_1")
        x_val = JuMP.value(x_var)
        if x_val > 1.0-1.0e-8
	        x_val = 1
            #model_info["x_soln_01"][l] = 1
        elseif x_val < 1.0e-8
	        x_val = 0
        end
        model_info["x_soln"][l] = x_val
    end
    ndata["x_soln"] = copy(model_info["x_soln"])
    ndata["bound_value"] = model_info["opt_val"]
end #end of function

function applyPrimalHeuristic(pm_data,pm_form)
    fp_pm = MaximinOPF.SolveFP(pm_data,pm_form,with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0))
    return JuMP.objective_value(fp_pm.model)
end

function findNextNode(E)
  nodekey=-1e20
  weakestUBVal = -1e20
  for (k,n) in E
    if n[2]["bound_value"] > weakestUBVal
    nodekey = n[1]
        weakestUBVal = n[2]["bound_value"]
    end
  end
  return nodekey,weakestUBVal
end

function findNextIndex(undecided_branches,x_vals)
  #return reduce((k1, k2) -> x_vals[k1] <= x_vals[k2] ? k1 : k2, undecided_branches)
  #return reduce((k1, k2) -> x_vals[k1] >= x_vals[k2] ? k1 : k2, undecided_branches)
  return reduce((k1, k2) -> abs(x_vals[k1]-0.5) <= abs(x_vals[k2]-0.5) ? k1 : k2, undecided_branches)
end

function solveMaxminViaBnB(pm_data,pm_form; use_dual_minmax=true)
    #global MAX_TIME
    K=pm_data["attacker_budget"]


    start_time = time_ns()
    init_node=Dict("inactive_branches"=>[], "protected_branches"=>[], "bound_value"=>1e20, "attacker_budget"=>K) ### CREATE ROOT NODE

    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false, rm_rsoc=true, rm_therm_line_lim=false)
    psd_base_maxmin = convertSOCtoPSD(base_maxmin)
    pm_optimizer=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=8)
    JuMP.set_optimizer(psd_base_maxmin,pm_optimizer)

    model_info=Dict{String,Any}()
    model_info["model"] = psd_base_maxmin
    model_info["branch_ids"] = pm_data["undecided_branches"]
    model_info["x_soln"] = Dict{Int64,Float64}()

    BnBTree = Dict()
    BnBTree[0] = init_node ### ADD ROOT NODE TO TREE
    nodekey=0

    feasXs = Dict()
    feasXs[1]=Dict("inactive_branches"=>[], "protected_branches"=>[], "bound_value"=>0, "attacker_budget"=>K) ### CREATE INITIAL INCUMBENT SOLN
    IncX = feasXs[1]
    nXs = 1
    bestUBVal=1e20

    nNodes=1
    maxidx=1
    maxval=-1
    E=enumerate(BnBTree)
    while true
        println("There are ",length(E)," nodes left out of a total of ", nNodes," generated.")
        currNode = pop!(BnBTree,nodekey)
        println("Current node info: ",currNode)
        if currNode["bound_value"] <= IncX["bound_value"]
            println("\tFathoming due to initial testing of bound ",currNode["bound_value"]," <= ",IncX["bound_value"])
        else
            solveNodeMinmaxSP(model_info, currNode)
	        idx=findNextIndex(model_info["branch_ids"],model_info["x_soln"])
            #pm=solveNodeMinmaxSP(pm_data,pm_form,currNode; use_dual=use_dual_minmax)
            println("Node bound value has been updated to: ",currNode["bound_value"])
            if currNode["bound_value"] <= IncX["bound_value"]
                println("\t\tFathoming due to bound ",currNode["bound_value"]," <= ",IncX["bound_value"])
            elseif model_info["x_soln"][idx] == 0 || model_info["x_soln"][idx] == 1
	            println("Fathoming due to optimality. Adding new attack solution: ",currNode["inactive_branches"])		
	            nXs += 1
	            feasXs[nXs] = currNode
	            if IncX["bound_value"] < currNode["bound_value"]
	                IncX = copy(currNode)
	                println("New incumbent solution: ",IncX)
	            end
            else 
	            # APPLY A PRIMAL HEURISTIC
                primal_val=applyPrimalHeuristic(pm_data,pm_form)
	            if IncX["bound_value"] < primal_val
	                nXs += 1
	                feasXs[nXs] = copy(currNode)
	                feasXs[nXs]["bound_value"] = primal_val
	                IncX = feasXs[nXs]
	                println("New incumbent solution: ",IncX)
	            end
	            # BRANCH

	            println("Adding two new nodes, the branching index is: ",idx)
                node0=Dict{String,Any}()
                node0["bound_value"] = currNode["bound_value"]
                node0["inactive_branches"]=Set(currNode["inactive_branches"])
                node0["protected_branches"]=Set(currNode["protected_branches"])
	            node0["attacker_budget"]=currNode["attacker_budget"]
                push!(node0["protected_branches"],idx)

                node1=Dict{String,Any}()
                node1["bound_value"]=currNode["bound_value"]
                node1["inactive_branches"]=Set(currNode["inactive_branches"])
                push!(node1["inactive_branches"],idx)
	            node1["attacker_budget"]=currNode["attacker_budget"] - 1 ### This node will have another branch inactive, debiting from the inherited budget
	            if node1["attacker_budget"] > 0
                    node1["protected_branches"]=Set(currNode["protected_branches"])
	            else
	                node1["protected_branches"]=setdiff(model_info["branch_ids"],node1["inactive_branches"])
	            end
                println(node0)
                println(node1)
	            BnBTree[nNodes]=node0
                nNodes += 1
	            BnBTree[nNodes]=node1
                nNodes += 1
            end
        end ### if initial fathoming due to bound  

        E=enumerate(BnBTree)
        if length(E) > 0
            nodekey,bestUBVal=findNextNode(E)
        else
            bestUBVal = IncX["bound_value"]
            break
        end
#=
        if (time_ns()-start_time)/1e9 > MAX_TIME
            break
        end
=#
    println("\t\t\tBest UB $bestUBVal versus incumbent value ",IncX["bound_value"])
    end # while
    end_time = time_ns()

    runtime = (end_time-start_time)/1e9
    println("Final best solution: ",IncX)
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
