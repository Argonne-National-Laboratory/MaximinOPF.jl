#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("../../MaximinOPF/src/MaximinOPF.jl")
include("../../MaximinOPF/src/utils.jl")
include("bnb_utils.jl")
include("psd_utils.jl")
using PowerModels
using Ipopt
using Mosek
using MosekTools
#using CPLEX
#using SCIP
PowerModels.silence()

global ac_optimizer=with_optimizer(Ipopt.Optimizer,
    print_level=0,
    linear_solver="ma57",
    max_iter=50000
)
global psd_optimizer=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=1)

function pad_string(str,len)
    while length(str)<len
        str=string(str," ")
    end
    return str
end

#  "nonconvex_ac=[ACPPowerModel, ACRPowerModel, ACTPowerModel]"
function solveNodeSP(pm_data, node_info; pm_form=ACRPowerModel, fix_all_x=false )
    global ac_optimizer, psd_optimizer
    pm_data["inactive_branches"]=node_info["inactive_branches"]
    if fix_all_x
        pm_data["protected_branches"] = filter(l->!(l in pm_data["inactive_branches"]), pm_data["branch_ids"])
    else
        pm_data["protected_branches"]=node_info["protected_branches"]
    end
    pm_data["undecided_branches"] = filter(l->!( (l in pm_data["inactive_branches"]) || (l in pm_data["protected_branches"]) ), pm_data["branch_ids"])
    pm_data["attacker_budget"] = min(node_info["K"] - length(pm_data["inactive_branches"]), length(pm_data["undecided_branches"]))

    model_pm = MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
    model = model_pm.model
    if pm_form==ACRPowerModel
        JuMP.set_optimizer(model,ac_optimizer)
    elseif pm_form==SparseSDPWRMPowerModel || pm_form==SDPWRMPowerModel
        JuMP.set_optimizer(model,psd_optimizer)
    end

    JuMP.optimize!(model)
    node_info["value"]=round(JuMP.objective_value(model);digits=3)
    if pm_form==ACRPowerModel
        node_info["ac_value"]=node_info["value"]
        if !fix_all_x
            node_info["bound_value"] = min(node_info["value"],node_info["bound_value"])
        end
    elseif pm_form==SparseSDPWRMPowerModel || pm_form==SDPWRMPowerModel
        node_info["sdp_value"]=node_info["value"]
    end
    node_info["solve_status"]=JuMP.termination_status(model)

    node_info["x_soln_str"]=""
    node_info["x_soln"]=Dict{Int64,Float64}()
    for l in pm_data["branch_ids"]
        if l in node_info["inactive_branches"]
            node_info["x_soln"][l]=1
            node_info["x_soln_str"] = string(node_info["x_soln_str"]," $l")
        elseif l in node_info["protected_branches"] || fix_all_x
            node_info["x_soln"][l]=0
        else
            @assert l in pm_data["undecided_branches"]
            node_info["x_soln"][l]=JuMP.dual(JuMP.constraint_by_name(model,"x[$l]"))
            if node_info["x_soln"][l] > 1.0-1.0e-5
                node_info["x_soln"][l]=1
                node_info["x_soln_str"] = string(node_info["x_soln_str"]," $l")
            elseif node_info["x_soln"][l] < 1.0e-5
                node_info["x_soln"][l]=0
            end
        end
    end
    return node_info
end


function solveMaxminACBnB(pm_data)
    #global MAX_TIME
    K=pm_data["attacker_budget"]
    @assert !haskey(pm_data,"branch_ids")
    pm_data["branch_ids"] = MaximinOPF.getBranchIds(pm_data)
    art_bd=ceil(MaximinOPF.getMaxShedding(pm_data);digits=6)

    start_time = time_ns()


    BnBTree = Dict{Int64,Dict{String,Any}}()
    BnBTree[0] = Dict("inactive_branches"=>[], "protected_branches"=>[], 
                    "bound_value"=>1e20, "K"=>K, "attacker_budget"=>K, "node_key"=>0) ### CREATE ROOT NODE
    bestUBVal=1e20

    feasXs = Dict{String,Dict{String,Any}}()
    feasXs[""]=Dict("x_soln"=>Dict{Int64,Float64}(),"x_soln_str"=>"", "bound_value"=>0, "value"=>0) ### CREATE INITIAL INCUMBENT SOLN
    IncX = feasXs[""]
    nXs = 1
    #println("Initial processing of root node yields an upper bound of : ",round(bestUBVal;digits=5))

    nNodes=1
    maxidx=1
    maxval=-1
    n_extra_cuts = 0
    while length(BnBTree) > 0
        nodekey,bestUBVal=findNextNode(BnBTree)
        incumbent_update=""
        currNode = pop!(BnBTree,nodekey)
        if currNode["bound_value"] <= IncX["value"]
            println("\nFathoming node $nodekey due to initial bound ",round(currNode["bound_value"];digits=5)," <= ",round(IncX["value"];digits=5))
            tree_report=string(" Nodes left: ",length(BnBTree)," out of ", nNodes)
            println("\tBound gap: [",round(IncX["value"];digits=5),", ",round(bestUBVal;round=5),"]",incumbent_update,",\t",tree_report)
            continue
        end
        printNode(currNode;pretext="\n",posttext="")
	    x_soln_str=""
    	cut_lines=sort(collect(currNode["inactive_branches"]))
    	for l in cut_lines
            x_soln_str=string(x_soln_str," $l")
    	end
        if !haskey(feasXs,x_soln_str)
	        nXs += 1
            feasXs[x_soln_str]=Dict{String,Any}("K"=>K,"x_soln"=>Dict{Int64,Float64}(),"x_soln_str"=>x_soln_str)
            solveNodeSP(pm_data, currNode; pm_form=ACRPowerModel, fix_all_x=true )
            solveNodeSP(pm_data, currNode; pm_form=SparseSDPWRMPowerModel, fix_all_x=true )
            feasXs[x_soln_str]=copy(currNode)
            println("\tNew solution ",x_soln_str," has known optimal value bounded between: [",
                feasXs[x_soln_str]["sdp_value"],",",
                feasXs[x_soln_str]["ac_value"],"]"
            )
	        if IncX["value"] < feasXs[x_soln_str]["sdp_value"]
	    	    IncX = feasXs[x_soln_str]
	            incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"]," with value ",IncX["value"])
                deleteNodesByBound(BnBTree,IncX["value"])
	        end
	    end

        while true
            solveNodeSP(pm_data, currNode; pm_form=ACRPowerModel, fix_all_x=false )

            if currNode["bound_value"] <= IncX["value"]
                println("\tFathoming node $nodekey due to updated bound ",currNode["bound_value"]," <= ",IncX["value"])
                break
            end
            println("\tNew node bound: ",currNode["bound_value"], 
                " with status: ", currNode["solve_status"])
            idx=findNextIndex(pm_data["branch_ids"],currNode["x_soln"])
            if currNode["x_soln"][idx] == 0 || currNode["x_soln"][idx] == 1
                x_soln_str = currNode["x_soln_str"]
                if !haskey(feasXs,x_soln_str)
                    solveNodeSP(pm_data, currNode; pm_form=SparseSDPWRMPowerModel, fix_all_x=false )
	                nXs += 1
                    feasXs[x_soln_str]=copy(currNode)
                    println("\tNew solution ",x_soln_str," has known optimal value bounded between: [",
                        feasXs[x_soln_str]["sdp_value"],",",
                        feasXs[x_soln_str]["ac_value"],"]"
                    )
	                if IncX["value"] < feasXs[x_soln_str]["value"]
	                    IncX = feasXs[x_soln_str]
	                    incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"]," with value ",IncX["value"])
                        deleteNodesByBound(BnBTree,IncX["value"])
	                end
                    continue
                else
	                println("\tFathoming due to optimality. ")
                    break
                end
            else
	            # "BRANCH"
                bd_val=currNode["bound_value"]
	            println("\tAdding two new nodes, the branching index is: ",idx)
                node0=Dict{String,Any}("node_key"=>nNodes, "bound_value"=>bd_val,"K"=>K,
                    "inactive_branches"=>copy(currNode["inactive_branches"]), 
                    "protected_branches"=>copy(currNode["protected_branches"])
                )
                push!(node0["protected_branches"],idx)
	            BnBTree[nNodes]=node0
                nNodes += 1

                node1=Dict{String,Any}("node_key"=>nNodes, "bound_value"=>bd_val,"K"=>K,
                    "inactive_branches"=>copy(currNode["inactive_branches"]), 
                    "protected_branches"=>copy(currNode["protected_branches"])
                )
                push!(node1["inactive_branches"],idx)
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
        println("\tBound gap: [",IncX["value"],", ",bestUBVal,
            "]",incumbent_update,",\t",tree_report)
    end # "while BnB tree not empty"
    bestUBVal = IncX["value"]
    end_time = time_ns()

    runtime = (end_time-start_time)/1e9
    println("Final best solution: ",IncX["x_soln_str"]," with value ",IncX["value"])
    println("\tFinal bound gap: [",IncX["value"],", ",bestUBVal,"], percent gap: ",100*(bestUBVal-IncX["value"])/IncX["value"])
    println("Runtime: ",runtime)
    println("Number of nodes processed: ",nNodes)
    #println("Number of feasibility problem solves: ",n_sdp_solves)
    println("case \t budget \t nFPSolves \t nNodes \t sec_per_node \t secs ")
    #println(" & $K & $n_sdp_solves & $nNodes & ",round(runtime/nNodes;digits=2),
	#" & ",round(runtime),"\\\\ \n")
    io = open(string("case",case_instance,"_acbnb_solns.txt"),"w")
    println(io,"Feasible X's and their value bounds:")
    for nn in keys(feasXs)
        if feasXs[nn]["ac_value"] > 0
            attack_str=string(feasXs[nn]["x_soln_str"])
            attack_str=pad_string(attack_str,16)
            acsdp_diff = round(feasXs[nn]["ac_value"] - feasXs[nn]["sdp_value"];digits=3)
            diff_str=""
            if acsdp_diff > 0
                diff_str=string("ac-sdp=",acsdp_diff)
            else
                diff_str="ac=sdp"
            end
            diff_str=pad_string(diff_str,16)
            flag_str=""
            if feasXs[nn]["ac_value"] > IncX["value"]
                flag_str="****"
            end
            gap_str=string("[",feasXs[nn]["sdp_value"],",",feasXs[nn]["ac_value"],"]")
            gap_str=pad_string(gap_str,16)
            println(io,attack_str,"\t",gap_str,"\t",diff_str,"\t",flag_str)
        end
    end
    close(io)
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
    elseif occursin("--form_str=",ARGS[aa])
            form_str=ARGS[aa][(length("--form_str=")+1):length(ARGS[aa])]
    elseif occursin("--use_rf",ARGS[aa])
    	if occursin("=no",ARGS[aa]) || occursin("=false",ARGS[aa])
            use_rf_relax=false
	else
            use_rf_relax=true
	end
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
println("Testing SolveMaxminViaOABnB for root node relaxation of Maxmin problem: ")
println(testcase)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry
pm_data["use_rf"]=testcase["use_rf"]

solveMaxminACBnB(pm_data)
#solveMaxminViaOABnB(pm_data,SDPWRMPowerModel,use_dual_minmax=true)
