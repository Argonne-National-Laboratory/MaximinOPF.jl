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
PowerModels.silence()



function solveMaxminViaOABnB(pm_data,pm_form; use_dual_minmax=true)
    #global MAX_TIME
    K=pm_data["attacker_budget"]
    art_bd=ceil(MaximinOPF.getMaxShedding(pm_data);digits=6)

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

    psd_base_maxmin = convertSOCtoPSD(base_maxmin)
    psd_optimizer=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=1)
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
    ### "CPX_ALG_PRIM=1" "CPX_ALG_DUAL=2"
    JuMP.set_parameter(cp_base_maxmin,"CPXPARAM_LPMethod",1)


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
    solveSP(psd_model_info; fix_x=true, compute_projection=false, compute_psd_dual=true)
    n_sdp_solves = 1
    add_cuts(cp_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
    feasXs[""]["value"] = psd_model_info["opt_val"]

    solveSP(psd_model_info; fix_x=false, compute_projection=false, compute_psd_dual=true)
    n_sdp_solves += 1
    add_cuts(cp_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
    
    #bound_obj(cp_model_info; bd_mag=art_bd)
    #println("Artifical bound on the objective function: ",art_bd)
	add_psd_initial_cuts(cp_model_info;io=devnull)
    ## "Initial solve"
    #solveNodeSP_ECP(cp_model_info, BnBTree[0]; compute_projection=false, compute_psd_dual=false)

    bestUBVal=1e20
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
            tree_report=string(" Nodes left: ",length(BnBTree)," out of ", nNodes,)
            println("\tBound gap: [",round(IncX["value"];digits=5),", ",round(bestUBVal;round=5),"]",incumbent_update,",\t",tree_report)
            continue
        end
        printNode(currNode;pretext="\n",posttext="")
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
            solveSP(psd_model_info; fix_x=true, compute_projection=false, compute_psd_dual=true)
            n_sdp_solves += 1
            add_cuts(cp_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
            feasXs[x_soln_str]["value"] = psd_model_info["opt_val"]
            println("\tNew solution ",x_soln_str," has known optimal value: ",round(feasXs[x_soln_str]["value"];digits=5))
	    if IncX["value"] < feasXs[x_soln_str]["value"]
	    	IncX = feasXs[x_soln_str]
	        incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"]," with value ",round(IncX["value"];digits=5))
                deleteNodesByBound(BnBTree,IncX["value"])
	    end
	end

        while true
            #JuMP.set_start_value.(all_variables(cp_model_info["model"]), value.(all_variables(cp_model_info["model"])))
            solveNodeSP_ECP(cp_model_info, currNode)
            for nn=1:n_extra_cuts
                PSDProjections(cp_model_info)
                add_cuts(cp_model_info, cp_model_info["psd_info"];cut_type="orth_expr_val")
                #add_cuts(psd_model_info, cp_model_info["psd_info"];cut_type="orth_expr_val")
                solveNodeSP_ECP(cp_model_info, currNode; compute_psd_dual=false)
                println("\tNew node bound: ",round(cp_model_info["opt_val"];digits=5), " with status: ", cp_model_info["solve_status"])
            end

            currNode["bound_value"] = min(cp_model_info["opt_val"],currNode["bound_value"])
            if currNode["bound_value"] <= IncX["value"]
                println("\tFathoming node $nodekey due to updated bound ",round(currNode["bound_value"];digits=5)," <= ",round(IncX["value"];digits=5))
                break
            end
            println("\tNew node bound: ",round(cp_model_info["opt_val"];digits=5), 
                " with status: ", cp_model_info["solve_status"])
            idx=findNextIndex(cp_model_info["branch_ids"],cp_model_info["x_soln"])
            if cp_model_info["x_soln"][idx] == 0 || cp_model_info["x_soln"][idx] == 1
                x_soln_str = cp_model_info["x_soln_str"]
                if !haskey(feasXs,x_soln_str)
	            nXs += 1
                    feasXs[x_soln_str]=Dict("bound_value"=>cp_model_info["opt_val"], "value"=>0)
                    feasXs[x_soln_str]["x_soln"] = copy(cp_model_info["x_soln"])
                    feasXs[x_soln_str]["x_soln_str"] = x_soln_str
                    solveSP(psd_model_info; fix_x=true, compute_projection=false, compute_psd_dual=true)
                    n_sdp_solves += 1
                    add_cuts(cp_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
                    #add_cuts(psd_model_info, psd_model_info["psd_info"]; cut_type="dual_val")
                    feasXs[x_soln_str]["value"] = psd_model_info["opt_val"]
                    println("\tNew solution ",x_soln_str," has known optimal value: ",round(feasXs[x_soln_str]["value"];digits=5))
	            if IncX["value"] < feasXs[x_soln_str]["value"]
	                IncX = feasXs[x_soln_str]
	                incumbent_update = string("\t***New incumbent***",IncX["x_soln_str"]," with value ",round(IncX["value"];digits=5))
                        deleteNodesByBound(BnBTree,IncX["value"])
	            end
                    continue
                else
                    println("\tSoln: ",x_soln_str, ": cp_val ",cp_model_info["opt_val"], " vs psd_val: ",feasXs[x_soln_str]["value"])
	                println("\tFathoming due to optimality. ")
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
    println("\tFinal bound gap: [",IncX["value"],", ",bestUBVal,"], percent gap: ",round(100*(bestUBVal-IncX["value"])/IncX["value"];digits=1))
    println("Runtime: ",runtime)
    println("Number of nodes processed: ",nNodes)
    println("Number of feasibility problem solves: ",n_sdp_solves)
    println("case \t budget \t nFPSolves \t nNodes \t sec_per_node \t secs ")
    println(" & $K & $n_sdp_solves & $nNodes & ",round(runtime/nNodes;digits=2),
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

solveMaxminViaOABnB(pm_data,SparseSDPWRMPowerModel,use_dual_minmax=true)
#solveMaxminViaOABnB(pm_data,SDPWRMPowerModel,use_dual_minmax=true)
