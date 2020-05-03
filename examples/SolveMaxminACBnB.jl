#=
May 2020
Branch-and-bound approach for solving the AC powerflow network maxmin vulnerability identification problem
Paper forthcoming with authors
Brian Dandurand
Kibaek Kim
Sven Leyffer
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

#global MAX_TIME = 24*60*60
global MAX_TIME = 300
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
function solveNodeSP(pm_data, node_info; pm_form=ACRPowerModel, fix_all_x=false, br_idx_rule="max_dual" )
    global ac_optimizer, psd_optimizer
    pm_data["inactive_branches"]=sort(node_info["inactive_branches"])
    if fix_all_x
        pm_data["protected_branches"] = getNonInactiveBranches(pm_data)
    else
        pm_data["protected_branches"]=sort(node_info["protected_branches"])
    end
    pm_data["undecided_branches"] = getUndecidedBranches(pm_data)

    pm_data["attacker_budget"] = max( node_info["K"] - length(pm_data["inactive_branches"]), 0)

    model_pm = MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
    model = model_pm.model
    if pm_form==ACRPowerModel
        JuMP.set_optimizer(model,ac_optimizer)
    elseif pm_form==SparseSDPWRMPowerModel || pm_form==SDPWRMPowerModel
        JuMP.set_optimizer(model,psd_optimizer)
    end

    JuMP.optimize!(model)
    node_info["value"] = 1.0e30
    node_info["solve_status"] = JuMP.termination_status(model)
    if pm_form==SparseSDPWRMPowerModel || pm_form==SDPWRMPowerModel
        node_info["value"] = round(JuMP.objective_value(model);digits=3)
        node_info["sdp_value"] = node_info["value"]
    elseif pm_form==ACRPowerModel
        acceptable_status = node_info["solve_status"]==LOCALLY_SOLVED || node_info["solve_status"]==ALMOST_LOCALLY_SOLVED 
        for nn=1:100
            if acceptable_status 
                node_info["value"] = round(JuMP.objective_value(model);digits=3)
                if node_info["value"] <= node_info["bound_value"]
                    break
                end
            end
            println("FLAGGING: Re-solving AC problem due to unsatisfactory status or local opt value....")
            for ii in ids(model_pm,:bus)
                vr=JuMP.variable_by_name(model,"0_vr[$ii]")
                vi=JuMP.variable_by_name(model,"0_vi[$ii]")
                vr_start = 1
                vi_start = rand()-0.5
                scal=1.0
                if vr_start^2 + vi_start^2 > 0
                    scal = 1.0/sqrt(vr_start^2+vi_start^2)
                end
                JuMP.set_start_value(vr,scal*vr_start)
                JuMP.set_start_value(vi,scal*vi_start)
            end
            JuMP.optimize!(model)
            node_info["solve_status"] = JuMP.termination_status(model)
            acceptable_status = node_info["solve_status"]==LOCALLY_SOLVED || node_info["solve_status"]==ALMOST_LOCALLY_SOLVED 
        end
        if !acceptable_status
            println("WARNING: satisfactory status never achieved for the AC problem.")
        end
        node_info["ac_value"] = min(node_info["value"],node_info["bound_value"])
    end

    if !fix_all_x
        node_info["bound_value"] = min(node_info["value"],node_info["bound_value"])
        node_info["uK"]=JuMP.value(var(model_pm,model_pm.cnw)[:u_K])
        node_info["u"]=Dict{Int64,Float64}()
    end

    node_info["x_soln_str"]=""
    node_info["x_soln"]=Dict{Int64,Float64}()
    best_dual_val=-1.0
    test_val=-1.0
    node_info["branch_idx"]=-1
    #println("dual values:")
    for l in pm_data["branch_ids"]
        if l in pm_data["inactive_branches"]
            node_info["x_soln"][l]=1
            node_info["x_soln_str"] = string(node_info["x_soln_str"]," $l")
        elseif l in pm_data["protected_branches"] 
            node_info["x_soln"][l]=0
        else
            @assert !fix_all_x
            node_info["x_soln"][l]=JuMP.dual(JuMP.constraint_by_name(model,"x[$l]"))

            pd_br = var(model_pm,model_pm.cnw,:pd_br)
            qd_br = var(model_pm,model_pm.cnw,:qd_br)
            u_ord_aux = var(model_pm, model_pm.cnw)[:u_ord_aux][l] 
            l_arcs = filter(a->(a[1] == l),ref(model_pm,model_pm.cnw,:arcs))
            obj_l_term = u_ord_aux + sum( pd_br[a,0] + qd_br[a,0] + pd_br[a,1] + qd_br[a,1] for a in l_arcs)
            node_info["u"][l]=JuMP.value(obj_l_term)

            if node_info["x_soln"][l] > 1.0
                node_info["x_soln"][l]=1
            elseif node_info["x_soln"][l] < 0.0
                node_info["x_soln"][l]=0
            end
            if br_idx_rule=="max_frac" 
                test_val = node_info["x_soln"][l]*(1-node_info["x_soln"][l])
            elseif br_idx_rule=="max_dual" 
                test_val = node_info["x_soln"][l]
            elseif br_idx_rule=="max_obj_term" 
                test_val = node_info["u"][l]
            end
            if best_dual_val < test_val
                best_dual_val = test_val
                node_info["branch_idx"]=l
            end
        end
    end
    return node_info
end


function solveMaxminACBnB(pm_data)
    global MAX_TIME
    K=pm_data["attacker_budget"]
    @assert !haskey(pm_data,"branch_ids")
    pm_data["branch_ids"] = MaximinOPF.getBranchIds(pm_data)
    eps_tol=1e-2

    start_time = time_ns()


    feasXs = Dict{String,Dict{String,Any}}()
    BnBTree = Dict{Int64,Dict{String,Any}}()
    BnBTree[0] = Dict("inactive_branches"=>[], "protected_branches"=>[], 
                    "bound_value"=>1e20, "K"=>K, "attacker_budget"=>K, "node_key"=>0) ### CREATE ROOT NODE
    treeUB=1e20
	nXs = 1
    BnBTree[0]=solveNodeSP(pm_data, BnBTree[0]; pm_form=SparseSDPWRMPowerModel, fix_all_x=true )
    #solveNodeSP(pm_data, BnBTree[0]; pm_form=ACRPowerModel, fix_all_x=true )
    n_sdp_solves=1

    treeLB = BnBTree[0]["sdp_value"]
    treeLB_soln = BnBTree[0]["x_soln_str"]
    feasXs[treeLB_soln]=Dict{String,Any}("sdp_value"=>treeLB,"bound_value"=>treeUB,"ac_value"=>treeUB, "status"=>"ACTIVE")
    println("\tInitial solution [",treeLB_soln,"] has known optimal value SDP lower bound ", treeLB)

    nNodes=1
    while length(BnBTree) > 0
        nodekey,treeUB=findNextNode(BnBTree)
        incumbent_update=""
        currNode = pop!(BnBTree,nodekey)
        end_time = time_ns()
        runtime = round((end_time-start_time)/1e9)
        if treeUB - treeLB <= eps_tol || runtime >= MAX_TIME
            if treeUB - treeLB <= eps_tol 
                println("Optimal solution found within tolerance $eps_tol, terminating.")
                if !haskey(feasXs,treeLB_soln)
                    feasXs[treeLB_soln]=Dict{String,Any}("bound_value"=>treeUB,"sdp_value"=>treeLB, "ac_value"=>treeUB, "x_soln_str"=>treeLB_soln, "status"=>"F_OPT")
                end
            else
                println("Maximum runtime of $MAX_TIME seconds reached, terminating.")
            end
            #currNode=solveNodeSP(pm_data, currNode; pm_form=ACRPowerModel, fix_all_x=true )
            #feasXs[currNode["x_soln_str"] ]["ac_value"] = currNode["ac_value"]
            #feasXs[currNode["x_soln_str"] ]["bound_value"] = currNode["ac_value"]
            #feasXs[ currNode["x_soln_str"] ]["status"]="F_TERM"
            break
        end
        printNode(currNode;pretext="\n",posttext="")

        currNode=solveNodeSP(pm_data, currNode; pm_form=ACRPowerModel, fix_all_x=false, br_idx_rule="max_dual")
        feasXs[currNode["x_soln_str"] ]["bound_value"] = currNode["bound_value"]
        feasXs[currNode["x_soln_str"] ]["ac_value"] = currNode["ac_value"]
        println("\tNew node bound: ",currNode["bound_value"], " with status: ", currNode["solve_status"])

        fathom_by_opt = (K == length(currNode["inactive_branches"]))
        fathom_by_opt = fathom_by_opt || length(currNode["protected_branches"]) + length(currNode["inactive_branches"]) == length(pm_data["branch_ids"])
        fathom_by_opt = fathom_by_opt || (currNode["bound_value"] - currNode["sdp_value"] <= eps_tol)
        if fathom_by_opt
	        println("\tFathoming node $nodekey due to optimality.")
            println("\tFathomed solution ", currNode["x_soln_str"]," has known optimal value bounded between: [",
                feasXs[ currNode["x_soln_str"] ]["sdp_value"],",",
                feasXs[ currNode["x_soln_str"] ]["ac_value"],",",
                feasXs[ currNode["x_soln_str"] ]["bound_value"],"]"
            )
            feasXs[ currNode["x_soln_str"] ]["status"]="F_OPT"
        elseif currNode["bound_value"] <= treeLB 
            println("\tFathoming node $nodekey due to its updated bound ",currNode["bound_value"]," <= ",treeLB)
            feasXs[ currNode["x_soln_str"] ]["status"]="F_BD"
        else
	        # "BRANCH"
            idx=currNode["branch_idx"]
            bd_val=currNode["bound_value"]
	        println("\tAdding two new nodes, the branching index is: ",idx," based on dual val: ",currNode["x_soln"][idx])
            node0=Dict{String,Any}("node_key"=>nNodes, "bound_value"=>bd_val, "K"=>K,"x_soln_str"=>currNode["x_soln_str"],
                "sdp_value"=>currNode["sdp_value"], "ac_value"=>currNode["ac_value"], 
                "inactive_branches"=>copy(currNode["inactive_branches"]), 
                "protected_branches"=>copy(currNode["protected_branches"])
            )
            push!(node0["protected_branches"],idx)
	        BnBTree[nNodes]=node0
            nNodes += 1

            node1=Dict{String,Any}("node_key"=>nNodes, "bound_value"=>bd_val,"K"=>K,"x_soln_str"=>currNode["x_soln_str"],
                "sdp_value"=>currNode["sdp_value"], "ac_value"=>currNode["ac_value"], 
                "inactive_branches"=>copy(currNode["inactive_branches"]), 
                "protected_branches"=>copy(currNode["protected_branches"])
            )
            push!(node1["inactive_branches"],idx)
            node1=solveNodeSP(pm_data, node1; pm_form=SparseSDPWRMPowerModel, fix_all_x=true )
            n_sdp_solves += 1
            feasXs[node1["x_soln_str"]]=Dict{String,Any}("sdp_value"=>node1["sdp_value"],"bound_value"=>node1["bound_value"],
                                                         "ac_value"=>node1["ac_value"], "status"=>"ACTIVE")
            println("\tSDP solve no. $n_sdp_solves, soln: ",node1["x_soln_str"],
                " has SDP lower bound value ", node1["sdp_value"])
	        if treeLB < node1["sdp_value"]
	            treeLB = node1["sdp_value"]
                treeLB_soln = node1["x_soln_str"]
	            incumbent_update = string("\t***New incumbent SDP LB***",node1["x_soln_str"]," with value ",treeLB)
                for bb in keys(BnBTree)
                    if BnBTree[bb]["bound_value"] < treeLB
                        feasXs[ BnBTree[bb]["x_soln_str"] ]["status"]="F_BD"
                        delete!(BnBTree,bb)
                    end
                end
	        end
            BnBTree[nNodes]=node1
            nNodes += 1
        end
        tree_report=string(" Nodes left: ",length(BnBTree)," out of ", nNodes,)
        println("\tBound gap: [",treeLB,", ",treeUB, "]",incumbent_update,",\t",tree_report)
    end # "while BnB tree not empty"

    treeUB = treeLB
    end_time = time_ns()
    runtime = round((end_time-start_time)/1e9)
    println("Final best SDP LB solution: ",treeLB_soln," with best SDP LB ",treeLB)
#=
    percent_gap_str=""
    if treeLB > 1e-2
        percent_gap_str=string("\tpercent gap: ",100*(treeUB-treeLB)/treeLB)
    end
    println("\tFinal bound gap: [",treeLB,", ",treeUB,"]",percent_gap_str)
=#
    println("Runtime: ",runtime)
    println("Number of nodes processed: ",nNodes)
    println("Number of feasibility problem solves: ",n_sdp_solves)
    println("case \t budget \t nFPSolves \t nNodes \t sec_per_node \t secs ")
    println(" & $K & $n_sdp_solves & $nNodes & ",round(runtime/nNodes;digits=2),
	" & ",round(runtime),"\\\\ \n")
    #io = open(string("case",case_instance,"K",K,"_acbnb_solns.txt"),"w")
    io=Base.stdout
    println(io,"Feasible X's and their value bounds:")
#=
    for kk in keys(BnBTree)
        BnBTree[kk]=solveNodeSP(pm_data, BnBTree[kk]; pm_form=ACRPowerModel, fix_all_x=true )
        feasXs[ BnBTree[kk]["x_soln_str"] ]["ac_value"] = BnBTree[kk]["ac_value"]
        feasXs[ BnBTree[kk]["x_soln_str"] ]["bound_value"] = BnBTree[kk]["ac_value"]
        feasXs[ BnBTree[kk]["x_soln_str"] ]["status"]="F_TERM"
    end
=#
    for nn in keys(feasXs)
        #if feasXs[nn]["status"] !="F_BD" 
        if true 
            attack_str=nn
            attack_str=pad_string(attack_str,16)
            acsdp_diff = round(feasXs[nn]["bound_value"] - feasXs[nn]["sdp_value"];digits=3)
            diff_str=""
            if acsdp_diff > 0
                diff_str=string("ac-sdp=",acsdp_diff)
            else
                diff_str="ac=sdp"
            end
            diff_str=pad_string(diff_str,16)
            flag_str=feasXs[nn]["status"]
            if feasXs[nn]["bound_value"] > treeLB
                flag_str=string(flag_str,"\t****")
            end
            gap_str=string("[",feasXs[nn]["sdp_value"],
                ",",feasXs[nn]["bound_value"],"]")
            gap_str=pad_string(gap_str,32)
            println(io,attack_str,"\t",gap_str,"\t",diff_str,"\t",flag_str)
        end
    end
    #close(io)
    return treeUB,nNodes,runtime
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
println("Testing SolveMaxminACBnB toward solving the AC Maxmin problem: ")
println(testcase)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry
pm_data["use_rf"]=testcase["use_rf"]

solveMaxminACBnB(pm_data)
