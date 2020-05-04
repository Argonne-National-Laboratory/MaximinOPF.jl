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

global MAX_TIME = 24*60*60
#global MAX_TIME = 180
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
function print_node(node_key::String, node::Dict{String,Any}; flag_str="",io=Base.stdout)
        attack_str=string("Inactive:",node["inactive_branches"],"\tProtected:",node["protected_branches"])
        println(io,"\n",attack_str)
        acsdp_diff = round(node["bound_value"] - node["sdp_value"];digits=3)
        diff_str=""
        if acsdp_diff > 0
            diff_str=string("ac-sdp=",acsdp_diff)
        else
            diff_str="ac=sdp"
        end
        diff_str=pad_string(diff_str,16)
        gap_str=string("[",node["sdp_value"],
            ",", node["bound_value"],"]")
        gap_str=pad_string(gap_str,16)
        println(io,"\t\t",gap_str,"\t",diff_str,"\t",flag_str)
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
        node_info["value"] = round( max(JuMP.objective_value(model),0.0);digits=3)
        node_info["sdp_value"] = node_info["value"]
    elseif pm_form==ACRPowerModel
        acceptable_status = node_info["solve_status"]==LOCALLY_SOLVED || node_info["solve_status"]==ALMOST_LOCALLY_SOLVED 
        for nn=1:100
            if acceptable_status 
                node_info["value"] = round( max(JuMP.objective_value(model),0.0);digits=3)
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

    node_info["x_soln"]=Dict{Int64,Float64}()
    best_dual_val=-1.0
    test_val=-1.0
    node_info["branch_idx"]=-1
    #println("dual values:")
    for l in pm_data["branch_ids"]
        if l in pm_data["inactive_branches"]
            node_info["x_soln"][l]=1
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
    BnBTree = Dict{String,Dict{String,Any}}()
    BnBTree[""] = Dict("inactive_branches"=>Int64[], "protected_branches"=>Int64[], 
                    "bound_value"=>1e20, "K"=>K, "sdp_value"=>0.0, "bound_value"=>1e20) ### CREATE ROOT NODE
    F_OPT = Dict{String,Dict{String,Any}}()
    F_BD = Dict{String,Dict{String,Any}}()

    treeUB=1e20
    nNodes=1

    BnBTree[""]=solveNodeSP(pm_data, BnBTree[""]; pm_form=SparseSDPWRMPowerModel, fix_all_x=true )
    n_sdp_solves=1

    treeLB = BnBTree[""]["sdp_value"]
    treeLB_soln = ""
    incumbent_node=BnBTree[""]

    while length(BnBTree) > 0
        x_soln_str,treeUB=findNextNode(BnBTree)
        node_soln_str=string(x_soln_str)
        node_soln_str=pad_string(node_soln_str,K*6)
        noteworthy_str=""
        end_time = time_ns()
        runtime = round((end_time-start_time)/1e9)
        if treeUB - treeLB <= eps_tol || runtime >= MAX_TIME
            if treeUB - treeLB <= eps_tol 
                println("Optimal solution found within tolerance $eps_tol, terminating.")
            else
                println("Maximum runtime of $MAX_TIME seconds reached, terminating.")
            end
            break
        end

        BnBTree[x_soln_str]=solveNodeSP(pm_data, BnBTree[x_soln_str]; pm_form=ACRPowerModel, fix_all_x=false, br_idx_rule="max_dual")
	    node_gap_str=string("[", BnBTree[x_soln_str]["sdp_value"],",",BnBTree[x_soln_str]["bound_value"],"]")
        node_gap_str=pad_string(node_gap_str,20)

        fathom_by_opt = (K == length(BnBTree[x_soln_str]["inactive_branches"]))
        fathom_by_opt = fathom_by_opt || length(BnBTree[x_soln_str]["protected_branches"]) + length(BnBTree[x_soln_str]["inactive_branches"]) == length(pm_data["branch_ids"])
        fathom_by_opt = fathom_by_opt || (BnBTree[x_soln_str]["bound_value"] - BnBTree[x_soln_str]["sdp_value"] <= eps_tol)
        if fathom_by_opt
	        noteworthy_str = string(noteworthy_str," ***Fathoming by optimality***")
            F_OPT[x_soln_str]=BnBTree[x_soln_str]
            delete!(BnBTree,x_soln_str)
        elseif BnBTree[x_soln_str]["bound_value"] <= treeLB 
	        noteworthy_str = string(noteworthy_str," ***Fathoming by treeLB $treeLB***")
            F_BD[x_soln_str]=BnBTree[x_soln_str]
            delete!(BnBTree,x_soln_str)
        else
	        # "BRANCH"
            idx = BnBTree[x_soln_str]["branch_idx"]
            up_br_str = string(x_soln_str," $idx")
            BnBTree[up_br_str] = Dict(  "inactive_branches"=>copy(BnBTree[x_soln_str]["inactive_branches"]), 
                                        "protected_branches"=>copy(BnBTree[x_soln_str]["protected_branches"]),"K"=>K, 
                                        "sdp_value"=>BnBTree[x_soln_str]["sdp_value"],
                                        "bound_value"=>BnBTree[x_soln_str]["bound_value"]) 
            nNodes += 2 ### "BnBTree[x_soln_str] is recycled as the 'down' branch node"

            push!(BnBTree[x_soln_str]["protected_branches"],idx)  ### "BnBTree[x_soln_str] is recycled as the 'down' branch node"
            push!(BnBTree[up_br_str]["inactive_branches"],idx)
            BnBTree[up_br_str]=solveNodeSP(pm_data, BnBTree[up_br_str]; pm_form=SparseSDPWRMPowerModel, fix_all_x=true )
            n_sdp_solves += 1
	        if treeLB < BnBTree[up_br_str]["sdp_value"]
	            treeLB = BnBTree[up_br_str]["sdp_value"]
                treeLB_soln = up_br_str
                incumbent_node=BnBTree[up_br_str]
	            noteworthy_str = string(noteworthy_str," ***New incumbent SDP LB $treeLB****")
                for bb in keys(BnBTree)
                    if BnBTree[bb]["bound_value"] < treeLB
                        F_BD[bb]=BnBTree[bb]
                        delete!(BnBTree,bb)
                    end
                end
	        end
        end
        tree_gap=string("[",treeLB,", ",treeUB, "]")
        tree_gap=pad_string(tree_gap,20)
        tree_report=string("(",length(BnBTree),",",length(F_OPT),",",length(F_BD),",",nNodes,")")
        tree_report=pad_string(tree_report,28)
        println("  Soln:",node_soln_str," Node:",node_gap_str," Tree:",tree_gap,"(ACT,OPT,BD,TOT)=",tree_report,noteworthy_str)
    end # "while BnB tree not empty"

    treeUB = max(treeUB,treeLB)
    end_time = time_ns()
    runtime = round((end_time-start_time)/1e9)
    println("Final best SDP LB solution: ",treeLB_soln," with best SDP LB ",treeLB)

    percent_gap_str=""
    if treeLB > 1e-2
        percent_gap_str=string("\tpercent gap: ",100*(treeUB-treeLB)/treeLB)
    end
    println("\tFinal bound gap: [",treeLB,", ",treeUB,"]",percent_gap_str)

    println("Runtime: ",runtime)
    println("Number of nodes processed: ",nNodes)
    println("Number of feasibility problem solves: ",n_sdp_solves)
    println("case \t budget \t nFPSolves \t nNodes \t sec_per_node \t secs ")
    println(" & $K & $n_sdp_solves & $nNodes & ",round(runtime/nNodes;digits=2),
	" & ",round(runtime),"\\\\ \n")
    #io = open(string("case",case_instance,"K",K,"_acbnb_solns.txt"),"w")
    io=Base.stdout

    println(io,"\nNodes fathomed by optimality: ")
    for nn in keys(F_OPT)
        post_str=""
        if F_OPT[nn]["bound_value"] >= treeLB
            post_str="\txxxx"
        end
        print_node(nn, F_OPT[nn];flag_str=post_str)
    end
    println(io,"\nNodes that were still active at termination: ")
    for nn in keys(BnBTree)
        print_node(nn, BnBTree[nn])
    end
    println(io,"\nNodes that were fathomed due to bound (and may warrant more processing): ")
    for nn in keys(F_BD)
        print_node(nn, F_BD[nn])
    end
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
