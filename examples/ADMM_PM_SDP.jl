include("../../MaximinOPF/src/MaximinOPF.jl")
include("../../MaximinOPF/src/utils.jl")
include("psd_utils.jl")
include("prox_pt_methods_psd.jl")
using JuMP
using PowerModels
using ProxSDP
using Ipopt

#using Mosek
#using MosekTools
#using SCIP
PowerModels.silence()

    testcase = Dict(
	    "file" => "data/case9.m", 
	    "PMOption" => SparseSDPWRMPowerModel,
 	    "name" => "case9K3SDP",  	
 	    "attack_budget" => 3,
 	    "inactive_indices" => [],
 	    "protected_indices" => []
	)

### "Initialize power model instance, (relaxed) maxmin problem for identifying critical contingencies"
    pm_data = PowerModels.parse_file(testcase["file"])
    pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
    pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
    pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

    ### "sdp_relax=[SDPWRMPowerModel, SparseSDPWRMPowerModel]"
    pm_form = SparseSDPWRMPowerModel
    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false, rm_rsoc=true, rm_therm_line_lim=false)
    branch_ids=sort(collect(pm_data["undecided_branches"]))
    psd_base_maxmin = convertSOCtoPSD(base_maxmin)
    JuMP.set_optimizer(psd_base_maxmin,with_optimizer(Ipopt.Optimizer))


    model_info=Dict{String,Any}()  ### "Essentially, this is an object type containing needed information to apply ADMM"
    model_info = prepare_psd_model_reformulation(psd_base_maxmin, model_info; io=devnull)

    solve_PSD_via_ADMM(model_info; max_n_iter=3000, prox_t=0.1, io=devnull)

    println("Optimal value using powerform ", pm_form, " is: ",model_info["opt_val"], " with status ",model_info["solve_status"])

