include("../../MaximinOPF/src/MaximinOPF.jl")
include("../../MaximinOPF/src/utils.jl")
include("psd_utils.jl")
include("prox_pt_methods_psd.jl")
using JuMP
using PowerModels
using ProxSDP
using Ipopt
using OSQP

#using Mosek
#using MosekTools
#using SCIP
PowerModels.silence()

    testcase = Dict(
	    "file" => "data/case118.m", 
	    "PMOption" => SparseSDPWRMPowerModel,
 	    "name" => "case118K4SDP",  	
 	    "attack_budget" => 4,
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
    base_maxmin = MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false, rm_rsoc=true, rm_therm_line_lim=true)
    branch_ids=sort(collect(pm_data["undecided_branches"]))
    psd_base_maxmin = convertSOCtoPSD(base_maxmin)
    JuMP.set_optimizer(psd_base_maxmin,with_optimizer(Ipopt.Optimizer,print_level=0,linear_solver="ma57"))
    #JuMP.set_optimizer(psd_base_maxmin,with_optimizer(OSQP.Optimizer,verbose=false,eps_abs=1e-9,max_iter=10000))
#=
options = Dict(:verbose => false,
                   :eps_abs => 1e-09,
                   :eps_rel => 1e-09,
                   :check_termination => 1,
                   :polish => false,
                   :max_iter => 4000,
                   :rho => 0.1,
                   :adaptive_rho => false,
                   :warm_start => true)
=#


    model_info=Dict{String,Any}()  ### "Essentially, this is an object type containing needed information to apply ADMM"
    model_info = prepare_psd_model_reformulation(psd_base_maxmin, model_info; io=devnull)

    start_time = time_ns()

    solve_PSD_via_ADMM(model_info; max_n_iter=10000, prox_t=0.01, rescale=true, display_freq=100,io=devnull)

    println("Testing final value with psd var fixed to center")
    fixPSDCtr(model_info)
    try
        JuMP.optimize!( model_info["model"] ) 
    catch exc
        println(exc)
	    println("Catching solve error, breaking from subiteration loop.")
    end
    model_info["opt_val"]=JuMP.value(model_info["model"][:linobj_expr])
    model_info["solve_status"]=JuMP.termination_status(model_info["model"])
    println("Optimal value using powerform ", pm_form, " is: ",model_info["opt_val"], " with status ",model_info["solve_status"])
    

    end_time = time_ns()
    runtime = (end_time-start_time)/1e9

    println("Runtime: ",runtime)

