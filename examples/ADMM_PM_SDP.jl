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

    global case_instance="30"
    global attack_budget=4
    global form_str="SDP"
    global use_h_relax="no"
    global init_prox_t=1
    global case_spec=false
    global budget_spec=false
    global with_cuts=false
    global update_rho=false
    for aa in 1:length(ARGS)
        global case_instance
        global attack_budget
        global init_prox_t
        global with_cuts
        global update_rho
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
        elseif occursin("--prox_t=",ARGS[aa])
            init_prox_t=parse(Float64,ARGS[aa][(length("--prox_t=")+1):length(ARGS[aa])])
            println("Initial prox parameter set to ",init_prox_t)
        elseif occursin("--use_cuts",ARGS[aa])
            if !occursin("--use_cuts=no",ARGS[aa])
                with_cuts=true
            end
        elseif occursin("--update_rho",ARGS[aa])
            if !occursin("--update_rho=no",ARGS[aa])
                update_rho=true
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
 	    "protected_indices" => []
	)
    println("Testing ADMM for root node relaxation of Maxmin problem: ")
    println(testcase)


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
    JuMP.set_optimizer(psd_base_maxmin,with_optimizer(Ipopt.Optimizer,
        print_level=0,
        jac_d_constant="yes",
        jac_c_constant="yes",
        hessian_constant="yes",
        linear_solver="ma57",
        max_iter=50000,
        tol=1e-10,
        constr_viol_tol=1e-10,
        dual_inf_tol=1e-10,
        compl_inf_tol=1e-10,
        acceptable_tol=1e-6,
        acceptable_constr_viol_tol=1e-6,
        acceptable_dual_inf_tol=1e-6,
        acceptable_compl_inf_tol=1e-6 )
    )


    #JuMP.set_optimizer(psd_base_maxmin,with_optimizer(OSQP.Optimizer,verbose=false,eps_abs=1e-9,max_iter=10000))
    #JuMP.set_optimizer(psd_base_maxmin,with_optimizer(CPLEX.Optimizer))
    #JuMP.set_parameter(psd_base_maxmin,"CPXPARAM_ScreenOutput",0)
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

    solve_PSD_via_ADMM(model_info; 
            max_n_iter=10000, prox_t=init_prox_t, prim_tol=1e-3, dual_tol=1e-3, 
            use_cuts=with_cuts, rescale=update_rho, display_freq=100,io=devnull) 
    

    end_time = time_ns()
    runtime = (end_time-start_time)/1e9

    println("Runtime: ",runtime)

