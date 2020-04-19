include("../../MaximinOPF/src/MaximinOPF.jl")
using JuMP
using PowerModels
using SCS
using ProxSDP
PowerModels.silence()

global case_instance="30"
global attack_budget=4
global form_str="SOC"
global pm_form=SOCWRConicPowerModel
for aa in 1:length(ARGS)
    global case_instance
    global attack_budget
    global pm_form
    if occursin("--case=",ARGS[aa])
            case_instance=ARGS[aa][(length("--case=")+1):length(ARGS[aa])]
            println("case being set to ",case_instance)
    elseif occursin("--K=",ARGS[aa])
            attack_budget=parse(Int64,ARGS[aa][(length("--K=")+1):length(ARGS[aa])])
            println("attack budget being set to ",attack_budget)
    elseif occursin("--form_str=",ARGS[aa])
            form_str=ARGS[aa][(length("--form_str=")+1):length(ARGS[aa])]
            if occursin("SDP",form_str) || occursin("PSD",form_str)
                pm_form=SparseSDPWRMPowerModel
            elseif occursin("SOC")
                pm_form=SOCWRConicPowerModel
            end
    else
            println(Base.stderr,"Argument ",ARGS[aa]," not recognized.")
    end
end

testcase = Dict(
    "file" => string("data/case",case_instance,".m"),
    "PMOption" => pm_form,
    "name" => string("case",case_instance,"K",attack_budget,form_str),
    "attack_budget" => attack_budget,
    "inactive_indices" => [],
    "protected_indices" => [],
)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

### "sdp_pm=[ SDPWRMPowerModel, SparseSDPWRMPowerModel ]"

conic_solver=with_optimizer(SCS.Optimizer,verbose=0,max_iters=30000)
sdp_solver=with_optimizer(SCS.Optimizer,verbose=0,max_iters=30000)
sparse_sdp_solver=with_optimizer(SCS.Optimizer,verbose=0,max_iters=30000)
#sdp_solver=with_optimizer(ProxSDP.Optimizer, max_iter=100000, log_verbose=false, initial_beta=0.005, full_eig_decomp=false )
#sparse_sdp_solver=with_optimizer(ProxSDP.Optimizer, max_iter=100000,log_verbose=false,initial_beta=0.005,full_eig_decomp=true)

for pm_form in MaximinOPF.conic_supported_pm
    println("\nFormulating and solving the dual minmax problem for the form ",pm_form)
    dual_minmax_model=MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=false)
    if pm_form == SDPWRMPowerModel
        JuMP.set_optimizer(dual_minmax_model,sdp_solver)
    elseif pm_form == SparseSDPWRMPowerModel
        JuMP.set_optimizer(dual_minmax_model,sparse_sdp_solver)
    else 
        JuMP.set_optimizer(dual_minmax_model,conic_solver)
    end
    JuMP.optimize!(dual_minmax_model)
    println("\tOptimal value using powerform ", pm_form, " is: ",JuMP.objective_value(dual_minmax_model), " with status ",JuMP.termination_status(dual_minmax_model))
end

