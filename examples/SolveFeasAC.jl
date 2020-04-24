include("../../MaximinOPF/src/MaximinOPF.jl")
using PowerModels
using JuMP
using Ipopt
PowerModels.silence()

global case_instance="30"
global attack_budget=4
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
    else
            println(Base.stderr,"Argument ",ARGS[aa]," not recognized.")
    end
end

testcase = Dict(
    "file" => string("data/case",case_instance,".m"),
    "name" => string("case",case_instance,"K",attack_budget,"AC"),
    "attack_budget" => attack_budget,
    "inactive_indices" => [6,7,29,36],
    "protected_indices" => [],
)

pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

ip_solver=with_optimizer(Ipopt.Optimizer,
    print_level=0,
    linear_solver="ma57",
    ma57_automatic_scaling="yes",
    max_iter=50000
)
#mu_strategy="adaptive",

# nonconvex AC forms
nonconvex_ac=[ACPPowerModel, ACRPowerModel, ACTPowerModel]
for pm_form in nonconvex_ac
    println("Formulating and solving the form ",pm_form)
    pm = MaximinOPF.PF_FeasModel(pm_data, pm_form)
    JuMP.set_optimizer(pm.model,ip_solver)
    JuMP.optimize!(pm.model)
    println("\tOptimal value using powerform ", pm_form, " is: ",JuMP.objective_value(pm.model), "with status ",JuMP.termination_status(pm.model))
end

