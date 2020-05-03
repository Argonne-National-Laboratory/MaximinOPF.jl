include("../../MaximinOPF/src/MaximinOPF.jl")
using PowerModels
using JuMP
using Ipopt
PowerModels.silence()

global case_instance="57"
global attack_budget=0
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
    #"inactive_indices" => [6,7,29,36],  # "For case 30 K=4"
    #"inactive_indices" => [57 60 65 66],  # "For case 57 K=4" 1.611
    #"inactive_indices" => [60 65 66 72],  # "For case 57 K=4" 1.752
    #"inactive_indices" => [58 60 65 66],  # "For case 57 K=4" 1.748
    #"inactive_indices" => [41,60,65,66], # "For case 57 K=4"
    #"inactive_indices" => [41 60 66 80],  # "For case 57 K=4"
    "inactive_indices" => [17 41 60 72],  # "For case 57 K=4"
    "protected_indices" => [8, 15, 59, 66, 58, 65, 27],
)
pm_data = PowerModels.parse_file(testcase["file"])
pm_data["attacker_budget"] = testcase["attack_budget"] ###Adding another key and entry
pm_data["inactive_branches"] = testcase["inactive_indices"] ###Adding another key and entry
pm_data["protected_branches"] = testcase["protected_indices"] ###Adding another key and entry

ip_solver=with_optimizer(Ipopt.Optimizer,
    print_level=0,
    linear_solver="ma57",
    max_iter=50000
)
    #mu_init=0.5,
    #mu_linear_decrease_factor=0.8,
    #mu_strategy="adaptive",
    #ma57_automatic_scaling="no",

# nonconvex AC forms
#nonconvex_ac=[ACPPowerModel, ACRPowerModel, ACTPowerModel]
    pm_form=ACRPowerModel
    println("Formulating and solving the form ",pm_form)
    #pm = MaximinOPF.PF_FeasModel(pm_data, pm_form)
    pm = MaximinOPF.MinimaxOPFModel(pm_data, pm_form)
    JuMP.set_optimizer(pm.model,ip_solver)
    #println("Starting voltage:")

    JuMP.optimize!(pm.model)
    opt_status = JuMP.termination_status(pm.model)
    opt_val = 1.0e20

    if opt_status == LOCALLY_SOLVED || opt_status == ALMOST_LOCALLY_SOLVED
        opt_val = round(JuMP.objective_value(pm.model);digits=2)
    end
    println("\tInit: Optimal value using powerform ", pm_form, " is: ",opt_val, " with status ",opt_status)


    loc_opt=Dict{Float64,Int64}(1.0e20=>0)
    for nn=1:100
        for ii in ids(pm,:bus)
            vr=JuMP.variable_by_name(pm.model,"0_vr[$ii]")
            vi=JuMP.variable_by_name(pm.model,"0_vi[$ii]")
            vr_start = 1
            vi_start = rand()-0.5
            scal=1.0
            if vr_start^2 + vi_start^2 > 0
                scal = 1.0/sqrt(vr_start^2+vi_start^2)
            end
            JuMP.set_start_value(vr,scal*vr_start)
            JuMP.set_start_value(vi,scal*vi_start)
            #print(" ($vr_start,$vi_start)")
        end
        #print("\n")
        JuMP.optimize!(pm.model)
        opt_status = JuMP.termination_status(pm.model)
        opt_val = 1.0e20

        if opt_status == LOCALLY_SOLVED || opt_status == ALMOST_LOCALLY_SOLVED
            opt_val = round(JuMP.objective_value(pm.model);digits=2)
        end
        println("\t$nn: Optimal value using powerform ", pm_form, " is: ",opt_val, " with status ",opt_status)
        if haskey(loc_opt,opt_val)
            loc_opt[opt_val] += 1
        else
            loc_opt[opt_val] = 1
        end
    end

    loc_opt_vals=sort(collect(keys(loc_opt)))
    println("Histogram for local optima: ")
    for vv in loc_opt_vals
        println("\t",vv," => ",loc_opt[vv])
    end

