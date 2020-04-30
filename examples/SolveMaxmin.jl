include("../../MaximinOPF/src/MaximinOPF.jl")
using PowerModels
using JuMP
using Mosek
using MosekTools
PowerModels.silence()

global case_instance="30"
global attack_budget=4
global form_str="SOC"
global relax_int=false
global pm_form=SOCWRConicPowerModel
for aa in 1:length(ARGS)
    global case_instance
    global attack_budget
    global pm_form
    global relax_int
    if occursin("--case=",ARGS[aa])
            case_instance=ARGS[aa][(length("--case=")+1):length(ARGS[aa])]
            println("case being set to ",case_instance)
    elseif occursin("--K=",ARGS[aa])
            attack_budget=parse(Int64,ARGS[aa][(length("--K=")+1):length(ARGS[aa])])
            println("attack budget being set to ",attack_budget)
    elseif occursin("--pm_form=",ARGS[aa])
            form_str=ARGS[aa][(length("--pm_form=")+1):length(ARGS[aa])]
            if occursin("SDP",form_str) || occursin("PSD",form_str)
                pm_form=SparseSDPWRMPowerModel
            elseif occursin("SOC",form_str)
                pm_form=SOCWRConicPowerModel
            end
    elseif occursin("--relax_int=yes",ARGS[aa])
        relax_int=true
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
 
soc_solver=with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0)

if (pm_form in MaximinOPF.sdp_pm) && !relax_int
    println(Base.stderr,"This script does not provide support for solving an SDP Maxmin with integrality enforced.")
else
    start_time = time_ns()
    println("Formulating and solving the form ",pm_form)
    maxmin_model=MaximinOPF.MaximinOPFModel(pm_data, pm_form; enforce_int=!relax_int,rm_therm_line_lim=false)
    JuMP.set_optimizer(maxmin_model,soc_solver)
    JuMP.optimize!(maxmin_model)
    branch_ids=sort(collect(pm_data["undecided_branches"]))
    println("Optimal value using powerform ", pm_form, " is: ",JuMP.objective_value(maxmin_model), " with status ",JuMP.termination_status(maxmin_model))
    global x_soln_str = ""
    x_soln = Dict{Int64,Float64}()
    for l in branch_ids
        global x_soln_str
        x_var = variable_by_name(maxmin_model,"x[$l]_1")
        x_val = JuMP.value(x_var)
        if x_val > 1.0-1.0e-8
            x_val = 1
            x_soln_str = string(x_soln_str," $l")
        elseif x_val < 1.0e-8
	        x_val = 0
        end
        x_soln[l] = x_val
    end
    end_time = time_ns()
    runtime = (end_time-start_time)/1e9
    if relax_int
        println("Final best solution: ")
        for l in branch_ids
            if round(x_soln[l];digits=4) > 0
                print(" ($l,",round(x_soln[l];digits=4),")")
            end
        end
        print("\n")
    else
        println("Final best solution: ",x_soln_str)
    end
    println("\twith value ",JuMP.objective_value(maxmin_model))
    println("Runtime: ",runtime)
end


