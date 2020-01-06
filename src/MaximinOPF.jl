module MaximinOPF
using PowerModels
using Dualization
using MathOptInterface
include("Variables.jl")
include("Objectives.jl")
include("Constraints.jl")
greet() = print("Hello World!")

function MaximinOPFModel(pm_data, powerform)
    println("Hello MaximinOPFModel")
    m = MinimaxOPFModel(pm_data, powerform)
    m = DualizeModel(m)
    return m
end

function DualizeModel(minmax_model_pm::AbstractPowerModel)
    io = open("minmax.txt","w")
    println(io,minmax_model_pm.model)
    close(io)
    maxmin_model = dualize(minmax_model_pm.model)
    for l in ids(minmax_model_pm, :branch)
        if !(l in minmax_model_pm.data["protected_branches"] || l in minmax_model_pm.data["inactive_branches"])
     	  if has_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
     	      delete_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
     	  end
     	  if has_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
     	      delete_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
     	  end
	  JuMP.set_integer(variable_by_name(maxmin_model,"x[$l]_1"))
	end
    end
    fn_base=string(minmax_model_pm.data["name"],".cbf")
println(fn_base)
    mathoptformat_model = MathOptInterface.FileFormats.CBF.Model()
    MOI.copy_to(MOI.Bridges.full_bridge_optimizer(mathoptformat_model,Float64), backend(maxmin_model))
    #MOI.copy_to(mathoptformat_model, backend(maxmin_model))
    MOI.write_to_file(mathoptformat_model, fn_base)

    return maxmin_model
end

function MinimaxOPFModel(pm_data, powerform)
    if powerform == SOCWRConicPowerModel || powerform == SDPWRMPowerModel || powerform == SparseSDPWRMPowerModel 
        println("Prototyping Algo with WRConic Forms")
        pm = build_model(pm_data, powerform, WRConicPost_PF_Minmax)
    else
	println("Do nothing")
    end

    return pm
end

function PF_FeasModel(pm_data, powerform)
    if powerform == SOCWRConicPowerModel || powerform == SDPWRMPowerModel || powerform == SparseSDPWRMPowerModel 
        println("Prototyping Algo with WRConic Forms")
        pm = build_model(pm_data, powerform, WRConicPost_PF)
        objective_feasibility_problem(pm)
    else
	println("Do nothing")
    end
    return pm
end

function WRConicPost_PF_Minmax(pm::AbstractPowerModel)
    WRConicPost_PF(pm)
    variable_ordering_auxiliary(pm)
    con(pm, pm.cnw, pm.ccnd)[:x] = Dict{Int,ConstraintRef}()
    undecided_branches = filter(l->!(l in pm.data["protected_branches"] || l in pm.data["inactive_branches"]), ids(pm,pm.cnw,:branch))
    for l in undecided_branches
        constraint_abs_branch_flow_ordering(pm, l)
    end
    objective_minmax_problem(pm)
end

function WRConicPost_PF(pm::AbstractPowerModel)
    variable_voltage(pm)
    variable_generation(pm)
    variable_branch_flow(pm)
    variable_dcline_flow(pm)
    remove_infinity_bnds(pm)

    # Add new variables
    #variable_bus_slacks(pm)
    variable_branch_flow_slacks(pm)

    # Post objective function for testing fisibility problem
    # objective_feasibility_problem(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        #constraint_power_balance_slacks(pm, i)
        constraint_power_balance(pm, i)
    end

    for l in ids(pm, :branch)
	#if l in pm.data["protected_branches"]
        #  constraint_ohms_yt_from(pm, l)
        #  constraint_ohms_yt_to(pm, l)
	#else
        constraint_ohms_yt_from_slacks(pm, l)
        constraint_ohms_yt_to_slacks(pm, l)
        constraint_def_abs_flow_values(pm, l)
	#end

        constraint_voltage_angle_difference(pm, l)

        constraint_thermal_limit_from(pm, l)
        constraint_thermal_limit_to(pm, l)
    end

    for l in ids(pm, :dcline)
        constraint_dcline(pm, l)  # DO WE NEED TO TREAT THESE CONSTRAINTS DIFFERENTLY
    end
end



end # module
