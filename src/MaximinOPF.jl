module MaximinOPF
using PowerModels
include("Variables.jl")
include("Objectives.jl")
include("Constraints.jl")
greet() = print("Hello World!")

function MaximinOPFModel(case, powerform, nLineAttacked)
	println("Hello MaximinOPFModel")
	println(case)
	println(powerform)
	println(nLineAttacked)
	pm = ""
	#Output Model from PowerModels

    if powerform == SOCWRConicPowerModel
        println("Prototyping Algo")
        pm = build_model(case, powerform, SOCWRConicPost_OPF)
	elseif (powerform == SDPWRMPowerModel)
		println("Brian Algo")
		pm = build_model(case, powerform, post_pf_feas)
	else
		println("Default Algo")
		pm = build_model(case, powerform, myPost_OPF)
	end

	return pm
end

function SOCWRConicPost_OPF(pm::AbstractPowerModel)
    variable_voltage(pm)
    variable_generation(pm)
    variable_branch_flow(pm)
    variable_dcline_flow(pm)

    # Add new variables
    variable_bus_slacks(pm)
    variable_branch_flow_slacks(pm)

    # Replace objective function
    # objective_min_fuel_and_flow_cost(pm)
    objective_feasibility_problem(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        # Replace power balnce constraint function
        # constraint_power_balance(pm, i)
        constraint_power_balance_slacks(pm, i)
    end

    for i in ids(pm, :branch)
        # Replace line flow constraint function
        # constraint_ohms_yt_from(pm, i)
        # constraint_ohms_yt_to(pm, i)
        constraint_ohms_yt_from_slacks(pm, i)
        constraint_ohms_yt_to_slacks(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline(pm, i)
    end
end

function myPost_OPF(pm::AbstractPowerModel)
	variable_voltage(pm)
    variable_generation(pm)
    variable_branch_flow(pm)
    variable_dcline_flow(pm)

    objective_min_fuel_and_flow_cost(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline(pm, i)
    end
end

function post_pf_feas(pm::AbstractPowerModel) # Only Working for MoSek Solver
    variable_voltage(pm, bounded = false)
    variable_generation(pm, bounded = false)
    variable_branch_flow(pm, bounded = false)
    variable_dcline_flow(pm, bounded = false)

    constraint_model_voltage(pm)
    
    # ADD SLACK VARS
    #@variable(pm.model,
    #@objective(pm.model, )

    # ADD OBJECTIVE FUNCTION THAT SUMS ALL SLACKS

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3
        constraint_theta_ref(pm, i)
        constraint_voltage_magnitude_setpoint(pm, i)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_power_balance(pm, i) # REPLACE

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_voltage_magnitude_setpoint(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_active_gen_setpoint(pm, j)
            end
        end
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i) # REPLACE
        constraint_ohms_yt_to(pm, i) # REPLACE
    end

    for (i,dcline) in ref(pm, :dcline)
        #constraint_dcline(pm, i) not needed, active power flow fully defined by dc line setpoints
        constraint_active_dcline_setpoint(pm, i)

        f_bus = ref(pm, :bus)[dcline["f_bus"]]
        if f_bus["bus_type"] == 1
            constraint_voltage_magnitude_setpoint(pm, f_bus["index"])
        end

        t_bus = ref(pm, :bus)[dcline["t_bus"]]
        if t_bus["bus_type"] == 1
            constraint_voltage_magnitude_setpoint(pm, t_bus["index"])
        end
    end
end

end # module