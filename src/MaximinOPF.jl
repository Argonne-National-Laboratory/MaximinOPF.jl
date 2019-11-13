module MaximinOPF
using PowerModels

greet() = print("Hello World!")

function MaximinOPFModel(case, powerform, nLineAttacked)
	println("Hello MaximinOPFModel")
	println(case)
	println(powerform)
	println(nLineAttacked)

	#Output Model from PowerModels
	pm = build_model(case, powerform, myPost_OPF)

	return pm
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

end # module