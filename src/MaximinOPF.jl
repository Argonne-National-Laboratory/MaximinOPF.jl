module MaximinOPF
using PowerModels
using JuMP

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

function variable_bus_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded::Bool = true)
    upbus = var(pm, nw, cnd)[:upbus] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_upbus",
	lower_bound = 0,
        start = 0
    )
    uqbus = var(pm, nw, cnd)[:uqbus] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_uqbus",
	lower_bound = 0,
        start = 0
    )
end

function variable_branch_flow_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    upf1 = var(pm, nw, cnd)[:upf1] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_upf1",
        lower_bound = 0,
        start = 0
    )
    upt1 = var(pm, nw, cnd)[:upt1] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_upt1",
        lower_bound = 0,
        start = 0
    )
    uqf1 = var(pm, nw, cnd)[:uqf1] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_uqf1",
        lower_bound = 0,
        start = 0
    )
    uqt1 = var(pm, nw, cnd)[:uqt1] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_uqt1",
        lower_bound = 0,
        start = 0
    )
end

function objective_feasibility_problem(pm::AbstractPowerModel,nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    upbus = var(pm,nw,cnd,:upbus)
    uqbus = var(pm,nw,cnd,:uqbus)
    upf1 = var(pm,nw,cnd,:upf1)
    upt1 = var(pm,nw,cnd,:upt1)
    uqf1 = var(pm,nw,cnd,:uqf1)
    uqt1 = var(pm,nw,cnd,:uqt1)
    return JuMP.@objective(pm.model, Min,
        sum( upbus[i] + uqbus[i] for i in ids(pm,nw,:bus)) + sum( upf1[l] + upt1[l] + uqf1[l] + uqt1[l] for l in ref(pm,nw,:arcs))
    )
end


function post_pf_feas(pm::AbstractPowerModel)

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
