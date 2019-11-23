using JuMP

function variable_bus_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
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

function variable_branch_flow_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
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
