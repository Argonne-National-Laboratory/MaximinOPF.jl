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
        [l in ids(pm, nw, :branch)], base_name="$(nw)_$(cnd)_upf1",
        lower_bound = 0,
        start = 0
    )
    upt1 = var(pm, nw, cnd)[:upt1] = JuMP.@variable(pm.model,
        [l in ids(pm, nw, :branch)], base_name="$(nw)_$(cnd)_upt1",
        lower_bound = 0,
        start = 0
    )
    uqf1 = var(pm, nw, cnd)[:uqf1] = JuMP.@variable(pm.model,
        [l in ids(pm, nw, :branch)], base_name="$(nw)_$(cnd)_uqf1",
        lower_bound = 0,
        start = 0
    )
    uqt1 = var(pm, nw, cnd)[:uqt1] = JuMP.@variable(pm.model,
        [l in ids(pm, nw, :branch)], base_name="$(nw)_$(cnd)_uqt1",
        lower_bound = 0,
        start = 0
    )
end

function variable_branch_flow_slacks0(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    upf0 = var(pm, nw, cnd)[:upf0] = JuMP.@variable(pm.model,
        [l in ids(pm, nw, :branch)], base_name="$(nw)_$(cnd)_upf0",
        lower_bound = 0,
        start = 0
    )
    upt0 = var(pm, nw, cnd)[:upt0] = JuMP.@variable(pm.model,
        [l in ids(pm, nw, :branch)], base_name="$(nw)_$(cnd)_upt0",
        lower_bound = 0,
        start = 0
    )
    uqf0 = var(pm, nw, cnd)[:uqf0] = JuMP.@variable(pm.model,
        [l in ids(pm, nw, :branch)], base_name="$(nw)_$(cnd)_uqf0",
        lower_bound = 0,
        start = 0
    )
    uqt0 = var(pm, nw, cnd)[:uqt0] = JuMP.@variable(pm.model,
        [l in ids(pm, nw, :branch)], base_name="$(nw)_$(cnd)_uqt0",
        lower_bound = 0,
        start = 0
    )
end

function variable_ordering_auxiliary(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    u_ord_aux = var(pm, nw, cnd)[:u_ord_aux] = JuMP.@variable(pm.model,
        [b in ids(pm, nw, :branch)], base_name="$(nw)_$(cnd)_u_ord_aux",
        lower_bound = 0,
        start = 0
    )
    u_K = var(pm,nw,cnd)[:u_K] = JuMP.@variable(pm.model, base_name="$(nw)_$(cnd)_u_K",lower_bound=0,start=0)
end

# Remove bound variables that are not used in MaxminOPF from PowerModel
function remove_infinity_bnds(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    p = var(pm, nw, cnd, :p)
    q = var(pm, nw, cnd, :q)
    for l in ref(pm,nw, :arcs)
      if lower_bound(p[l])==-Inf
        delete_lower_bound(p[l])
      end      
      if lower_bound(q[l])==-Inf
        delete_lower_bound(q[l])
      end      
      if upper_bound(p[l])==Inf
        delete_upper_bound(p[l])
      end      
      if upper_bound(q[l])==Inf
        delete_upper_bound(q[l])
      end      
    end
end
