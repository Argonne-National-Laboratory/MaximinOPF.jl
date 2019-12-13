using JuMP

function variable_bus_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    up_bus = var(pm, nw, cnd)[:up_bus] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_up_bus",
	lower_bound = 0,
        start = 0
    )
    uq_bus = var(pm, nw, cnd)[:uq_bus] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_uq_bus",
	lower_bound = 0,
        start = 0
    )
end

function variable_branch_flow_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    up_br = var(pm, nw, cnd)[:up_br] = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs),k=0:1], base_name="$(nw)_$(cnd)_up_br",
            lower_bound = 0,
            start = 0
    )
    uq_br = var(pm, nw, cnd)[:uq_br] = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs),k=0:1], base_name="$(nw)_$(cnd)_uq_br",
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
     if has_lower_bound(p[l])
      if lower_bound(p[l])==-Inf 
        delete_lower_bound(p[l])
      end      
     end
     if has_lower_bound(q[l])
      if lower_bound(q[l])==-Inf 
        delete_lower_bound(q[l])
      end      
     end
     if has_upper_bound(p[l])
      if upper_bound(p[l])==Inf 
        delete_upper_bound(p[l])
      end      
     end
     if has_upper_bound(q[l])
      if upper_bound(q[l])==Inf 
        delete_upper_bound(q[l])
      end      
     end
    end
end
