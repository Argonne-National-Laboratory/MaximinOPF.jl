using JuMP

function objective_feasibility_problem(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    upbus = var(pm,nw,cnd,:up_bus)
    uqbus = var(pm,nw,cnd,:uq_bus)
    up_br = var(pm,nw,cnd,:up_br)
    uq_br = var(pm,nw,cnd,:uq_br)
    return JuMP.@objective(pm.model, Min,
        sum( upbus[i] + uqbus[i] for i in ids(pm,nw,:bus)) 
	+ sum( up_br[l,1] + uq_br[l,1] for l in ref(pm,nw,:arcs))
    )
end

function objective_minmax_problem(pm::AbstractPowerModel, K::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    upbus = var(pm,nw,cnd,:up_bus)
    uqbus = var(pm,nw,cnd,:uq_bus)
    up_br = var(pm,nw,cnd,:up_br)
    uq_br = var(pm,nw,cnd,:uq_br)
    u_ord_aux = var(pm,nw,cnd,:u_ord_aux)
    u_K = var(pm,nw,cnd,:u_K)
    return JuMP.@objective(pm.model, Min, K*u_K + sum( u_ord_aux[l] for l in ids(pm,nw,:branch) )
        + sum( upbus[i] + uqbus[i] for i in ids(pm,nw,:bus)) 
	+ sum( up_br[l,1] + uq_br[l,1] for l in ref(pm,nw,:arcs))
    )
end
