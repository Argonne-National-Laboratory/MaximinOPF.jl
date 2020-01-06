using JuMP

function objective_feasibility_problem(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    up_br = var(pm,nw,cnd,:up_br)
    uq_br = var(pm,nw,cnd,:uq_br)
    active_arcs = filter(a->!(a[1] in pm.data["inactive_branches"]),ref(pm,nw,:arcs))
    inactive_arcs = filter(a->(a[1] in pm.data["inactive_branches"]),ref(pm,nw,:arcs))
    return JuMP.@objective(pm.model, Min,
	+ sum( up_br[a,1] + uq_br[a,1] for a in active_arcs)
	+ sum( up_br[a,0] + uq_br[a,0] for a in inactive_arcs)
    )
end

function objective_minmax_problem(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    K = pm.data["attacker_budget"]
    #upbus = var(pm,nw,cnd,:up_bus)
    #uqbus = var(pm,nw,cnd,:uq_bus)
    up_br = var(pm,nw,cnd,:up_br)
    uq_br = var(pm,nw,cnd,:uq_br)
    u_ord_aux = var(pm,nw,cnd,:u_ord_aux)
    u_K = var(pm,nw,cnd,:u_K)
    undecided_branches = filter(l->!(l in pm.data["protected_branches"] || l in pm.data["inactive_branches"]), ids(pm,nw,:branch))
    undecided_arcs = filter(a->(a[1] in undecided_branches),ref(pm,nw,:arcs))
    protected_arcs = filter(a->(a[1] in pm.data["protected_branches"]),ref(pm,nw,:arcs))
    attacked_arcs = filter(a->(a[1] in pm.data["inactive_branches"]),ref(pm,nw,:arcs))
    return JuMP.@objective(pm.model, Min, K*u_K + sum( u_ord_aux[l] for l in undecided_branches ) + sum( up_br[a,1] + uq_br[a,1] for a in undecided_arcs)
	+ sum( up_br[a,1] + uq_br[a,1] for a in protected_arcs)
	+ sum( up_br[a,0] + uq_br[a,0] for a in attacked_arcs)
    )
end
