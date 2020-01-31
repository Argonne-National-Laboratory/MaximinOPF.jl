using JuMP

function objective_feasibility_problem(pm::AbstractPowerModel, x_vals::Dict{Int64,Float64}=Dict{Int64,Float64}(); nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    up_br1 = var(pm,nw,cnd,:up_br1)
    up_br0 = var(pm,nw,cnd,:up_br0)
    uq_br1 = var(pm,nw,cnd,:uq_br1)
    uq_br0 = var(pm,nw,cnd,:uq_br0)
    protected_arcs = filter(a->(a[1] in pm.data["protected_branches"]),ref(pm,nw,:arcs))
    inactive_arcs = filter(a->(a[1] in pm.data["inactive_branches"]),ref(pm,nw,:arcs))
    undecided_arcs = filter(a->!(a in protected_arcs || a in inactive_arcs),ref(pm,nw,:arcs))
    return JuMP.@objective(pm.model, Min, 
	sum( (1-x_vals[a[1]])*(up_br1[a,0] + uq_br1[a,0] + up_br1[a,1] + uq_br1[a,1]) 
            + x_vals[a[1]]*(up_br0[a,0] + uq_br0[a,0] + up_br0[a,1] + uq_br0[a,1]) for a in undecided_arcs)
            + sum( (up_br1[a,0] + uq_br1[a,0] + up_br1[a,1] + uq_br1[a,1])  for a in protected_arcs)
            + sum( (up_br0[a,0] + uq_br0[a,0] + up_br0[a,1] + uq_br0[a,1])  for a in inactive_arcs)
    )
end

function objective_minmax_problem(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    K = pm.data["attacker_budget"]
    up_br1 = var(pm,nw,cnd,:up_br1)
    up_br0 = var(pm,nw,cnd,:up_br0)
    uq_br1 = var(pm,nw,cnd,:uq_br1)
    uq_br0 = var(pm,nw,cnd,:uq_br0)
    u_ord_aux = var(pm,nw,cnd,:u_ord_aux)
    u_K = var(pm,nw,cnd,:u_K)
    undecided_branches = filter(l->!(l in pm.data["protected_branches"] || l in pm.data["inactive_branches"]), ids(pm,nw,:branch))
    undecided_arcs = filter(a->(a[1] in undecided_branches),ref(pm,nw,:arcs))
    protected_arcs = filter(a->(a[1] in pm.data["protected_branches"]),ref(pm,nw,:arcs))
    attacked_arcs = filter(a->(a[1] in pm.data["inactive_branches"]),ref(pm,nw,:arcs))
    return JuMP.@objective(pm.model, Min, K*u_K + sum( u_ord_aux[l] for l in undecided_branches ) 
            + sum( up_br1[a,0] + uq_br1[a,0] + up_br1[a,1] + uq_br1[a,1] for a in undecided_arcs)
            + sum( up_br1[a,0] + uq_br1[a,0] + up_br1[a,1] + uq_br1[a,1] for a in protected_arcs)
            + sum( up_br0[a,0] + uq_br0[a,0] + up_br0[a,1] + uq_br0[a,1] for a in attacked_arcs)
    )
end
