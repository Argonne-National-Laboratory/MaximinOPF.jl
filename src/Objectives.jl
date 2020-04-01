using JuMP

function objective_feasibility_problem(pm::AbstractPowerModel, x_vals::Dict{Int64,Float64}=Dict{Int64,Float64}(); nw::Int=pm.cnw)
    for l in ids(pm,pm.cnw,:branch)
        if !haskey(x_vals,l)
	        if l in pm.data["inactive_branches"]
		        x_vals[l]=1
	        else
		        x_vals[l]=0
	        end
        end
    end
    pd_br = var(pm,nw,:pd_br)
    pt_br = var(pm,nw,:pt_br)
    qd_br = var(pm,nw,:qd_br)
    qt_br = var(pm,nw,:qt_br)
    protected_arcs = filter(a->(a[1] in pm.data["protected_branches"]),ref(pm,nw,:arcs))
    inactive_arcs = filter(a->(a[1] in pm.data["inactive_branches"]),ref(pm,nw,:arcs))
    undecided_arcs = filter(a->!(a in protected_arcs || a in inactive_arcs),ref(pm,nw,:arcs))


    return JuMP.@objective(pm.model, Min, 
	sum( (1-x_vals[a[1]])*(pd_br[a,0] + qd_br[a,0] + pd_br[a,1] + qd_br[a,1]) 
            + x_vals[a[1]]*(pt_br[a,0] + qt_br[a,0] + pt_br[a,1] + qt_br[a,1]) for a in undecided_arcs)
            + sum( (pd_br[a,0] + qd_br[a,0] + pd_br[a,1] + qd_br[a,1])  for a in protected_arcs)
            + sum( (pt_br[a,0] + qt_br[a,0] + pt_br[a,1] + qt_br[a,1])  for a in inactive_arcs)
    )
end

function objective_minmax_problem(pm::AbstractPowerModel; nw::Int=pm.cnw)
    K = pm.data["attacker_budget"]
    pd_br = var(pm,nw,:pd_br)
    pt_br = var(pm,nw,:pt_br)
    qd_br = var(pm,nw,:qd_br)
    qt_br = var(pm,nw,:qt_br)
    u_ord_aux = var(pm,nw,:u_ord_aux)
    u_K = var(pm,nw,:u_K)
    undecided_branches = filter(l->!(l in pm.data["protected_branches"] || l in pm.data["inactive_branches"]), ids(pm,nw,:branch))
    undecided_arcs = filter(a->(a[1] in undecided_branches),ref(pm,nw,:arcs))
    protected_arcs = filter(a->(a[1] in pm.data["protected_branches"]),ref(pm,nw,:arcs))
    attacked_arcs = filter(a->(a[1] in pm.data["inactive_branches"]),ref(pm,nw,:arcs))

    return JuMP.@objective(pm.model, Min, K*u_K + sum( u_ord_aux[l] for l in undecided_branches ) 
            + sum( pd_br[a,0] + qd_br[a,0] + pd_br[a,1] + qd_br[a,1] for a in undecided_arcs)
            + sum( pd_br[a,0] + qd_br[a,0] + pd_br[a,1] + qd_br[a,1] for a in protected_arcs)
            + sum( pt_br[a,0] + qt_br[a,0] + pt_br[a,1] + qt_br[a,1] for a in attacked_arcs)
    )
end
