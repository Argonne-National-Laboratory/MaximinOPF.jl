using JuMP

function objective_feasibility_problem(pm::AbstractPowerModel,nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    upbus = var(pm,nw,cnd,:upbus)
    uqbus = var(pm,nw,cnd,:uqbus)
    upf1 = var(pm,nw,cnd,:upf1)
    upt1 = var(pm,nw,cnd,:upt1)
    uqf1 = var(pm,nw,cnd,:uqf1)
    uqt1 = var(pm,nw,cnd,:uqt1)
    return JuMP.@objective(pm.model, Min,
        sum( upbus[i] + uqbus[i] for i in ids(pm,nw,:bus)) 
	+ sum( upf1[l] + uqf1[l] for l in ref(pm,nw,:arcs_from))
	+ sum( upt1[l] + uqt1[l] for l in ref(pm,nw,:arcs_to))
    )
end

        [b in ref(pm, nw, :branch)], base_name="$(nw)_$(cnd)_u_ord_aux",
function objective_feasibility_problem(pm::AbstractPowerModel,nw::Int=pm.cnw, cnd::Int=pm.ccnd, K::Int)
    upbus = var(pm,nw,cnd,:upbus)
    uqbus = var(pm,nw,cnd,:uqbus)
    upf1 = var(pm,nw,cnd,:upf1)
    upt1 = var(pm,nw,cnd,:upt1)
    uqf1 = var(pm,nw,cnd,:uqf1)
    uqt1 = var(pm,nw,cnd,:uqt1)
    u_ord_aux = var(pm,nw,cnd,:u_ord_aux)
    u_K = var(pm,nw,cnd,:u_K)
    return JuMP.@objective(pm.model, Min, K*u_K + sum( u_ord_aux[l] for l in ids(pm,nw,:branch) )
        + sum( upbus[i] + uqbus[i] for i in ids(pm,nw,:bus)) 
	+ sum( upf1[l] + uqf1[l] for l in ref(pm,nw,:arcs_from))
	+ sum( upt1[l] + uqt1[l] for l in ref(pm,nw,:arcs_to))
    )
end
