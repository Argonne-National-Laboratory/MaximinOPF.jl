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
