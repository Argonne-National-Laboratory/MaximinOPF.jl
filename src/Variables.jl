using JuMP

function variable_branch_flow_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw)
    ### "line flow discrepancies (d), with separated positive and negative components"
    pd_br = var(pm, nw)[:pd_br] = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm,nw,:arcs),k=0:1], base_name="$(nw)_pd_br",
            lower_bound = 0,
            start = 0
    )
    qd_br = var(pm, nw)[:qd_br] = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm,nw,:arcs),k=0:1], base_name="$(nw)_qd_br",
            lower_bound = 0,
            start = 0
    )

    ### "line flow targets (t), with separated positive and negative components"
    pt_br = var(pm, nw)[:pt_br] = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm,nw,:arcs),k=0:1], base_name="$(nw)_pt_br",
            lower_bound = 0,
            start = 0
    )
    qt_br = var(pm, nw)[:qt_br] = JuMP.@variable(pm.model,
            [(l,i,j) in ref(pm,nw,:arcs),k=0:1], base_name="$(nw)_qt_br",
            lower_bound = 0,
            start = 0
    )
end

function variable_ordering_auxiliary(pm::AbstractPowerModel; nw::Int=pm.cnw)
    undecided_branches = filter(l->!(l in pm.data["protected_branches"] || l in pm.data["inactive_branches"]), ids(pm,nw,:branch))
    u_ord_aux = var(pm, nw)[:u_ord_aux] = JuMP.@variable(pm.model,
        [b in undecided_branches], base_name="$(nw)_u_ord_aux",
        lower_bound = 0,
        start = 0
    )
    u_K = var(pm,nw)[:u_K] = JuMP.@variable(pm.model, base_name="$(nw)_u_K",lower_bound=0,start=0)
end

