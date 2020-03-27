using JuMP

function variable_branch_flow_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw)
    #u1_arcs = filter(a->!(a[1] in pm.data["inactive_branches"]),ref(pm,nw,:arcs)) #Exclude inactive arcs
    u1_arcs = ref(pm,nw,:arcs)
    up_br1 = var(pm, nw)[:up_br1] = JuMP.@variable(pm.model,
            [(l,i,j) in u1_arcs,k=0:1], base_name="$(nw)_up_br1",
            lower_bound = 0,
            start = 0
    )
    uq_br1 = var(pm, nw)[:uq_br1] = JuMP.@variable(pm.model,
            [(l,i,j) in u1_arcs,k=0:1], base_name="$(nw)_uq_br1",
            lower_bound = 0,
            start = 0
    )

    #u0_arcs = filter(a->!(a[1] in pm.data["protected_branches"]),ref(pm,nw,:arcs)) #Exclude protected arcs
    u0_arcs = ref(pm,nw,:arcs)
    up_br0 = var(pm, nw)[:up_br0] = JuMP.@variable(pm.model,
            [(l,i,j) in u0_arcs,k=0:1], base_name="$(nw)_up_br0",
            lower_bound = 0,
            start = 0
    )
    uq_br0 = var(pm, nw)[:uq_br0] = JuMP.@variable(pm.model,
            [(l,i,j) in u0_arcs,k=0:1], base_name="$(nw)_uq_br0",
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
