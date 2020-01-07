using JuMP

function constraint_power_balance_slacks(pm::AbstractWModels, n::Int, c::Int, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    w    = var(pm, n, c, :w, i)
    p    = get(var(pm, n, c),    :p, Dict()); _check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(var(pm, n, c),    :q, Dict()); _check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(var(pm, n, c),   :pg, Dict()); _check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(var(pm, n, c),   :qg, Dict()); _check_var_keys(qg, bus_gens, "reactive power", "generator")
    ps   = get(var(pm, n, c),   :ps, Dict()); _check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(var(pm, n, c),   :qs, Dict()); _check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(var(pm, n, c),  :psw, Dict()); _check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(var(pm, n, c),  :qsw, Dict()); _check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(var(pm, n, c), :p_dc, Dict()); _check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(var(pm, n, c), :q_dc, Dict()); _check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")
    upbus = var(pm, n, c, :up_bus, i)
    uqbus = var(pm, n, c, :uq_bus, i)


    ### Active power balance -up <= ... <= up
    con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        - sum(pg[g] for g in bus_gens)
        + sum(ps[s] for s in bus_storage)
        + sum(pd for pd in values(bus_pd))
        + sum(gs for gs in values(bus_gs))*w
	- upbus <= 0
    )
    con(pm, n, c, :kcl_p)[i] = JuMP.@constraint(pm.model,
        -sum(p[a] for a in bus_arcs)
        - sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        - sum(psw[a_sw] for a_sw in bus_arcs_sw)
        + sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd for pd in values(bus_pd))
        - sum(gs for gs in values(bus_gs))*w
	- upbus <= 0 
    )

    ### Reactive power balance -uq <= ... <= uq
    con(pm, n, c, :kcl_q)[i] = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        - sum(qg[g] for g in bus_gens)
        + sum(qs[s] for s in bus_storage)
        + sum(qd for qd in values(bus_qd))
        - sum(bs for bs in values(bus_bs))*w
	- uqbus <= 0
    )
    con(pm, n, c, :kcl_q)[i] = JuMP.@constraint(pm.model,
        -sum(q[a] for a in bus_arcs)
        - sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        - sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        + sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd for qd in values(bus_qd))
        + sum(bs for bs in values(bus_bs))*w
	- uqbus <= 0
    )
end

###NOTE: This is modeled on what is in PowerModels.jl/blob/master/src/core/constraint_template.jl
function constraint_power_balance_slacks(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(con(pm, nw, cnd), :kcl_p)
        con(pm, nw, cnd)[:kcl_p] = Dict{Int,JuMP.ConstraintRef}()
    end
    if !haskey(con(pm, nw, cnd), :kcl_q)
        con(pm, nw, cnd)[:kcl_q] = Dict{Int,JuMP.ConstraintRef}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)
    bus_storage = ref(pm, nw, :bus_storage, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_power_balance_slacks(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_yt_from_slacks(pm::AbstractWRModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_fr, b_fr, tr, ti, tm)
    l=f_idx[1]
    p_fr = var(pm, n, c, :p, f_idx)
    q_fr = var(pm, n, c, :q, f_idx)
    w_fr = var(pm, n, c, :w, f_bus)
    wr   = var(pm, n, c, :wr, (f_bus, t_bus))
    wi   = var(pm, n, c, :wi, (f_bus, t_bus))
    upf1 = var(pm,n,c, :up_br)[f_idx,1]
    uqf1 = var(pm,n,c, :uq_br)[f_idx,1]

    cref=JuMP.@constraint(pm.model,  (g+g_fr)/tm^2*w_fr + (-g*tr+b*ti)/tm^2*wr + (-b*tr-g*ti)/tm^2*wi - p_fr - upf1 <= 0)
      JuMP.set_name(cref,"lambda_f+[$l]")  
    cref=JuMP.@constraint(pm.model, -(g+g_fr)/tm^2*w_fr - (-g*tr+b*ti)/tm^2*wr - (-b*tr-g*ti)/tm^2*wi + p_fr - upf1 <= 0)
      JuMP.set_name(cref,"lambda_f-[$l]")  
    cref=JuMP.@constraint(pm.model, -(b+b_fr)/tm^2*w_fr - (-b*tr-g*ti)/tm^2*wr + (-g*tr+b*ti)/tm^2*wi - q_fr - uqf1 <= 0)
      JuMP.set_name(cref,"mu_f+[$l]")  
    cref=JuMP.@constraint(pm.model, (b+b_fr)/tm^2*w_fr + (-b*tr-g*ti)/tm^2*wr - (-g*tr+b*ti)/tm^2*wi + q_fr - uqf1 <= 0)
      JuMP.set_name(cref,"mu_f-[$l]")  
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_yt_to_slacks(pm::AbstractWRModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    l=t_idx[1]
    q_to = var(pm, n, c, :q, t_idx)
    p_to = var(pm, n, c, :p, t_idx)
    w_to = var(pm, n, c, :w, t_bus)
    wr   = var(pm, n, c, :wr, (f_bus, t_bus))
    wi   = var(pm, n, c, :wi, (f_bus, t_bus))
    upt1 = var(pm,n,c, :up_br)[t_idx,1]
    uqt1 = var(pm,n,c, :uq_br)[t_idx,1]

    cref=JuMP.@constraint(pm.model,  (g+g_to)*w_to + (-g*tr-b*ti)/tm^2*wr + (-b*tr+g*ti)/tm^2*-wi - p_to - upt1 <= 0)
      JuMP.set_name(cref,"lambda_t+[$l]")  
    cref=JuMP.@constraint(pm.model,  -(g+g_to)*w_to - (-g*tr-b*ti)/tm^2*wr - (-b*tr+g*ti)/tm^2*-wi + p_to - upt1 <= 0)
      JuMP.set_name(cref,"lambda_t-[$l]")  
    cref=JuMP.@constraint(pm.model, -(b+b_to)*w_to - (-b*tr+g*ti)/tm^2*wr + (-g*tr-b*ti)/tm^2*-wi - q_to - uqt1 <= 0)
      JuMP.set_name(cref,"mu_t+[$l]")  
    cref=JuMP.@constraint(pm.model, (b+b_to)*w_to + (-b*tr+g*ti)/tm^2*wr - (-g*tr-b*ti)/tm^2*-wi + q_to - uqt1 <= 0)
      JuMP.set_name(cref,"mu_t-[$l]")  
end

### Branch - Ohm's Law Constraints ###

""
function constraint_ohms_yt_from_slacks(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_fr = branch["g_fr"][cnd]
    b_fr = branch["b_fr"][cnd]
    tm = branch["tap"][cnd]

    constraint_ohms_yt_from_slacks(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g[cnd,cnd], b[cnd,cnd], g_fr, b_fr, tr[cnd], ti[cnd], tm)
end


""
function constraint_ohms_yt_to_slacks(pm::AbstractPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)

    g, b = calc_branch_y(branch)
    tr, ti = calc_branch_t(branch)
    g_to = branch["g_to"][cnd]
    b_to = branch["b_to"][cnd]
    tm = branch["tap"][cnd]

    constraint_ohms_yt_to_slacks(pm, nw, cnd, f_bus, t_bus, f_idx, t_idx, g[cnd,cnd], b[cnd,cnd], g_to, b_to, tr[cnd], ti[cnd], tm)
end

function constraint_def_abs_flow_values(pm::AbstractPowerModel, l::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, l)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (l, f_bus, t_bus)
    t_idx = (l, t_bus, f_bus)

    p_fr = var(pm, nw, cnd, :p, f_idx)
    p_to = var(pm, nw, cnd, :p, t_idx)
    q_fr = var(pm, nw, cnd, :q, f_idx)
    q_to = var(pm, nw, cnd, :q, t_idx)
    upf0 = var(pm,nw,cnd, :up_br)[f_idx,0]
    upt0 = var(pm,nw,cnd, :up_br)[t_idx,0]
    uqf0 = var(pm,nw,cnd, :uq_br)[f_idx,0]
    uqt0 = var(pm,nw,cnd, :uq_br)[t_idx,0]

    cref=JuMP.@constraint(pm.model, p_fr - upf0 <= 0)
      JuMP.set_name(cref,"pi_f+[$l]")  
    cref=JuMP.@constraint(pm.model, -p_fr - upf0 <= 0)
      JuMP.set_name(cref,"pi_f-[$l]")  
    cref=JuMP.@constraint(pm.model, p_to - upt0 <= 0)
      JuMP.set_name(cref,"pi_t+[$l]")  
    cref=JuMP.@constraint(pm.model, -p_to - upt0 <= 0)
      JuMP.set_name(cref,"pi_t-[$l]")  

    cref=JuMP.@constraint(pm.model, q_fr - uqf0 <= 0)
      JuMP.set_name(cref,"phi_f+[$l]")  
    cref=JuMP.@constraint(pm.model, -q_fr - uqf0 <= 0)
      JuMP.set_name(cref,"phi_f-[$l]")  
    cref=JuMP.@constraint(pm.model, q_to - uqt0 <= 0)
      JuMP.set_name(cref,"phi_t+[$l]")  
    cref=JuMP.@constraint(pm.model, -q_to - uqt0 <= 0)
      JuMP.set_name(cref,"phi_t-[$l]")  
end

function constraint_abs_branch_flow_ordering(pm::AbstractPowerModel, l::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, l)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (l, f_bus, t_bus)
    t_idx = (l, t_bus, f_bus)

    upf0 = var(pm,nw,cnd, :up_br)[f_idx,0]
    upt0 = var(pm,nw,cnd, :up_br)[t_idx,0]
    uqf0 = var(pm,nw,cnd, :uq_br)[f_idx,0]
    uqt0 = var(pm,nw,cnd, :uq_br)[t_idx,0]
    upf1 = var(pm,nw,cnd, :up_br)[f_idx,1]
    upt1 = var(pm,nw,cnd, :up_br)[t_idx,1]
    uqf1 = var(pm,nw,cnd, :uq_br)[f_idx,1]
    uqt1 = var(pm,nw,cnd, :uq_br)[t_idx,1]
    u_ord_aux = var(pm,nw,cnd,:u_ord_aux,l)
    u_K = var(pm,nw,cnd,:u_K)

    con(pm, nw, cnd)[:x][l] = 
    JuMP.@constraint(pm.model, -(upf0 + upt0 + uqf0 + uqt0) + (upf1 + upt1 + uqf1 + uqt1) + u_ord_aux + u_K >= 0)

    JuMP.set_name(con(pm, nw, cnd)[:x][l],"x[$l]")  
end

"checks if a sufficient number of variables exist for the given keys collection"
function _check_var_keys(vars, keys, var_name, comp_name)
    if length(vars) < length(keys)
        error(_LOGGER, "$(var_name) decision variables appear to be missing for $(comp_name) components")
    end
end

