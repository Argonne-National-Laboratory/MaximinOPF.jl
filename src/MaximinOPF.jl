module MaximinOPF
using PowerModels
using JuMP

greet() = print("Hello World!")

function MaximinOPFModel(case, powerform, nLineAttacked)
	println("Hello MaximinOPFModel")
	println(case)
	println(powerform)
	println(nLineAttacked)

	#Output Model from PowerModels
	pm = build_model(case, powerform, myPost_OPF)

	return pm
end

function variable_bus_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded::Bool = true)
    upbus = var(pm, nw, cnd)[:upbus] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_upbus",
	lower_bound = 0,
        start = 0
    )
    uqbus = var(pm, nw, cnd)[:uqbus] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :bus)], base_name="$(nw)_$(cnd)_uqbus",
	lower_bound = 0,
        start = 0
    )
end

function variable_branch_flow_slacks(pm::AbstractPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    upf1 = var(pm, nw, cnd)[:upf1] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_upf1",
        lower_bound = 0,
        start = 0
    )
    upt1 = var(pm, nw, cnd)[:upt1] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_upt1",
        lower_bound = 0,
        start = 0
    )
    uqf1 = var(pm, nw, cnd)[:uqf1] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_uqf1",
        lower_bound = 0,
        start = 0
    )
    uqt1 = var(pm, nw, cnd)[:uqt1] = JuMP.@variable(pm.model,
        [(l,i,j) in ref(pm, nw, :arcs)], base_name="$(nw)_$(cnd)_uqt1",
        lower_bound = 0,
        start = 0
    )
end

function objective_feasibility_problem(pm::AbstractPowerModel,nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    upbus = var(pm,nw,cnd,:upbus)
    uqbus = var(pm,nw,cnd,:uqbus)
    upf1 = var(pm,nw,cnd,:upf1)
    upt1 = var(pm,nw,cnd,:upt1)
    uqf1 = var(pm,nw,cnd,:uqf1)
    uqt1 = var(pm,nw,cnd,:uqt1)
    return JuMP.@objective(pm.model, Min,
        sum( upbus[i] + uqbus[i] for i in ids(pm,nw,:bus)) + sum( upf1[l] + upt1[l] + uqf1[l] + uqt1[l] for l in ref(pm,nw,:arcs))
    )
end

### This is based on the function from PowerModels/src/form/shared.jl
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
    upbus = get(var(pm,n,c,i), :upbus, Dict())
    uqbus = get(var(pm,n,c,i), :uqbus, Dict())


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
    p_fr = var(pm, n, c, :p, f_idx)
    q_fr = var(pm, n, c, :q, f_idx)
    w_fr = var(pm, n, c, :w, f_bus)
    wr   = var(pm, n, c, :wr, (f_bus, t_bus))
    wi   = var(pm, n, c, :wi, (f_bus, t_bus))
    upf1 = var(pm,n,c, :upf1, f_idx)
    uqf1 = var(pm,n,c, :uqf1, f_idx)

    JuMP.@constraint(pm.model,  (g+g_fr)/tm^2*w_fr + (-g*tr+b*ti)/tm^2*wr + (-b*tr-g*ti)/tm^2*wi - p_fr - upf1 <= 0)
    JuMP.@constraint(pm.model, -(g+g_fr)/tm^2*w_fr - (-g*tr+b*ti)/tm^2*wr - (-b*tr-g*ti)/tm^2*wi + p_fr - upf1 <= 0)
    JuMP.@constraint(pm.model, -(b+b_fr)/tm^2*w_fr - (-b*tr-g*ti)/tm^2*wr + (-g*tr+b*ti)/tm^2*wi - uqf1 <= 0)
    JuMP.@constraint(pm.model, (b+b_fr)/tm^2*w_fr + (-b*tr-g*ti)/tm^2*wr - (-g*tr+b*ti)/tm^2*wi - uqf1 <= 0)
end


"""
Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)
"""
function constraint_ohms_yt_to_slacks(pm::AbstractWRModels, n::Int, c::Int, f_bus, t_bus, f_idx, t_idx, g, b, g_to, b_to, tr, ti, tm)
    q_to = var(pm, n, c, :q, t_idx)
    p_to = var(pm, n, c, :p, t_idx)
    w_to = var(pm, n, c, :w, t_bus)
    wr   = var(pm, n, c, :wr, (f_bus, t_bus))
    wi   = var(pm, n, c, :wi, (f_bus, t_bus))
    upt1 = var(pm,n,c, :upt1, t_idx)
    uqt1 = var(pm,n,c, :uqt1, t_idx)

    JuMP.@constraint(pm.model,  (g+g_to)*w_to + (-g*tr-b*ti)/tm^2*wr + (-b*tr+g*ti)/tm^2*-wi - p_to - upt1 <= 0)
    JuMP.@constraint(pm.model,  -(g+g_to)*w_to - (-g*tr-b*ti)/tm^2*wr - (-b*tr+g*ti)/tm^2*-wi + p_to - upt1 <= 0)
    JuMP.@constraint(pm.model, -(b+b_to)*w_to - (-b*tr+g*ti)/tm^2*wr + (-g*tr-b*ti)/tm^2*-wi - q_to - uqt1 <= 0)
    JuMP.@constraint(pm.model, (b+b_to)*w_to + (-b*tr+g*ti)/tm^2*wr - (-g*tr-b*ti)/tm^2*-wi + q_to - uqt1 <= 0)
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


function myPost_OPF(pm::AbstractPowerModel)
	variable_voltage(pm)
    variable_generation(pm)
    variable_branch_flow(pm)
    variable_dcline_flow(pm)

    objective_min_fuel_and_flow_cost(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline(pm, i)
    end
end

end # module
