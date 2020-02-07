using JuMP
using LinearAlgebra


function constraint_def_abs_flow_values(pm::AbstractPowerModel, l::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, l)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (l, f_bus, t_bus)
    t_idx = (l, t_bus, f_bus)

    cref_pf,cref_pt=nothing,nothing
    if haskey( var(pm,nw), :p )
      p_fr = var(pm, nw,:p, f_idx)
      p_to = var(pm, nw,:p, t_idx)

      upf0m = var(pm,nw, :up_br0)[f_idx,0]
      upf0p = var(pm,nw, :up_br0)[f_idx,1]
      upt0m = var(pm,nw, :up_br0)[t_idx,0]
      upt0p = var(pm,nw, :up_br0)[t_idx,1]

      upf1m = var(pm,nw, :up_br1)[f_idx,0]
      upf1p = var(pm,nw, :up_br1)[f_idx,1]
      upt1m = var(pm,nw, :up_br1)[t_idx,0]
      upt1p = var(pm,nw, :up_br1)[t_idx,1]

      cref_pf=JuMP.@constraint(pm.model, p_fr - upf1m + upf1p - upf0m + upf0p == 0)
      cref_pt=JuMP.@constraint(pm.model, p_to - upt1m + upt1p - upt0m + upt0p == 0)
    end
    
    cref_qf,cref_qt=nothing,nothing
    if haskey( var(pm,nw), :q )
      q_fr = var(pm, nw, :q, f_idx)
      q_to = var(pm, nw, :q, t_idx)

      uqf0m = var(pm,nw, :uq_br0)[f_idx,0]
      uqf0p = var(pm,nw, :uq_br0)[f_idx,1]
      uqt0m = var(pm,nw, :uq_br0)[t_idx,0]
      uqt0p = var(pm,nw, :uq_br0)[t_idx,1]

      uqf1m = var(pm,nw, :uq_br1)[f_idx,0]
      uqf1p = var(pm,nw, :uq_br1)[f_idx,1]
      uqt1m = var(pm,nw, :uq_br1)[t_idx,0]
      uqt1p = var(pm,nw, :uq_br1)[t_idx,1]

      cref_qf=JuMP.@constraint(pm.model, q_fr - uqf1m + uqf1p - uqf0m + uqf0p == 0)
      cref_qt=JuMP.@constraint(pm.model, q_to - uqt1m + uqt1p - uqt0m + uqt0p == 0)
    end
    return cref_pf,cref_pt,cref_qf,cref_qt
end

function constraint_abs_branch_flow_ordering(pm::AbstractPowerModel, l::Int; nw::Int=pm.cnw)
    branch = ref(pm, nw, :branch, l)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (l, f_bus, t_bus)
    t_idx = (l, t_bus, f_bus)
    u_K = var(pm,nw,:u_K)
    u_ord_aux = var(pm,nw,:u_ord_aux,l)

    upf0m = var(pm,nw, :up_br0)[f_idx,0]
    upt0m = var(pm,nw, :up_br0)[t_idx,0]
    upf0p = var(pm,nw, :up_br0)[f_idx,1]
    upt0p = var(pm,nw, :up_br0)[t_idx,1]
    upf1m = var(pm,nw, :up_br1)[f_idx,0]
    upt1m = var(pm,nw, :up_br1)[t_idx,0]
    upf1p = var(pm,nw, :up_br1)[f_idx,1]
    upt1p = var(pm,nw, :up_br1)[t_idx,1]

    cref=nothing

    
    if haskey( var(pm,nw), :p ) && haskey( var(pm,nw), :q )
      uqf0m = var(pm,nw, :uq_br0)[f_idx,0]
      uqt0m = var(pm,nw, :uq_br0)[t_idx,0]
      uqf0p = var(pm,nw, :uq_br0)[f_idx,1]
      uqt0p = var(pm,nw, :uq_br0)[t_idx,1]
      uqf1m = var(pm,nw, :uq_br1)[f_idx,0]
      uqt1m = var(pm,nw, :uq_br1)[t_idx,0]
      uqf1p = var(pm,nw, :uq_br1)[f_idx,1]
      uqt1p = var(pm,nw, :uq_br1)[t_idx,1]
      cref=JuMP.@constraint(pm.model, -(upf0m + upf0p + upt0m + upt0p + uqf0m + uqf0p + uqt0m + uqt0p) 
		+ (upf1m + upf1p + upt1m + upt1p + uqf1m + uqf1p + uqt1m + uqt1p) + u_ord_aux + u_K >= 0)
    elseif haskey( var(pm,nw), :p )
      cref=JuMP.@constraint(pm.model, -(upf0m + upf0p + upt0m + upt0p) + (upf1m + upf1p + upt1m + upt1p ) + u_ord_aux + u_K >= 0)
    end
    return cref
end

#=
"`[rate_a, p[f_idx], q[f_idx]] in SecondOrderCone`"
function constraint_thermal_limit_from_psd(pm::AbstractConicModels, i::Int)
    if !haskey(con(pm, pm.cnw), :sm_fr)
        con(pm, pm.cnw)[:sm_fr] = Dict{Int,Any}() # note this can be a constraint or a variable bound
    end

    branch = ref(pm, pm.cnw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    if haskey(branch, "rate_a")
        rate_a = branch["rate_a"]
        p_fr = var(pm, pm.cnw, :p, f_idx)
        q_fr = var(pm, pm.cnw, :q, f_idx)
        JuMP.@constraint(pm.model, Symmetric([rate_a p_fr q_fr; p_fr rate_a 0; q_fr 0 rate_a]) in JuMP.PSDCone())
    end
end

"`[rate_a, p[t_idx], q[t_idx]] in SecondOrderCone`"
function constraint_thermal_limit_to_psd(pm::AbstractConicModels, i::Int)
    if !haskey(con(pm, pm.cnw), :sm_to)
        con(pm, pm.cnw)[:sm_to] = Dict{Int,Any}() # note this can be a constraint or a variable bound
    end

    branch = ref(pm, pm.cnw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)
    if haskey(branch, "rate_a")
        rate_a = branch["rate_a"]
        p_to = var(pm, pm.cnw, :p, t_idx)
        q_to = var(pm, pm.cnw, :q, t_idx)
        JuMP.@constraint(pm.model, Symmetric([rate_a p_to q_to; p_to rate_a 0; q_to 0 rate_a]) in JuMP.PSDCone())
    end
end

"checks if a sufficient number of variables exist for the given keys collection"
function _check_var_keys(vars, keys, var_name, comp_name)
    if length(vars) < length(keys)
        error(_LOGGER, "$(var_name) decision variables appear to be missing for $(comp_name) components")
    end
end

=#
