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

