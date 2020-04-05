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

      pd_fm = var(pm,nw, :pd_br)[f_idx,0]
      pd_fp = var(pm,nw, :pd_br)[f_idx,1]
      pd_tm = var(pm,nw, :pd_br)[t_idx,0]
      pd_tp = var(pm,nw, :pd_br)[t_idx,1]

      pt_fm = var(pm,nw, :pt_br)[f_idx,0]
      pt_fp = var(pm,nw, :pt_br)[f_idx,1]
      pt_tm = var(pm,nw, :pt_br)[t_idx,0]
      pt_tp = var(pm,nw, :pt_br)[t_idx,1]

      cref_pf=JuMP.@constraint(pm.model, p_fr - pd_fm + pd_fp + pt_fm - pt_fp == 0)
      cref_pt=JuMP.@constraint(pm.model, p_to - pd_tm + pd_tp + pt_tm - pt_tp == 0)
    end
    
    cref_qf,cref_qt=nothing,nothing
    if haskey( var(pm,nw), :q )
      q_fr = var(pm, nw, :q, f_idx)
      q_to = var(pm, nw, :q, t_idx)

      qd_fm = var(pm,nw, :qd_br)[f_idx,0]
      qd_fp = var(pm,nw, :qd_br)[f_idx,1]
      qd_tm = var(pm,nw, :qd_br)[t_idx,0]
      qd_tp = var(pm,nw, :qd_br)[t_idx,1]

      qt_fm = var(pm,nw, :qt_br)[f_idx,0]
      qt_fp = var(pm,nw, :qt_br)[f_idx,1]
      qt_tm = var(pm,nw, :qt_br)[t_idx,0]
      qt_tp = var(pm,nw, :qt_br)[t_idx,1]

      cref_qf=JuMP.@constraint(pm.model, q_fr - qd_fm + qd_fp + qt_fm - qt_fp == 0)
      cref_qt=JuMP.@constraint(pm.model, q_to - qd_tm + qd_tp + qt_tm - qt_tp == 0)
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

    pt_fm = var(pm,nw, :pt_br)[f_idx,0]
    pt_tm = var(pm,nw, :pt_br)[t_idx,0]
    pt_fp = var(pm,nw, :pt_br)[f_idx,1]
    pt_tp = var(pm,nw, :pt_br)[t_idx,1]
    pd_fm = var(pm,nw, :pd_br)[f_idx,0]
    pd_tm = var(pm,nw, :pd_br)[t_idx,0]
    pd_fp = var(pm,nw, :pd_br)[f_idx,1]
    pd_tp = var(pm,nw, :pd_br)[t_idx,1]

    cref=nothing

    
    if haskey( var(pm,nw), :p ) && haskey( var(pm,nw), :q )
        qt_fm = var(pm,nw, :qt_br)[f_idx,0]
        qt_tm = var(pm,nw, :qt_br)[t_idx,0]
        qt_fp = var(pm,nw, :qt_br)[f_idx,1]
        qt_tp = var(pm,nw, :qt_br)[t_idx,1]
        qd_fm = var(pm,nw, :qd_br)[f_idx,0]
        qd_tm = var(pm,nw, :qd_br)[t_idx,0]
        qd_fp = var(pm,nw, :qd_br)[f_idx,1]
        qd_tp = var(pm,nw, :qd_br)[t_idx,1]
        cref=JuMP.@constraint(pm.model, -(pt_fm + pt_fp + pt_tm + pt_tp + qt_fm + qt_fp + qt_tm + qt_tp) 
		    + (pd_fm + pd_fp + pd_tm + pd_tp + qd_fm + qd_fp + qd_tm + qd_tp) + u_ord_aux + u_K >= 0)
    elseif haskey( var(pm,nw), :p )
        cref=JuMP.@constraint(pm.model, -(pt_fm + pt_fp + pt_tm + pt_tp) + (pd_fm + pd_fp + pd_tm + pd_tp ) + u_ord_aux + u_K >= 0)
    end
    return cref
end

