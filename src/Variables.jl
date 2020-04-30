using JuMP
include("utils.jl")

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
    for l in ids(pm,nw,:branch)
        branch = ref(pm, nw, :branch, l)
        f_bus,t_bus = branch["f_bus"],branch["t_bus"]
        f_idx,t_idx = (l, f_bus, t_bus),(l, t_bus, f_bus)
        if haskey(var(pm,nw),:p)
            p_fr = var(pm, nw,:p, f_idx)
            p_to = var(pm, nw,:p, t_idx)
            if typeof(p_fr)==VariableRef
                if JuMP.has_upper_bound(p_fr)
                    @assert JuMP.upper_bound(p_fr) >= 0
                    JuMP.set_upper_bound(pt_br[f_idx,1], JuMP.upper_bound(p_fr))
                    JuMP.delete_upper_bound(p_fr)
                end
                if JuMP.has_lower_bound(p_fr)
                    @assert JuMP.lower_bound(p_fr) <= 0
                    JuMP.set_upper_bound(pt_br[f_idx,0], abs(JuMP.lower_bound(p_fr)))
                    JuMP.delete_lower_bound(p_fr)
                end
            else
                #TODO
            end
            if typeof(p_to)==VariableRef
                if JuMP.has_upper_bound(p_to)
                    @assert JuMP.upper_bound(p_to) >= 0
                    JuMP.set_upper_bound(pt_br[t_idx,1], JuMP.upper_bound(p_to))
                    JuMP.delete_upper_bound(p_to)
                end
                if JuMP.has_lower_bound(p_to)
                    @assert JuMP.lower_bound(p_to) <= 0
                    JuMP.set_upper_bound(pt_br[t_idx,0], abs(JuMP.lower_bound(p_to)))
                    JuMP.delete_lower_bound(p_to)
                end
            else
                #TODO
            end
        end
        if haskey(var(pm,nw),:q)
            q_fr = var(pm, nw, :q, f_idx)
            q_to = var(pm, nw, :q, t_idx)
            if typeof(q_fr)==VariableRef
                if JuMP.has_upper_bound(q_fr)
                    @assert JuMP.upper_bound(q_fr) >= 0
                    JuMP.set_upper_bound(qt_br[f_idx,1], JuMP.upper_bound(q_fr))
                    JuMP.delete_upper_bound(q_fr)
                end
                if JuMP.has_lower_bound(q_fr)
                    @assert JuMP.lower_bound(q_fr) <= 0
                    JuMP.set_upper_bound(qt_br[f_idx,0], abs(JuMP.lower_bound(q_fr)))
                    JuMP.delete_lower_bound(q_fr)
                end
            else
                #TODO
            end
            if typeof(q_to)==VariableRef
                if JuMP.has_upper_bound(q_to)
                    @assert JuMP.upper_bound(q_to) >= 0
                    JuMP.set_upper_bound(qt_br[t_idx,1], JuMP.upper_bound(q_to))
                    JuMP.delete_upper_bound(q_to)
                end
                if JuMP.has_lower_bound(q_to)
                    @assert JuMP.lower_bound(q_to) <= 0
                    JuMP.set_upper_bound(qt_br[t_idx,0], abs(JuMP.lower_bound(q_to)))
                    JuMP.delete_lower_bound(q_to)
                end
            else
                #TODO
            end
        end
    end
    replaceThermalLineLimits(pm)
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

