#=
"
Utility functions for branch-and-bound
9 April 2020
Brian Dandurand
"
=#

include("psd_utils.jl")
using JuMP, MathOptInterface
using Mosek, MosekTools
using CPLEX
using JuMP

function solveNodeSP_ECP(model_info, ndata)
    empty!(model_info["x_soln"])
    for l in ndata["inactive_branches"]
        model_info["x_soln"][l]=1
    end
    for l in ndata["protected_branches"]
        model_info["x_soln"][l]=0
    end
    solveSP(model_info; fix_x=true, record_x=true, compute_projection=true, compute_psd_dual=false, io=devnull)

    ndata["x_soln"] = copy(model_info["x_soln"])
    ndata["bound_value"] = model_info["opt_val"]
end

function applyPrimalHeuristic(model_info)
    model = model_info["model"]
    branch_ids = model_info["branch_ids"]
    unfix_vars(model,branch_ids)
    idx=sortperm(sort(collect(keys(model_info["x_soln"]))), by=kk->model_info["x_soln"][kk]; rev=true)
    n_idx = length(idx)
    n_attack = min(model_info["attacker_budget"],n_idx)
    model_info["x_soln_str"]=""
    for kk in idx
        l=idx[kk]
        if kk <= n_attack && model_info["x_soln"][l] > 0
            model_info["heur_x_soln"][l] = 1
            fix(variable_by_name(model,"x[$l]_1"), 1; force=true)
            model_info["x_soln_str"] = string(model_info["x_soln_str"]," $l")
        else
            model_info["heur_x_soln"][l] = 0
            fix(variable_by_name(model,"x[$l]_1"), 0; force=true)
        end
    end
    JuMP.optimize!(model)
    model_info["opt_val"]=JuMP.objective_value(model)
    model_info["solve_status"]=JuMP.termination_status(model)
    return model_info["opt_val"]
end

function deleteNodesByBound(BnBTree,bdval)
    for bb in keys(BnBTree)
        if BnBTree[bb]["bound_value"] <= bdval
            delete!(BnBTree,bb)
        end
    end
end

function findNextNode(BnBTree::Dict{Int64,Dict{String,Any}})
  nodekey=-1
  weakestUBVal = -1e20
  for kk in keys(BnBTree) 
    if BnBTree[kk]["bound_value"] > weakestUBVal
        nodekey = kk
        weakestUBVal = BnBTree[kk]["bound_value"]
    end
  end
  return nodekey,weakestUBVal
end

function findNextIndex(undecided_branches,x_vals)
  return reduce((k1, k2) -> abs(x_vals[k1]-0.5) <= abs(x_vals[k2]-0.5) ? k1 : k2, undecided_branches)
end

function printNode(node::Dict{String,Any}; pretext="",posttext="")
    println(pretext,"Node ",node["node_key"]," has bound ",round(node["bound_value"];digits=5),
        ": inactive=",node["inactive_branches"]," protected: ",node["protected_branches"],posttext)
end
