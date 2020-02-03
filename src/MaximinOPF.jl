module MaximinOPF
using PowerModels
using Dualization
using MathOptInterface

include("Variables.jl")
include("Objectives.jl")
include("Constraints.jl")
greet() = print("Hello World!")

function MaximinOPFModel(pm_data, powerform)
    println("Hello MaximinOPFModel")
    m = MinimaxOPFModel(pm_data, powerform)

    #Test Out
    io = open(string(pm_data["name"],".out"), "w")
    println(io,"name: ", pm_data["name"])
    println(io,"attacker_budget: ", pm_data["attacker_budget"])
    println(io,"inactive_branches: ", pm_data["inactive_branches"])
    println(io,"protected_branches: ", pm_data["protected_branches"])
    println(io, "Model:")
    println(io, m.model)
    close(io)

    if powerform == ACRPowerModel
        println("ACRPowerModel is returned only MinMaxModel")
    else
        m = DualizeModel(m)  
    end

    
    return m
end


function MinimaxOPFModel(pm_data, powerform)
    if powerform == SOCWRConicPowerModel || powerform == SDPWRMPowerModel || powerform == SparseSDPWRMPowerModel || powerform == ACRPowerModel
      pm = instantiate_model(pm_data, powerform, WRConicPost_PF_Minmax)
    else
        println("Not Supported Power Model Option")
    end
    return pm
end

function WRConicPost_PF_Minmax(pm::AbstractPowerModel)
    WRConicPost_PF(pm)
    variable_ordering_auxiliary(pm)
    con(pm, pm.cnw, pm.ccnd)[:x] = Dict{Int,JuMP.ConstraintRef}()
    for l in pm.data["undecided_branches"]
        con(pm, pm.cnw, pm.ccnd)[:x][l] = constraint_abs_branch_flow_ordering(pm, l)
        JuMP.set_name(con(pm, pm.cnw, pm.ccnd)[:x][l],"x[$l]")  
    end
    objective_minmax_problem(pm)
end

function PF_FeasModel(pm_data, powerform, x_vals=Dict{Int64,Float64}() )
    if powerform == SOCWRConicPowerModel || powerform == SDPWRMPowerModel || powerform == SparseSDPWRMPowerModel || powerform == ACRPowerModel
        pm = instantiate_model(pm_data, powerform, WRConicPost_PF)
        for l in ids(pm,pm.cnw,:branch)
          if !haskey(x_vals,l)
	    if l in pm_data["inactive_branches"]
		x_vals[l]=1
	    else
		x_vals[l]=0
	    end
          end
        end
        objective_feasibility_problem(pm,x_vals)
    else
	println("Do nothing")
    end
    return pm
end

function WRConicPost_PF(pm::AbstractPowerModel)
    if !haskey(pm.data,"inactive_branches")
        pm.data["inactive_branches"] = []
    end
    if !haskey(pm.data,"protected_branches")
        pm.data["protected_branches"] = []
    end
    pm.data["undecided_branches"] = filter(l->!(l in pm.data["protected_branches"] || l in pm.data["inactive_branches"]), ids(pm,pm.cnw,:branch))
    variable_voltage(pm)
    variable_generation(pm)
    variable_branch_flow(pm)
    variable_dcline_flow(pm)
    constraint_model_voltage(pm)
    remove_infinity_bnds(pm)

    # Add new variables
    variable_branch_flow_slacks(pm)


    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_disc] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_disc] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_disc] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_disc] = Dict{Int,JuMP.ConstraintRef}()

    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to] = Dict{Int,JuMP.ConstraintRef}()

    for l in setdiff(ids(pm, :branch),pm.data["inactive_branches"])
        cref1,cref2 = constraint_ohms_yt_from_slacks(pm, l)
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_disc][l] = cref1
        JuMP.set_name(cref1,"lambda_f[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_disc][l] = cref2
        JuMP.set_name(cref2,"mu_f[$l]")  

        cref1,cref2 = constraint_ohms_yt_to_slacks(pm, l)
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_disc][l] = cref1
        JuMP.set_name(cref1,"lambda_t[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_disc][l] = cref2
        JuMP.set_name(cref2,"mu_t[$l]")  
    end
#equal_to_constraints = all_constraints(pm.model, GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64})  ###USE FOR VALIDATION?
#println(equal_to_constraints)

    for l in setdiff(ids(pm, :branch),pm.data["protected_branches"])
        ref_p1,ref_p3,ref_q1,ref_q3 = constraint_def_abs_flow_values(pm, l)
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr][l] = ref_p1
        JuMP.set_name(ref_p1,"pi_f[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to][l] = ref_p3
        JuMP.set_name(ref_p3,"pi_t[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr][l] = ref_q1
        JuMP.set_name(ref_q1,"phi_f[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to][l] = ref_q3
        JuMP.set_name(ref_q3,"phi_t[$l]")  
    end

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end


    for l in ids(pm, :branch)
        constraint_voltage_angle_difference(pm, l)
        constraint_thermal_limit_from(pm, l)
        constraint_thermal_limit_to(pm, l)
    end

    for l in ids(pm, :dcline)
        constraint_dcline(pm, l)  # DO WE NEED TO TREAT THESE CONSTRAINTS DIFFERENTLY
    end
end

function DualizeModel(minmax_model_pm::AbstractPowerModel)
    maxmin_model=DualizeMinmaxModel(minmax_model_pm::AbstractPowerModel)
    for l in ids(minmax_model_pm, :branch)
        if !(l in minmax_model_pm.data["protected_branches"] || l in minmax_model_pm.data["inactive_branches"])
            if has_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                delete_lower_bound(variable_by_name(maxmin_model,"x[$l]_1"))
            end
            if has_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
                delete_upper_bound(variable_by_name(maxmin_model,"x[$l]_1"))
            end
            JuMP.set_integer(variable_by_name(maxmin_model,"x[$l]_1"))
        end
    end
    return maxmin_model
end
function DualizeMinmaxModel(minmax_model_pm::AbstractPowerModel)
    dualizable_minmax_model = MOI.Utilities.Model{Float64}()
    bridged_model = MOI.Bridges.Constraint.Square{Float64}(dualizable_minmax_model)
    MOI.copy_to(bridged_model,backend(minmax_model_pm.model))    
    dualized_minmax_problem = dualize(dualizable_minmax_model)
    dualized_minmax_model = JuMP.Model() 
    MOI.copy_to(dualized_minmax_model,dualized_minmax_problem.dual_model)
    return dualized_minmax_model
end


function write_to_cbf(model,fn_base::String)
    JuMP.write_to_file( model, string(fn_base,".cbf"), format = MOI.FileFormats.FORMAT_CBF)
end
function write_to_cbf_scip(model,fn_base::String)
    #BRIDGE SOC CONSTRAINTS
    model_psd_moi = MOI.Utilities.Model{Float64}()
    bridged_model = MOI.Bridges.Constraint.SOCtoPSD{Float64}(model_psd_moi)
    MOI.copy_to(bridged_model,backend(model))
    model_psd = JuMP.Model()
    MOI.copy_to(model_psd,model_psd_moi)
    JuMP.write_to_file( model_psd, string(fn_base,"_scip",".cbf"), format = MOI.FileFormats.FORMAT_CBF)
end

end # module
