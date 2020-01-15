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
    m = DualizeModel(m)
    return m
end

function DualizeModel(minmax_model_pm::AbstractPowerModel)
    io = open("minmax.txt","w")
    println(io,minmax_model_pm.model)
    close(io)
    maxmin_model = dualize(minmax_model_pm.model)
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
    fn_base=string(minmax_model_pm.data["name"],".cbf")
println(fn_base)
    mathoptformat_model = MathOptInterface.FileFormats.CBF.Model()
    MOI.copy_to(MOI.Bridges.full_bridge_optimizer(mathoptformat_model,Float64), backend(maxmin_model))
    #MOI.copy_to(mathoptformat_model, backend(maxmin_model))
    MOI.write_to_file(mathoptformat_model, fn_base)
#=
    fn_base=string(minmax_model_pm.data["name"],".")
println(fn_base)
    mathoptformat_model = MathOptInterface.FileFormats.CBF.Model()
    MOI.copy_to(MOI.Bridges.full_bridge_optimizer(mathoptformat_model,Float64), backend(maxmin_model))
    #MOI.copy_to(mathoptformat_model, backend(maxmin_model))
    MOI.write_to_file(mathoptformat_model, fn_base)
=#

    return maxmin_model
end

function MinimaxOPFModel(pm_data, powerform)
    if powerform == SOCWRConicPowerModel || powerform == SDPWRMPowerModel || powerform == SparseSDPWRMPowerModel 
        println("Prototyping Algo with WRConic Forms")
        pm = build_model(pm_data, powerform, WRConicPost_PF_Minmax)
    else
	println("Do nothing")
    end

    return pm
end

function PF_FeasModel(pm_data, powerform)
    if powerform == SOCWRConicPowerModel || powerform == SDPWRMPowerModel || powerform == SparseSDPWRMPowerModel 
        println("Prototyping Algo with WRConic Forms")
        pm = build_model(pm_data, powerform, WRConicPost_PF)
        objective_feasibility_problem(pm)
    else
	println("Do nothing")
    end
    return pm
end

function WRConicPost_PF_Minmax(pm::AbstractPowerModel)
    WRConicPost_PF(pm)
    variable_ordering_auxiliary(pm)
    con(pm, pm.cnw, pm.ccnd)[:x] = Dict{Int,JuMP.ConstraintRef}()
    undecided_branches = filter(l->!(l in pm.data["protected_branches"] || l in pm.data["inactive_branches"]), ids(pm,pm.cnw,:branch))
    for l in undecided_branches
        con(pm, pm.cnw, pm.ccnd)[:x][l] = constraint_abs_branch_flow_ordering(pm, l)
        JuMP.set_name(con(pm, pm.cnw, pm.ccnd)[:x][l],"x[$l]")  
    end
    objective_minmax_problem(pm)
end

function WRConicPost_PF(pm::AbstractPowerModel)
    variable_voltage(pm)
    variable_generation(pm)
    variable_branch_flow(pm)
    variable_dcline_flow(pm)
    remove_infinity_bnds(pm)

    # Add new variables
    #variable_bus_slacks(pm)
    variable_branch_flow_slacks(pm)

    # Post objective function for testing fisibility problem
    # objective_feasibility_problem(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        #constraint_power_balance_slacks(pm, i)
        constraint_power_balance(pm, i)
    end

    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_disc_ub] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_disc_lb] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_disc_ub] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_disc_lb] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_disc_ub] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_disc_lb] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_disc_ub] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_disc_lb] = Dict{Int,JuMP.ConstraintRef}()

    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_ub] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_lb] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_ub] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_lb] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_ub] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_lb] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_ub] = Dict{Int,JuMP.ConstraintRef}()
    con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_lb] = Dict{Int,JuMP.ConstraintRef}()
    for l in ids(pm, :branch)
        cref1,cref2,cref3,cref4 = constraint_ohms_yt_from_slacks(pm, l)
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_disc_ub][l] = cref1
        JuMP.set_name(cref1,"lambda_f+[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_disc_lb][l] = cref2
        JuMP.set_name(cref2,"lambda_f-[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_disc_ub][l] = cref3
        JuMP.set_name(cref3,"mu_f+[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_disc_lb][l] = cref4
        JuMP.set_name(cref4,"mu_f-[$l]")  

        cref1,cref2,cref3,cref4 = constraint_ohms_yt_to_slacks(pm, l)
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_disc_ub][l] = cref1
        JuMP.set_name(cref1,"lambda_t+[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_disc_lb][l] = cref2
        JuMP.set_name(cref2,"lambda_t-[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_disc_ub][l] = cref3
        JuMP.set_name(cref3,"mu_t+[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_disc_lb][l] = cref4
        JuMP.set_name(cref4,"mu_t-[$l]")  

        ref_p1,ref_p2,ref_p3,ref_p4,ref_q1,ref_q2,ref_q3,ref_q4 = constraint_def_abs_flow_values(pm, l)
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_ub][l] = ref_p1
        JuMP.set_name(ref_p2,"pi_f+[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_fr_lb][l] = ref_p2
        JuMP.set_name(ref_p1,"pi_f-[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_ub][l] = ref_p3
        JuMP.set_name(ref_p4,"pi_t+[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_pflow_to_lb][l] = ref_p4
        JuMP.set_name(ref_p3,"pi_t-[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_ub][l] = ref_q1
        JuMP.set_name(ref_q1,"phi_f+[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_fr_lb][l] = ref_q2
        JuMP.set_name(ref_q2,"phi_f-[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_ub][l] = ref_q3
        JuMP.set_name(ref_q4,"phi_t+[$l]")  
        con(pm, pm.cnw, pm.ccnd)[:abs_qflow_to_lb][l] = ref_q4
        JuMP.set_name(ref_q3,"phi_t-[$l]")  


        constraint_voltage_angle_difference(pm, l)
	if pm isa AbstractSDPWRMModel
	  ## leave only the bounds on the branch flows set during initialization when there is need to avoid SOC or quadratic constraints
           #constraint_thermal_limit_from(pm, l)
           #constraint_thermal_limit_to(pm, l)
	else
          constraint_thermal_limit_from(pm, l)
          constraint_thermal_limit_to(pm, l)
	end
    end

    for l in ids(pm, :dcline)
        constraint_dcline(pm, l)  # DO WE NEED TO TREAT THESE CONSTRAINTS DIFFERENTLY
    end
end



end # module
