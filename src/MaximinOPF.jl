module MaximinOPF
using PowerModels
using Dualization
using MathOptInterface
using SCS
include("Variables.jl")
include("Objectives.jl")
include("Constraints.jl")
greet() = print("Hello World!")

function MaximinOPFModel(pm_data, powerform)
    println("Hello MaximinOPFModel")
    m = MinimaxOPFModel(pm_data, powerform)
    m = Minmax_to_Maxmin(m)
    return m
end

function Minmax_to_Maxmin(minmax_model_pm::AbstractPowerModel)
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
    
    io = open("dualized_minmax.txt","w")
    println(io,dualized_minmax_model)
    close(io)
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
    io = open("temp","w")
    println(io,model_psd)
    close(io)
end

function MinimaxOPFModel(pm_data, powerform)
    if powerform == SOCWRConicPowerModel || powerform == SDPWRMPowerModel || powerform == SparseSDPWRMPowerModel || powerform == ACRPowerModel
        pm = instantiate_model(pm_data, powerform, WRConicPost_PF_Minmax)
    else
	println("Do nothing")
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

"After the minmax model is solved, we can query various measures of sensitivity"
#=
function computeSensitivityData(pm::AbstractPowerModel)
    pm.data["x_vals"]=Dict{Int64,Float64}()
    pm.data["sg"]=Dict{Int64,Float64}()
    pm.data["uoa_vals"]=Dict{Int64,Float64}()
    if haskey( var(pm,pm.cnw,pm.ccnd), :u_K )
      pm.data["uK_val"] = JuMP.value(var(pm,pm.cnw,pm.ccnd)[:u_K])
    end
    for l in ids(pm, :branch)
      branch = ref(pm, pm.cnw, :branch, l)
      f_bus = branch["f_bus"]
      t_bus = branch["t_bus"]
      f_idx = (l, f_bus, t_bus)
      t_idx = (l, t_bus, f_bus)
      upf0 = abs(JuMP.value(var(pm, pm.cnw, pm.ccnd,:p)[f_idx]))
      upt0 = abs(JuMP.value(var(pm, pm.cnw, pm.ccnd,:p)[t_idx]))
      uqf0 = abs(JuMP.value(var(pm, pm.cnw, pm.ccnd,:q)[f_idx]))
      uqt0 = abs(JuMP.value(var(pm, pm.cnw, pm.ccnd,:q)[t_idx]))
      upf1 = abs(JuMP.value(var(pm, pm.cnw, pm.ccnd,:p_expr)[f_idx]) - JuMP.value(var(pm, pm.cnw, pm.ccnd,:p)[f_idx]))
      upt1 = abs(JuMP.value(var(pm, pm.cnw, pm.ccnd,:p_expr)[t_idx]) - JuMP.value(var(pm, pm.cnw, pm.ccnd,:p)[t_idx]))
      uqf1 = abs(JuMP.value(var(pm, pm.cnw, pm.ccnd,:q_expr)[f_idx]) - JuMP.value(var(pm, pm.cnw, pm.ccnd,:q)[f_idx]))
      uqt1 = abs(JuMP.value(var(pm, pm.cnw, pm.ccnd,:q_expr)[t_idx]) - JuMP.value(var(pm, pm.cnw, pm.ccnd,:q)[t_idx]))
      pm.data["sg"][l] = upf0+upt0+uqf0+uqt0-upf1-upt1-uqf1-uqt1	        

      if !(l in pm.data["undecided_branches"])
        if l in pm.data["protected_branches"]
          pm.data["x_vals"][l] = 0
        end
        if l in pm.data["inactive_branches"]
          pm.data["x_vals"][l] = 1
        end
      else
        if haskey(con(pm,pm.cnw,pm.ccnd), :x )
          pm.data["x_vals"][l]=abs(JuMP.dual(con(pm, pm.cnw, pm.ccnd)[:x][l])) 
	end
        pm.data["x_vals"][l] = max(0,pm.data["x_vals"][l])
        pm.data["x_vals"][l] = min(pm.data["x_vals"][l],1)
	if haskey(var(pm,pm.cnw,pm.ccnd), :u_ord_aux)
          pm.data["uoa_vals"][l] = JuMP.value(var(pm,pm.cnw,pm.ccnd)[:u_ord_aux][l])
	end
      end
    end
    println("sg: ",pm.data["sg"])
    println("x_vals: ",pm.data["x_vals"])
    println("uoa_vals: ",pm.data["uoa_vals"])
    println("uK_val: ",pm.data["uK_val"])
end
=#

function PF_FeasModel(pm_data, powerform, x_vals=Dict{Int64,Float64}() )
    if powerform == SOCWRConicPowerModel || powerform == SDPWRMPowerModel || powerform == SparseSDPWRMPowerModel 
        pm = instantiate_model(pm_data, powerform, WRConicPost_PF)
        for l in ids(pm,pm.cnw,:branch)
          if !haskey(x_vals,l)
	    x_vals[l]=0
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

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

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


    for l in ids(pm, :branch)
        constraint_voltage_angle_difference(pm, l)
	if pm isa AbstractSDPWRMModel
	  ## leave only the bounds on the branch flows set during initialization when there is need to avoid SOC or quadratic constraints
	   #println("Using psd form of thermal line limits")
           #constraint_thermal_limit_from_psd(pm, l)
           #constraint_thermal_limit_to_psd(pm, l)
           constraint_thermal_limit_from(pm, l)
           constraint_thermal_limit_to(pm, l)
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
