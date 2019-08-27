
function setupDualModel(opfdata, model, useBin=false, HEUR=0)
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

    # each x[l] is either 0 (line l active) or 1 (line l is cut)
      if useBin
        @variable(model, x[l=L], Bin, start=0)
      else
        @variable(model, 0 <= x[l=L] <= 1, start=0)
      end
    # dual multipliers associated with active power flow balance constraints
      @variable(model, -1 <= α[i=N] <= 1, start=0)
    # dual multipliers associated with active power flow balance constraints
      @variable(model, -1 <= β[i=N] <= 1, start=0)
    # dual multipliers associated with the voltage magnitude bounds
      @variable(model, δ[i=N] >= 0); @variable(model, γ[i=N] >= 0); @constraint(model, [i=N], δ[i]+γ[i] <= 1)
    # dual multipliers associated with the power generation bounds
      @variable(model, ζpUB[g=G] >=0); @variable(model, ζpLB[g=G] >=0); @variable(model, ζqUB[g=G] >=0); @variable(model, ζqLB[g=G] >=0)

    # constraints needed for terms in pG[g] and qG[g] in the Lagrangian to vanish (needed for dual feasibility)
      for i in N
        for g in BusGeners[i]
          @constraint(model, -α[i] + ζpUB[g] - ζpLB[g] == 0 )
          @constraint(model, -β[i] + ζqUB[g] - ζqLB[g] == 0 )
        end
      end

  @constraint(model, sum(x[l] for l in L) <= K)

 # McCormick inequalities enforcing bilinear equalities
    # auxiliary dual variables due to McCormick reformulation of cross terms appearing in the Lagrangian
      @variable(model, λF[l=L], start=0); @variable(model, λT[l=L], start=0); @variable(model, μF[l=L], start=0); @variable(model, μT[l=L], start=0)
    @constraint(model, AMcf1[l in L], α[fromBus[l]] - x[l] <= λF[l]); @constraint(model, AMcf2[l in L], α[fromBus[l]] + x[l] >= λF[l])
    @constraint(model, AMcf3[l in L], -(1 - x[l]) <= λF[l]); @constraint(model, AMcf4[l in L], (1 - x[l]) >= λF[l])
    @constraint(model, AMct1[l in L], α[toBus[l]] - x[l] <= λT[l]); @constraint(model, AMct2[l in L], α[toBus[l]] + x[l] >= λT[l])
    @constraint(model, AMct3[l in L], -(1 - x[l]) <= λT[l]); @constraint(model, AMct4[l in L], (1 - x[l]) >= λT[l])

    @constraint(model, BMcf1[l in L], β[fromBus[l]] - x[l] <= μF[l]); @constraint(model, BMcf2[l in L], β[fromBus[l]] + x[l] >= μF[l])
    @constraint(model, BMcf3[l in L], -(1 - x[l]) <= μF[l]); @constraint(model, BMcf4[l in L], (1 - x[l]) >= μF[l])
    @constraint(model, BMct1[l in L], β[toBus[l]] - x[l] <= μT[l]); @constraint(model, BMct2[l in L], β[toBus[l]] + x[l] >= μT[l])
    @constraint(model, BMct3[l in L], -(1 - x[l]) <= μT[l]); @constraint(model, BMct4[l in L], (1 - x[l]) >= μT[l])

  # Lagrangian objective, after vanishing terms are removed under the assumption of dual feasibility
    @objective(model, Max, sum(ζpLB[g]*opfdata.Pmin[g] - ζpUB[g]*opfdata.Pmax[g] + ζqLB[g]*opfdata.Qmin[g] - ζqUB[g]*opfdata.Qmax[g]  for g in G)
    	+ sum( γ[i]*opfdata.Wmin[i]-δ[i]*opfdata.Wmax[i] + α[i]*opfdata.PD[i] + β[i]*opfdata.QD[i] for i in N) )
  ### END DEFINING THE LaGRANGIAN DUAL PROBLEM

    if HEUR == 1
      @constraint(mMP, LambdaMuConstr1[l in L], λF[l]*Y["ftI"][l] - λT[l]*Y["tfI"][l] + μF[l]*Y["ftR"][l] - μT[l]*Y["tfR"][l] == 0.0)
    elseif HEUR == 2
      @constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l]  == 0)
      @constraint(mMP, muFequalsT[l in L], μF[l] - μT[l]  == 0)
    elseif HEUR == 3
      lRelax=rand(Bool,nlines)
      @show lRelax
      for l in L
        if lRelax[l]
          @constraint(mMP, λF[l] - λT[l]  == 0)
          @constraint(mMP, μF[l] - μT[l]  == 0)
        end
      end
    end
end

function setupDualModelObj(opfdata, model)
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata
  # Lagrangian objective, after vanishing terms are removed under the assumption of dual feasibility
    @objective(model, Max, sum(ζpLB[g]*opfdata.Pmin[g] - ζpUB[g]*opfdata.Pmax[g] + ζqLB[g]*opfdata.Qmin[g] - ζqUB[g]*opfdata.Qmax[g]  for g in G)
    	+ sum( γ[i]*opfdata.Wmin[i]-δ[i]*opfdata.Wmax[i] + α[i]*opfdata.PD[i] + β[i]*opfdata.QD[i] for i in N) )
  ### END DEFINING THE LaGRANGIAN DUAL PROBLEM
end
