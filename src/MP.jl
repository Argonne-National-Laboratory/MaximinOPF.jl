#=
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")

CP,PROX0,PROX,LVL1,LVL2,LVLINF,FEAS=0,1,2,3,4,5,6


function createBasicMP(opfdata,nodeinfo,ctr,K,CTR_PARAM)
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = 1:opfdata.nbuses, 1:opfdata.nlines, 1:opfdata.ngens
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC

  # CREATE MODEL
    mMP = Model(with_optimizer(Ipopt.Optimizer))

    # each x[l] is either 0 (line l active) or 1 (line l is cut)
      @variable(mMP, nodeinfo.x_lbs[l] <= x[l=L] <= nodeinfo.x_ubs[l])
    # dual multipliers associated with active power flow balance constraints
      @variable(mMP, -1 <= α[i=N] <= 1)
    # dual multipliers associated with active power flow balance constraints
      @variable(mMP, -1 <= β[i=N] <= 1)
    # dual multipliers associated with the voltage magnitude bounds
      @variable(mMP, δ[i=N] >= 0)
      @variable(mMP, γ[i=N] >= 0)
      @constraint(mMP, [i=N], δ[i]+γ[i] <= 1)
    # dual multipliers associated with the power generation bounds
      @variable(mMP, 0 <= ζpUB[g=G] <= 1)
      @variable(mMP, 0 <= ζpLB[g=G] <= 1)
      @variable(mMP, 0 <= ζqUB[g=G] <= 1)
      @variable(mMP, 0 <= ζqLB[g=G] <= 1)

    # constraints needed for terms in pG[g] and qG[g] in the Lagrangian to vanish (needed for dual feasibility)
      for i in N
        for g in BusGeners[i]
          @constraint(mMP, -α[i] + ζpUB[g] - ζpLB[g] == 0 )
          @constraint(mMP, -β[i] + ζqUB[g] - ζqLB[g] == 0 )
        end
      end

  @constraint(mMP, sum(x[l] for l in L) <= K)

 # McCormick inequalities enforcing bilinear equalities
    # auxiliary dual variables due to McCormick reformulation of cross terms appearing in the Lagrangian
      @variable(mMP, λF[l=L]); @variable(mMP, λT[l=L]); @variable(mMP, μF[l=L]); @variable(mMP, μT[l=L])
    for l in L
      if nodeinfo.x_lbs[l] > 0.9999
         set_lower_bound(λF[l],0)
         set_upper_bound(λF[l],0)
         set_lower_bound(λT[l],0)
         set_upper_bound(λT[l],0)
         set_lower_bound(μF[l],0)
         set_upper_bound(μF[l],0)
         set_lower_bound(μT[l],0)
         set_upper_bound(μT[l],0)
      elseif nodeinfo.x_ubs[l] < 0.0001
	 @constraint(mMP, λF[l] - α[fromBus[l]] == 0)
	 @constraint(mMP, λT[l] - α[toBus[l]] == 0)
	 @constraint(mMP, μF[l] - β[fromBus[l]] == 0)
	 @constraint(mMP, μT[l] - β[toBus[l]] == 0)
      else
        @constraint(mMP, α[fromBus[l]] - x[l] <= λF[l]) 
        @constraint(mMP, α[fromBus[l]] + x[l] >= λF[l])
        @constraint(mMP, -(1 - x[l]) <= λF[l]) 
        @constraint(mMP,  (1 - x[l]) >= λF[l])
        @constraint(mMP, α[toBus[l]] - x[l] <= λT[l]) 
        @constraint(mMP, α[toBus[l]] + x[l] >= λT[l])
        @constraint(mMP, -(1 - x[l]) <= λT[l])
        @constraint(mMP,  (1 - x[l]) >= λT[l])

        @constraint(mMP, β[fromBus[l]] - x[l] <= μF[l])
        @constraint(mMP, β[fromBus[l]] + x[l] >= μF[l])
        @constraint(mMP, -(1 - x[l]) <= μF[l])
        @constraint(mMP,  (1 - x[l]) >= μF[l])
        @constraint(mMP, β[toBus[l]] - x[l] <= μT[l])
        @constraint(mMP, β[toBus[l]] + x[l] >= μT[l])
        @constraint(mMP, -(1 - x[l]) <= μT[l])
        @constraint(mMP,  (1 - x[l]) >= μT[l])
      end
    end
    if CTR_PARAM == PROX 
        @variable(mMP, psi)
    end

  # Lagrangian objective, after vanishing terms are removed under the assumption of dual feasibility
    lin_cfs = create_soln(opfdata)
    for i in N
      lin_cfs.α[i],lin_cfs.β[i],lin_cfs.γ[i],lin_cfs.δ[i] = opfdata.PD[i],opfdata.QD[i],opfdata.Wmin[i],-opfdata.Wmax[i]
    end
    for g in G
      lin_cfs.ζpLB[g],lin_cfs.ζpUB[g],lin_cfs.ζqLB[g],lin_cfs.ζqUB[g] = opfdata.Pmin[g],-opfdata.Pmax[g],opfdata.Qmin[g],-opfdata.Qmax[g]
    end
    @expression(mMP, linobj, sum(lin_cfs.ζpLB[g]*ζpLB[g] + lin_cfs.ζpUB[g]*ζpUB[g] + lin_cfs.ζqLB[g]*ζqLB[g] + lin_cfs.ζqUB[g]*ζqUB[g]  for g in G)
    			+ sum( lin_cfs.γ[i]*γ[i] + lin_cfs.δ[i]*δ[i] + lin_cfs.α[i]*α[i] + lin_cfs.β[i]*β[i] for i in N) 
    )

    return mMP
end

function setObjMP(opfdata,mMP,nodeinfo,ctr,CTR_PARAM)
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = 1:nbuses,1:nlines,1:ngens 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC

    α,β,γ,δ=mMP[:α],mMP[:β],mMP[:γ],mMP[:δ]
    x,λF,λT,μF,μT=mMP[:x],mMP[:λF],mMP[:λT],mMP[:μF],mMP[:μT]
    ζpLB,ζpUB,ζqLB,ζqUB=mMP[:ζpLB],mMP[:ζpUB],mMP[:ζqLB],mMP[:ζqUB]
    linobj=mMP[:linobj]


   # METHOD-SPECIFIC SETUP
    if CTR_PARAM == PROX0
      @objective(mMP, Max, linobj - 0.5*nodeinfo.tVal*(
	  sum( (ctr.soln.α[i] - α[i])^2 + (ctr.soln.β[i] - β[i])^2 + (ctr.soln.γ[i] - γ[i])^2 + (ctr.soln.δ[i] - δ[i])^2 for i in N)
	  +sum( (ctr.soln.ζpLB[g] - ζpLB[g])^2 + (ctr.soln.ζpUB[g] - ζpUB[g])^2 + (ctr.soln.ζqLB[g] - ζqLB[g])^2 + (ctr.soln.ζqUB[g] - ζqUB[g])^2 for g in G)
	  +sum( (ctr.soln.x[l]-x[l])^2 + (ctr.soln.λF[l] - λF[l])^2 + (ctr.soln.λT[l] - λT[l])^2 + (ctr.soln.μF[l] - μF[l])^2 + (ctr.soln.μT[l] - μT[l])^2 for l in L)
	)
      )
    elseif CTR_PARAM == PROX 
        @objective(mMP, Min, mMP[:psi] + 0.5*nodeinfo.tVal*(
	  sum( (ctr.soln.α[i] - α[i])^2 + (ctr.soln.β[i] - β[i])^2 + (ctr.soln.γ[i] - γ[i])^2 + (ctr.soln.δ[i] - δ[i])^2 for i in N)
	  +sum( (ctr.soln.ζpLB[g] - ζpLB[g])^2 + (ctr.soln.ζpUB[g] - ζpUB[g])^2 + (ctr.soln.ζqLB[g] - ζqLB[g])^2 + (ctr.soln.ζqUB[g] - ζqUB[g])^2 for g in G)
	  +sum( (ctr.soln.x[l]-x[l])^2 + (ctr.soln.λF[l] - λF[l])^2 + (ctr.soln.λT[l] - λT[l])^2 + (ctr.soln.μF[l] - μF[l])^2 + (ctr.soln.μT[l] - μT[l])^2 for l in L)
	  )
	)
    elseif CTR_PARAM == LVL2 || CTR_PARAM == LVL2 || CTR_PARAM == LVLINF
      @constraint(mMP, LVLConstr, linobj >= nodeinfo.nodeBd)
      if CTR_PARAM==LVL1
        @variable(mMP, alphaSlack[i=N] >= 0)
        @variable(mMP, betaSlack[i=N] >= 0)
        @variable(mMP, gammaSlack[i=N] >= 0)
        @variable(mMP, deltaSlack[i=N] >= 0)
        @constraint(mMP, alphaUBSlack[i=N], ctr.soln.α[i] - α[i] <= alphaSlack[i])
        @constraint(mMP, alphaLBSlack[i=N], ctr.soln.α[i] - α[i] >= -alphaSlack[i])
        @constraint(mMP, betaUBSlack[i=N], ctr.soln.β[i] - β[i] <= betaSlack[i])
        @constraint(mMP, betaLBSlack[i=N], ctr.soln.β[i] - β[i] >= -betaSlack[i])
        @constraint(mMP, gammaUBSlack[i=N], ctr.soln.γ[i] - γ[i] <= gammaSlack[i])
        @constraint(mMP, gammaLBSlack[i=N], ctr.soln.γ[i] - γ[i] >= -gammaSlack[i])
        @constraint(mMP, deltaUBSlack[i=N], ctr.soln.δ[i] - δ[i] <= deltaSlack[i])
        @constraint(mMP, deltaLBSlack[i=N], ctr.soln.δ[i] - δ[i] >= -deltaSlack[i])
        @objective(mMP, Max, linobj - nodeinfo.tVal*sum(alphaSlack[i] + betaSlack[i] + gammaSlack[i] + deltaSlack[i] for i in N))
      elseif CTR_PARAM==LVLINF
        @variable(mMP, Slack >= 0)
        @constraint(mMP, alphaUBSlack[i=N], ctr.soln.α[i] - α[i] <= Slack)
        @constraint(mMP, alphaLBSlack[i=N], ctr.soln.α[i] - α[i] >= -Slack)
        @constraint(mMP, betaUBSlack[i=N], ctr.soln.β[i] - β[i] <= Slack)
        @constraint(mMP, betaLBSlack[i=N], ctr.soln.β[i] - β[i] >= -Slack)
        @constraint(mMP, gammaUBSlack[i=N], ctr.soln.γ[i] - γ[i] <= Slack)
        @constraint(mMP, gammaLBSlack[i=N], ctr.soln.γ[i] - γ[i] >= -Slack)
        @constraint(mMP, deltaUBSlack[i=N], ctr.soln.δ[i] - δ[i] <= Slack)
        @constraint(mMP, deltaLBSlack[i=N], ctr.soln.δ[i] - δ[i] >= -Slack)
        @objective(mMP, Max, linobj - nodeinfo.tVal*Slack)
      else
        @objective(mMP, Max, linobj - 0.5*nodeinfo.tVal*sum( ( ctr.soln.α[i] - α[i])^2 + (ctr.soln.β[i] - β[i])^2 + (ctr.soln.γ[i] - γ[i])^2 + (ctr.soln.δ[i] - δ[i])^2 for i in N))
      end
    elseif CTR_PARAM == CP
      @objective(mMP, Max, linobj)
    elseif CTR_PARAM == FEAS
      @variable(mMP, sLev >= 0)
      @constraint(mMP, LVLConstr, linobj + sLev >= nodeinfo.nodeBd)
      @objective(mMP, Min, sLev)
    else
      @objective(mMP, Min, 0)
    end
end

function addHeurConstr(opfdata,mMP,HEUR)
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = 1:nbuses,1:nlines,1:ngens 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
    λF,λT,μF,μT=mMP[:λF],mMP[:λT],mMP[:μF],mMP[:μT]
  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.
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

function solveNodeMP(opfdata,mMP,nodeinfo,trl_bundles,ctr_bundles,agg_bundles,ctr,CTR_PARAM,
			mpsoln)
  init_time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = 1:nbuses,1:nlines,1:ngens 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC

    α,β,γ,δ=mMP[:α],mMP[:β],mMP[:γ],mMP[:δ]
    x,λF,λT,μF,μT=mMP[:x],mMP[:λF],mMP[:λT],mMP[:μF],mMP[:μT]
    ζpLB,ζpUB,ζqLB,ζqUB=mMP[:ζpLB],mMP[:ζpUB],mMP[:ζqLB],mMP[:ζqUB]
    linobj=mMP[:linobj]
    psi = 0
    if CTR_PARAM == PROX
      psi = mMP[:psi]
      @constraint(mMP, PsiLinConstr,  ctr.linobjval - mMP[:linobj] - mMP[:psi] <= 0)
    end

  # Adding the extra cuts
    if length(agg_bundles) > 0
      @constraint(mMP, CutPlanesAgg[n=1:length(agg_bundles)], 0 >= -psi + agg_bundles[n].psival + agg_bundles[n].etahat
	+ sum( agg_bundles[n].eta_sg.α[i]*(α[i]-agg_bundles[n].soln.α[i]) + agg_bundles[n].eta_sg.β[i]*(β[i]-agg_bundles[n].soln.β[i])
	+ agg_bundles[n].eta_sg.γ[i]*(γ[i]-agg_bundles[n].soln.γ[i]) + agg_bundles[n].eta_sg.δ[i]*(δ[i]-agg_bundles[n].soln.δ[i]) for i in N)
        + sum( agg_bundles[n].eta_sg.λF[l]*(λF[l]-agg_bundles[n].soln.λF[l]) + agg_bundles[n].eta_sg.λT[l]*(λT[l]-agg_bundles[n].soln.λT[l]) 
	+ agg_bundles[n].eta_sg.μF[l]*(μF[l]-agg_bundles[n].soln.μF[l]) + agg_bundles[n].eta_sg.μT[l]*(μT[l]-agg_bundles[n].soln.μT[l]) for l in L)
      )
    end
    if length(trl_bundles) > 0
      @constraint(mMP, CutPlanesTrl[n=1:length(trl_bundles)], 0 >= -psi + trl_bundles[n].eta 
	+ sum( trl_bundles[n].eta_sg.α[i]*(α[i]-trl_bundles[n].soln.α[i]) + trl_bundles[n].eta_sg.β[i]*(β[i]-trl_bundles[n].soln.β[i])
	+ trl_bundles[n].eta_sg.γ[i]*(γ[i]-trl_bundles[n].soln.γ[i]) + trl_bundles[n].eta_sg.δ[i]*(δ[i]-trl_bundles[n].soln.δ[i]) for i in N)
        + sum( trl_bundles[n].eta_sg.λF[l]*(λF[l]-trl_bundles[n].soln.λF[l]) + trl_bundles[n].eta_sg.λT[l]*(λT[l]-trl_bundles[n].soln.λT[l]) 
	+ trl_bundles[n].eta_sg.μF[l]*(μF[l]-trl_bundles[n].soln.μF[l]) + trl_bundles[n].eta_sg.μT[l]*(μT[l]-trl_bundles[n].soln.μT[l]) for l in L)
      )
    end
    if length(ctr_bundles) > 0
      @constraint(mMP, CutPlanesCtr[n=1:length(ctr_bundles)], 0 >= -psi + ctr_bundles[n].eta 
	+ sum( ctr_bundles[n].eta_sg.α[i]*(α[i]-ctr_bundles[n].soln.α[i]) + ctr_bundles[n].eta_sg.β[i]*(β[i]-ctr_bundles[n].soln.β[i])
	+ ctr_bundles[n].eta_sg.γ[i]*(γ[i]-ctr_bundles[n].soln.γ[i]) + ctr_bundles[n].eta_sg.δ[i]*(δ[i]-ctr_bundles[n].soln.δ[i]) for i in N)
        + sum( ctr_bundles[n].eta_sg.λF[l]*(λF[l]-ctr_bundles[n].soln.λF[l]) + ctr_bundles[n].eta_sg.λT[l]*(λT[l]-ctr_bundles[n].soln.λT[l]) 
	+ ctr_bundles[n].eta_sg.μF[l]*(μF[l]-ctr_bundles[n].soln.μF[l]) + ctr_bundles[n].eta_sg.μT[l]*(μT[l]-ctr_bundles[n].soln.μT[l]) for l in L)
        )
    end
  mpsoln.init_time += (time_ns()-init_time_Start)/1e9

  ### END DEFINING THE LaGRANGIAN DUAL PROBLEM

  solve_time_Start = time_ns()
  JuMP.optimize!(mMP)
  mpsoln.solve_time += (time_ns()-solve_time_Start)/1e9

  mpsoln.status=JuMP.termination_status(mMP)
  if mpsoln.status == MOI.OPTIMAL || mpsoln.status == MOI.LOCALLY_SOLVED || mpsoln.status == MOI.ALMOST_LOCALLY_SOLVED
      pp_time_Start = time_ns()

      for i in N
          mpsoln.soln.α[i],mpsoln.soln.β[i],mpsoln.soln.γ[i],mpsoln.soln.δ[i] = JuMP.value(α[i]), JuMP.value(β[i]), JuMP.value(γ[i]), JuMP.value(δ[i])
      end
      for g in G
          mpsoln.soln.ζpLB[g],mpsoln.soln.ζpUB[g],mpsoln.soln.ζqLB[g],mpsoln.soln.ζqUB[g] = JuMP.value(ζpLB[g]), JuMP.value(ζpUB[g]), JuMP.value(ζqLB[g]), JuMP.value(ζqUB[g])
      end
      for l in L
        mpsoln.soln.x[l] = JuMP.value(x[l])
        mpsoln.soln.λF[l], mpsoln.soln.λT[l], mpsoln.soln.μF[l], mpsoln.soln.μT[l] = JuMP.value(λF[l]), JuMP.value(λT[l]), JuMP.value(μF[l]), JuMP.value(μT[l])
      end
      mpsoln.objval = JuMP.objective_value(mMP)
      mpsoln.linobjval = JuMP.value(linobj)
      if CTR_PARAM == PROX
        mpsoln.psival = JuMP.value(psi)
      else
        mpsoln.psival = 0
      end
      mpsoln.pp_time += (time_ns()-pp_time_Start)/1e9

      sg_time_Start = time_ns()
      mpsoln.eta,stat = computeSG(opfdata,mpsoln) #This computes mpsoln.eta 
      mpsoln.sg_time += (time_ns()-sg_time_Start)/1e9

      pp_time_Start = time_ns()


      for n=1:length(trl_bundles)
        trl_bundles[n].cut_dual = abs(JuMP.dual(CutPlanesTrl[n]))
      end
      for n=1:length(ctr_bundles)
        ctr_bundles[n].cut_dual = abs(JuMP.dual(CutPlanesCtr[n]))
      end
      for n=1:length(agg_bundles)
        agg_bundles[n].cut_dual = abs(JuMP.dual(CutPlanesAgg[n]))
      end
      psiLCDual=0
      if CTR_PARAM == PROX
        psiLCDual = abs(JuMP.dual(PsiLinConstr))
      end

      update_rho(nodeinfo,trl_bundles,ctr_bundles,agg_bundles)
      if CTR_PARAM == PROX
	#nodeinfo.rho += psiLCDual
      end
      nodeinfo.rhoUB = nodeinfo.rho
      rho_est=node_data.rho
      update_agg(opfdata,node_data,ctr,mpsoln,node_data.sg_agg)
      node_data.agg_sg_norm=comp_norm(opfdata,node_data.sg_agg)
      if CTR_PARAM == PROX0
        node_data.epshat = mpsoln.linobjval - (ctr.linobjval - node_data.rho*ctr.eta) - (1.0/node_data.tVal)*node_data.agg_sg_norm^2
        nodeinfo.descent_est = mpsoln.linobjval-(ctr.linobjval - rho_est*ctr.eta) 
        nodeinfo.descent = (mpsoln.linobjval - rho_est*mpsoln.eta)-(ctr.linobjval - rho_est*ctr.eta) 
      elseif CTR_PARAM == PROX
        node_data.epshat = max(0,ctr.eta)-mpsoln.psival - (1.0/node_data.tVal)*node_data.agg_sg_norm^2
        nodeinfo.descent_est = node_data.epshat + (0.5/node_data.tVal)*node_data.agg_sg_norm^2 
        nodeinfo.descent = max(0,ctr.eta) - max(ctr.linobjval - mpsoln.linobjval, mpsoln.eta) 
      end
      node_data.linerr = nodeinfo.descent  ###TO BE MODIFIED BELOW

     # CONTINUE COMPUTING LINEAR ERRORS
      node_data.linerr -= dot( node_data.sg_agg.α[N], (ctr.soln.α[N]-mpsoln.soln.α[N]) ) 
      node_data.linerr -= dot( node_data.sg_agg.β[N], (ctr.soln.β[N]-mpsoln.soln.β[N]) ) 
      node_data.linerr -= dot( node_data.sg_agg.γ[N], (ctr.soln.γ[N]-mpsoln.soln.γ[N]) ) 
      node_data.linerr -= dot( node_data.sg_agg.δ[N], (ctr.soln.δ[N]-mpsoln.soln.δ[N]) ) 
      node_data.linerr -= dot( node_data.sg_agg.ζpLB[G], (ctr.soln.ζpLB[G]-mpsoln.soln.ζpLB[G]) ) 
      node_data.linerr -= dot( node_data.sg_agg.ζpUB[G], (ctr.soln.ζpUB[G]-mpsoln.soln.ζpUB[G]) ) 
      node_data.linerr -= dot( node_data.sg_agg.ζqLB[G], (ctr.soln.ζqLB[G]-mpsoln.soln.ζqLB[G]) ) 
      node_data.linerr -= dot( node_data.sg_agg.ζqUB[G], (ctr.soln.ζqUB[G]-mpsoln.soln.ζqUB[G]) ) 
      node_data.linerr -= dot( node_data.sg_agg.x[L],(ctr.soln.x[L]-mpsoln.soln.x[L]) ) 
      node_data.linerr -= dot( node_data.sg_agg.λF[L],(ctr.soln.λF[L]-mpsoln.soln.λF[L]) ) 
      node_data.linerr -= dot( node_data.sg_agg.λT[L],(ctr.soln.λT[L]-mpsoln.soln.λT[L]) ) 
      node_data.linerr -= dot( node_data.sg_agg.μF[L],(ctr.soln.μF[L]-mpsoln.soln.μF[L]) ) 
      node_data.linerr -= dot( node_data.sg_agg.μT[L],(ctr.soln.μT[L]-mpsoln.soln.μT[L]) )
      for n=1:length(trl_bundles)
        node_data.linerr += trl_bundles[n].cut_dual*dot( trl_bundles[n].eta_sg.α[N],(ctr.soln.α[N]-mpsoln.soln.α[N])) 
        node_data.linerr += trl_bundles[n].cut_dual*dot( trl_bundles[n].eta_sg.β[N],(ctr.soln.β[N]-mpsoln.soln.β[N]))
        node_data.linerr += trl_bundles[n].cut_dual*dot( trl_bundles[n].eta_sg.γ[N],(ctr.soln.γ[N]-mpsoln.soln.γ[N]))
        node_data.linerr += trl_bundles[n].cut_dual*dot( trl_bundles[n].eta_sg.δ[N],(ctr.soln.δ[N]-mpsoln.soln.δ[N]))
        node_data.linerr += trl_bundles[n].cut_dual*dot( trl_bundles[n].eta_sg.λF[L],(ctr.soln.λF[L]-mpsoln.soln.λF[L])) 
        node_data.linerr += trl_bundles[n].cut_dual*dot( trl_bundles[n].eta_sg.λT[L],(ctr.soln.λT[L]-mpsoln.soln.λT[L])) 
        node_data.linerr += trl_bundles[n].cut_dual*dot( trl_bundles[n].eta_sg.μF[L],(ctr.soln.μF[L]-mpsoln.soln.μF[L])) 
        node_data.linerr += trl_bundles[n].cut_dual*dot( trl_bundles[n].eta_sg.μT[L],(ctr.soln.μT[L]-mpsoln.soln.μT[L])) 
      end
      for n=1:length(ctr_bundles)
        node_data.linerr += ctr_bundles[n].cut_dual*dot( ctr_bundles[n].eta_sg.α[N],(ctr.soln.α[N]-mpsoln.soln.α[N])) 
        node_data.linerr += ctr_bundles[n].cut_dual*dot( ctr_bundles[n].eta_sg.β[N],(ctr.soln.β[N]-mpsoln.soln.β[N]))
        node_data.linerr += ctr_bundles[n].cut_dual*dot( ctr_bundles[n].eta_sg.γ[N],(ctr.soln.γ[N]-mpsoln.soln.γ[N]))
        node_data.linerr += ctr_bundles[n].cut_dual*dot( ctr_bundles[n].eta_sg.δ[N],(ctr.soln.δ[N]-mpsoln.soln.δ[N]))
        node_data.linerr += ctr_bundles[n].cut_dual*dot( ctr_bundles[n].eta_sg.λF[L],(ctr.soln.λF[L]-mpsoln.soln.λF[L])) 
        node_data.linerr += ctr_bundles[n].cut_dual*dot( ctr_bundles[n].eta_sg.λT[L],(ctr.soln.λT[L]-mpsoln.soln.λT[L])) 
        node_data.linerr += ctr_bundles[n].cut_dual*dot( ctr_bundles[n].eta_sg.μF[L],(ctr.soln.μF[L]-mpsoln.soln.μF[L])) 
        node_data.linerr += ctr_bundles[n].cut_dual*dot( ctr_bundles[n].eta_sg.μT[L],(ctr.soln.μT[L]-mpsoln.soln.μT[L])) 
      end
      for n=1:length(agg_bundles)
        node_data.linerr += agg_bundles[n].cut_dual*dot( agg_bundles[n].eta_sg.α[N],(ctr.soln.α[N]-mpsoln.soln.α[N])) 
        node_data.linerr += agg_bundles[n].cut_dual*dot( agg_bundles[n].eta_sg.β[N],(ctr.soln.β[N]-mpsoln.soln.β[N]))
        node_data.linerr += agg_bundles[n].cut_dual*dot( agg_bundles[n].eta_sg.γ[N],(ctr.soln.γ[N]-mpsoln.soln.γ[N]))
        node_data.linerr += agg_bundles[n].cut_dual*dot( agg_bundles[n].eta_sg.δ[N],(ctr.soln.δ[N]-mpsoln.soln.δ[N]))
        node_data.linerr += agg_bundles[n].cut_dual*dot( agg_bundles[n].eta_sg.λF[L],(ctr.soln.λF[L]-mpsoln.soln.λF[L])) 
        node_data.linerr += agg_bundles[n].cut_dual*dot( agg_bundles[n].eta_sg.λT[L],(ctr.soln.λT[L]-mpsoln.soln.λT[L])) 
        node_data.linerr += agg_bundles[n].cut_dual*dot( agg_bundles[n].eta_sg.μF[L],(ctr.soln.μF[L]-mpsoln.soln.μF[L])) 
        node_data.linerr += agg_bundles[n].cut_dual*dot( agg_bundles[n].eta_sg.μT[L],(ctr.soln.μT[L]-mpsoln.soln.μT[L])) 
      end
      if CTR_PARAM == PROX0
        node_data.linerr -= node_data.rho*dot( mpsoln.eta_sg.α[N], (ctr.soln.α[N]-mpsoln.soln.α[N]) ) 
        node_data.linerr -= node_data.rho*dot( mpsoln.eta_sg.β[N], (ctr.soln.β[N]-mpsoln.soln.β[N]) ) 
        node_data.linerr -= node_data.rho*dot( mpsoln.eta_sg.γ[N], (ctr.soln.γ[N]-mpsoln.soln.γ[N]) ) 
        node_data.linerr -= node_data.rho*dot( mpsoln.eta_sg.δ[N], (ctr.soln.δ[N]-mpsoln.soln.δ[N]) ) 
        node_data.linerr -= node_data.rho*dot( mpsoln.eta_sg.λF[L],(ctr.soln.λF[L]-mpsoln.soln.λF[L]) ) 
        node_data.linerr -= node_data.rho*dot( mpsoln.eta_sg.λT[L],(ctr.soln.λT[L]-mpsoln.soln.λT[L]) ) 
        node_data.linerr -= node_data.rho*dot( mpsoln.eta_sg.μF[L],(ctr.soln.μF[L]-mpsoln.soln.μF[L]) ) 
        node_data.linerr -= node_data.rho*dot( mpsoln.eta_sg.μT[L],(ctr.soln.μT[L]-mpsoln.soln.μT[L]) )
      elseif CTR_PARAM == PROX
	node_data.linerr -= psiLCDual*(ctr.linobjval-mpsoln.linobjval)
	if ctr.linobjval - mpsoln.linobjval <= mpsoln.eta
          node_data.linerr -= dot( mpsoln.eta_sg.α[N], (ctr.soln.α[N]-mpsoln.soln.α[N]) ) 
          node_data.linerr -= dot( mpsoln.eta_sg.β[N], (ctr.soln.β[N]-mpsoln.soln.β[N]) ) 
          node_data.linerr -= dot( mpsoln.eta_sg.γ[N], (ctr.soln.γ[N]-mpsoln.soln.γ[N]) ) 
          node_data.linerr -= dot( mpsoln.eta_sg.δ[N], (ctr.soln.δ[N]-mpsoln.soln.δ[N]) ) 
          node_data.linerr -= dot( mpsoln.eta_sg.λF[L],(ctr.soln.λF[L]-mpsoln.soln.λF[L]) ) 
          node_data.linerr -= dot( mpsoln.eta_sg.λT[L],(ctr.soln.λT[L]-mpsoln.soln.λT[L]) ) 
          node_data.linerr -= dot( mpsoln.eta_sg.μF[L],(ctr.soln.μF[L]-mpsoln.soln.μF[L]) ) 
          node_data.linerr -= dot( mpsoln.eta_sg.μT[L],(ctr.soln.μT[L]-mpsoln.soln.μT[L]) )
	else
	  node_data.linerr -= (ctr.linobjval-mpsoln.linobjval)
	end
      end
      if node_data.linerr < -1e-8
        println("FLAGGING: linear error = ",node_data.linerr," is negative.")
      end
      if mpsoln.status == MOI.OPTIMAL && (CTR_PARAM==LVL1 || CTR_PARAM == LVL2 || CTR_PARAM == LVLINF )
	mpsoln.lvl_dual = -getdual(mMP[:LVLConstr])
      end
      mpsoln.pp_time += (time_ns()-pp_time_Start)/1e9
  else
    #println("solveNodeMP: Return status ",mpsoln.status)
  end
  return mpsoln.status
end

