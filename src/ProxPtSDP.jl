#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")

CP,PROX0,PROX,LVL1,LVL2,LVLINF,FEAS=0,1,2,3,4,5,6



function solveNodeProxPt(opfdata,nodeinfo,params,bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,CTR_PARAM,
			mpsoln)
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata


  # The master problem MP
      #mMP = Model(solver=CplexSolver(
#=
      mMP = Model(with_optimizer(CPLEX.Optimizer,
        CPX_PARAM_SCRIND=0,
        CPX_PARAM_TILIM=MAX_TIME,
      # CPX_PARAM_MIPDISPLAY=4,
      # CPX_PARAM_MIPINTERVAL=1,
      # CPX_PARAM_NODELIM=1,
      # CPX_PARAM_HEURFREQ=-1,
        CPX_PARAM_THREADS=1
      #CPXPARAM_LPMethod=4  ###BARRIER ALG
      #CPX_PARAM_ADVIND=0
	)
      )
=#
    mMP = Model(with_optimizer(Ipopt.Optimizer))
    #mMP = Model(with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=4))
	### Got error indicating lack of support for quadratic objective ???

    # each x[l] is either 0 (line l active) or 1 (line l is cut)
      @variable(mMP, nodeinfo.x_lbs[l] <= x[l=L] <= nodeinfo.x_ubs[l])
    # dual multipliers associated with active power flow balance constraints
      @variable(mMP, -1 <= α[i=N] <= 1, start=0)
    # dual multipliers associated with active power flow balance constraints
      @variable(mMP, -1 <= β[i=N] <= 1, start=0)
    # dual multipliers associated with the voltage magnitude bounds
      @variable(mMP, δ[i=N] >= 0); @variable(mMP, γ[i=N] >= 0); @constraint(mMP, [i=N], δ[i]+γ[i] <= 1)
    # dual multipliers associated with the power generation bounds
      @variable(mMP, 0 <= ζpUB[g=G] <= 1); @variable(mMP, 0 <= ζpLB[g=G] <= 1); @variable(mMP, 0 <= ζqUB[g=G] <= 1); @variable(mMP, 0 <= ζqLB[g=G] <= 1)

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
      @variable(mMP, λF[l=L], start=0); @variable(mMP, λT[l=L], start=0); @variable(mMP, μF[l=L], start=0); @variable(mMP, μT[l=L], start=0)
    @constraint(mMP, AMcf1[l in L], α[fromBus[l]] - x[l] <= λF[l]) 
    @constraint(mMP, AMcf2[l in L], α[fromBus[l]] + x[l] >= λF[l])
    @constraint(mMP, AMcf3[l in L], -(1 - x[l]) <= λF[l]) 
    @constraint(mMP, AMcf4[l in L], (1 - x[l]) >= λF[l])
    @constraint(mMP, AMct1[l in L], α[toBus[l]] - x[l] <= λT[l]) 
    @constraint(mMP, AMct2[l in L], α[toBus[l]] + x[l] >= λT[l])
    @constraint(mMP, AMct3[l in L], -(1 - x[l]) <= λT[l])
    @constraint(mMP, AMct4[l in L], (1 - x[l]) >= λT[l])

    @constraint(mMP, BMcf1[l in L], β[fromBus[l]] - x[l] <= μF[l])
    @constraint(mMP, BMcf2[l in L], β[fromBus[l]] + x[l] >= μF[l])
    @constraint(mMP, BMcf3[l in L], -(1 - x[l]) <= μF[l])
    @constraint(mMP, BMcf4[l in L], (1 - x[l]) >= μF[l])
    @constraint(mMP, BMct1[l in L], β[toBus[l]] - x[l] <= μT[l])
    @constraint(mMP, BMct2[l in L], β[toBus[l]] + x[l] >= μT[l])
    @constraint(mMP, BMct3[l in L], -(1 - x[l]) <= μT[l])
    @constraint(mMP, BMct4[l in L], (1 - x[l]) >= μT[l])

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
    
    @variable(mMP, psi)
    if CTR_PARAM != PROX
      setlowerbound(psi,0)
      setupperbound(psi,0)
    end
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
   # METHOD-SPECIFIC SETUP
    if CTR_PARAM == PROX0
      @objective(mMP, Max, linobj - 0.5*params.tVal*(
	  sum( (ctr.soln.α[i] - α[i])^2 + (ctr.soln.β[i] - β[i])^2 + (ctr.soln.γ[i] - γ[i])^2 + (ctr.soln.δ[i] - δ[i])^2 for i in N)
	  +sum( (ctr.soln.ζpLB[g] - ζpLB[g])^2 + (ctr.soln.ζpUB[g] - ζpUB[g])^2 + (ctr.soln.ζqLB[g] - ζqLB[g])^2 + (ctr.soln.ζqUB[g] - ζqUB[g])^2 for g in G)
	  +sum( (ctr.soln.x[l]-x[l])^2 + (ctr.soln.λF[l] - λF[l])^2 + (ctr.soln.λT[l] - λT[l])^2 + (ctr.soln.μF[l] - μF[l])^2 + (ctr.soln.μT[l] - μT[l])^2 for l in L)
	)
      )
    elseif CTR_PARAM == PROX 
        @objective(mMP, Min, psi + 0.5*params.tVal*(
	  sum( (ctr.soln.α[i] - α[i])^2 + (ctr.soln.β[i] - β[i])^2 + (ctr.soln.γ[i] - γ[i])^2 + (ctr.soln.δ[i] - δ[i])^2 for i in N)
	  +sum( (ctr.soln.ζpLB[g] - ζpLB[g])^2 + (ctr.soln.ζpUB[g] - ζpUB[g])^2 + (ctr.soln.ζqLB[g] - ζqLB[g])^2 + (ctr.soln.ζqUB[g] - ζqUB[g])^2 for g in G)
	  +sum( (ctr.soln.x[l]-x[l])^2 + (ctr.soln.λF[l] - λF[l])^2 + (ctr.soln.λT[l] - λT[l])^2 + (ctr.soln.μF[l] - μF[l])^2 + (ctr.soln.μT[l] - μT[l])^2 for l in L)
	  )
	)
	@constraint(mMP, ctr.linobjval - linobj - psi <= 0)
	if length(bundles) > 0
          @constraint(mMP, CutPlanes[n=1:length(bundles)],bundles[n].eta - sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i]) 
	    + bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
            - sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
	    + bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L) - psi <= 0.0
          )
	end
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
        @objective(mMP, Max, linobj - params.tVal*sum(alphaSlack[i] + betaSlack[i] + gammaSlack[i] + deltaSlack[i] for i in N))
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
        @objective(mMP, Max, linobj - params.tVal*Slack)
      else
        @objective(mMP, Max, linobj - 0.5*params.tVal*sum( ( ctr.soln.α[i] - α[i])^2 + (ctr.soln.β[i] - β[i])^2 + (ctr.soln.γ[i] - γ[i])^2 + (ctr.soln.δ[i] - δ[i])^2 for i in N))
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

  # Adding the extra cuts
    if length(agg_bundles) > 0
      @constraint(mMP, CutPlanesAgg[n=1:length(agg_bundles)], 0 <= psi - agg_bundles[n].etahat
	+ sum( agg_bundles[n].eta_sg.α[i]*(α[i]-agg_bundles[n].soln.α[i]) + agg_bundles[n].eta_sg.β[i]*(β[i]-agg_bundles[n].soln.β[i])
	+ agg_bundles[n].eta_sg.γ[i]*(γ[i]-agg_bundles[n].soln.γ[i]) + agg_bundles[n].eta_sg.δ[i]*(δ[i]-agg_bundles[n].soln.δ[i]) for i in N)
        + sum( agg_bundles[n].eta_sg.λF[l]*(λF[l]-agg_bundles[n].soln.λF[l]) + agg_bundles[n].eta_sg.λT[l]*(λT[l]-agg_bundles[n].soln.λT[l]) 
	+ agg_bundles[n].eta_sg.μF[l]*(μF[l]-agg_bundles[n].soln.μF[l]) + agg_bundles[n].eta_sg.μT[l]*(μT[l]-agg_bundles[n].soln.μT[l]) for l in L)
      )
    end
    if length(bundles) > 0
      @constraint(mMP, CutPlanes[n=1:length(bundles)], 0 <= psi - bundles[n].eta 
	+ sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i])
	+ bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
        + sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
	+ bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L)
      )
    end
    if length(ctr_bundles) > 0
      @constraint(mMP, CutPlanesCtr[n=1:length(ctr_bundles)], 0 <= psi - ctr_bundles[n].eta 
	+ sum( ctr_bundles[n].eta_sg.α[i]*(α[i]-ctr_bundles[n].soln.α[i]) + ctr_bundles[n].eta_sg.β[i]*(β[i]-ctr_bundles[n].soln.β[i])
	+ ctr_bundles[n].eta_sg.γ[i]*(γ[i]-ctr_bundles[n].soln.γ[i]) + ctr_bundles[n].eta_sg.δ[i]*(δ[i]-ctr_bundles[n].soln.δ[i]) for i in N)
        + sum( ctr_bundles[n].eta_sg.λF[l]*(λF[l]-ctr_bundles[n].soln.λF[l]) + ctr_bundles[n].eta_sg.λT[l]*(λT[l]-ctr_bundles[n].soln.λT[l]) 
	+ ctr_bundles[n].eta_sg.μF[l]*(μF[l]-ctr_bundles[n].soln.μF[l]) + ctr_bundles[n].eta_sg.μT[l]*(μT[l]-ctr_bundles[n].soln.μT[l]) for l in L)
        )
    end

  ### END DEFINING THE LaGRANGIAN DUAL PROBLEM

  JuMP.optimize!(mMP)
  mpsoln.status=JuMP.termination_status(mMP)
  if mpsoln.status == MOI.OPTIMAL || mpsoln.status == MOI.LOCALLY_SOLVED
      for i in N
          mpsoln.soln.α[i],mpsoln.soln.β[i],mpsoln.soln.γ[i],mpsoln.soln.δ[i] = getvalue(α[i]), getvalue(β[i]), getvalue(γ[i]), getvalue(δ[i])
      end
      for g in G
          mpsoln.soln.ζpLB[g],mpsoln.soln.ζpUB[g],mpsoln.soln.ζqLB[g],mpsoln.soln.ζqUB[g] = getvalue(ζpLB[g]), getvalue(ζpUB[g]), getvalue(ζqLB[g]), getvalue(ζqUB[g])
      end
      for l in L
        mpsoln.soln.x[l] = getvalue(x[l])
        mpsoln.soln.λF[l], mpsoln.soln.λT[l], mpsoln.soln.μF[l], mpsoln.soln.μT[l] = getvalue(λF[l]), getvalue(λT[l]), getvalue(μF[l]), getvalue(μT[l])
      end
      #mpsoln.objval,mpsoln.solvetime = getobjectivevalue(mMP), getsolvetime(mMP)
      mpsoln.objval = getobjectivevalue(mMP)
      mpsoln.linobjval = getvalue(linobj)
      if CTR_PARAM == PROX
        mpsoln.psival = getvalue(psi)
      end

      computeSG(opfdata,mpsoln) #This computes mpsoln.eta 

      mpsoln.linerr = ctr.eta - mpsoln.eta 
      mpsoln.linerr += dot( mpsoln.eta_sg.α[N], (ctr.soln.α[N]-mpsoln.soln.α[N]) ) 
      mpsoln.linerr += dot( mpsoln.eta_sg.β[N], (ctr.soln.β[N]-mpsoln.soln.β[N]) ) 
      mpsoln.linerr += dot( mpsoln.eta_sg.γ[N], (ctr.soln.γ[N]-mpsoln.soln.γ[N]) ) 
      mpsoln.linerr += dot( mpsoln.eta_sg.δ[N], (ctr.soln.δ[N]-mpsoln.soln.δ[N]) ) 
      mpsoln.linerr += dot( mpsoln.eta_sg.λF[L],(ctr.soln.λF[L]-mpsoln.soln.λF[L]) ) 
      mpsoln.linerr += dot( mpsoln.eta_sg.λT[L],(ctr.soln.λT[L]-mpsoln.soln.λT[L]) ) 
      mpsoln.linerr += dot( mpsoln.eta_sg.μF[L],(ctr.soln.μF[L]-mpsoln.soln.μF[L]) ) 
      mpsoln.linerr += dot( mpsoln.eta_sg.μT[L],(ctr.soln.μT[L]-mpsoln.soln.μT[L]) )

      for n=1:length(bundles)
        etaval=-getvalue(-bundles[n].eta 
	  + sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i])
	    + bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
          + sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
	    + bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L)
        )
	if mpsoln.etahat < etaval
	  mpsoln.etahat = etaval
	end
      end
      for n=1:length(ctr_bundles)
        etaval=-getvalue(-ctr_bundles[n].eta 
	  + sum( ctr_bundles[n].eta_sg.α[i]*(α[i]-ctr_bundles[n].soln.α[i]) + ctr_bundles[n].eta_sg.β[i]*(β[i]-ctr_bundles[n].soln.β[i])
	    + ctr_bundles[n].eta_sg.γ[i]*(γ[i]-ctr_bundles[n].soln.γ[i]) + ctr_bundles[n].eta_sg.δ[i]*(δ[i]-ctr_bundles[n].soln.δ[i]) for i in N)
          + sum( ctr_bundles[n].eta_sg.λF[l]*(λF[l]-ctr_bundles[n].soln.λF[l]) + ctr_bundles[n].eta_sg.λT[l]*(λT[l]-ctr_bundles[n].soln.λT[l]) 
	    + ctr_bundles[n].eta_sg.μF[l]*(μF[l]-ctr_bundles[n].soln.μF[l]) + ctr_bundles[n].eta_sg.μT[l]*(μT[l]-ctr_bundles[n].soln.μT[l]) for l in L)
        )
	if mpsoln.etahat < etaval
	  mpsoln.etahat = etaval
	end
      end
      for n=1:length(agg_bundles)
        etaval=-getvalue(-agg_bundles[n].etahat 
	  + sum( agg_bundles[n].eta_sg.α[i]*(α[i]-agg_bundles[n].soln.α[i]) + agg_bundles[n].eta_sg.β[i]*(β[i]-agg_bundles[n].soln.β[i])
	    + agg_bundles[n].eta_sg.γ[i]*(γ[i]-agg_bundles[n].soln.γ[i]) + agg_bundles[n].eta_sg.δ[i]*(δ[i]-agg_bundles[n].soln.δ[i]) for i in N)
          + sum( agg_bundles[n].eta_sg.λF[l]*(λF[l]-agg_bundles[n].soln.λF[l]) + agg_bundles[n].eta_sg.λT[l]*(λT[l]-agg_bundles[n].soln.λT[l]) 
	    + agg_bundles[n].eta_sg.μF[l]*(μF[l]-agg_bundles[n].soln.μF[l]) + agg_bundles[n].eta_sg.μT[l]*(μT[l]-agg_bundles[n].soln.μT[l]) for l in L)
        )
	if mpsoln.etahat < etaval
	  mpsoln.etahat = etaval
	end
      end

      if mpsoln.status == MOI.OPTIMAL || mpsoln.status == MOI.LOCALLY_SOLVED
	for n=1:length(bundles)
	  bundles[n].cut_dual = abs(getdual(CutPlanes[n]))
	end
	for n=1:length(ctr_bundles)
	  ctr_bundles[n].cut_dual = abs(getdual(CutPlanesCtr[n]))
	end
	for n=1:length(agg_bundles)
	  agg_bundles[n].cut_dual = abs(getdual(CutPlanesAgg[n]))
	end
      end
      if mpsoln.status == MOI.OPTIMAL && (CTR_PARAM==LVL1 || CTR_PARAM == LVL2 || CTR_PARAM == LVLINF )
	mpsoln.lvl_dual = -getdual(LVLConstr)
      end
  else
    #println("solveNodeProxPt: Return status ",mpsoln.status)
  end
  return mpsoln.status
end

function testProxTraj(opfdata,params,K,HEUR)
    TOL = 1e-5
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    bundles=Dict()
    ncuts=0
    params.rho = 0
    mpsoln=create_bundle(opfdata)
    ctr=create_bundle(opfdata)
    agg_bundles=Dict()
    agg_bundle=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)


    v_est,ssc_cntr,agg_norm = 1e20,0,0.0
  # MAIN LOOP
    for kk=1:params.maxNSG
     # PHASE I
      tMax = max(tMin,0.5*(comp_norm(opfdata,ctr.eta_sg)^2)/(1+abs(ctr.eta)))
      println("Penalty will be in [",tMin,",",tMax,"].")
      params.tVal = (tMax+tMin)/2
     # PHASE II
      while true
     # STEP 1
        mpsoln=create_bundle(opfdata)
        while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
	  params.tVal /= 2
          println("Resolving with reduced prox parameter value: ",params.tVal)
          status = solveNodeProxPt(opfdata,node_data,params,bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
        end	
        if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
	  ncuts=length(bundles)
	  naggcuts=length(bundles)
	  if ncuts > 0
	    params.rho = sum(bundles[n].cut_dual for n in 1:ncuts) + sum(agg_bundles[n].cut_dual for n in 1:naggcuts) 
           # COMPUTING LINEARIZATION ERRORS
	    agg_norm,epshat,params.rho=update_agg(opfdata,bundles,ctr,mpsoln,sg_agg,ctr_bundles,agg_bundles,agg_bundle)
            if params.rho*mpsoln.eta <= params.ssc*(1/params.tVal)*agg_norm^2 
	      break
	    end
	  end
          agg_bundles[1]=agg_bundle
	  ncuts=purgeSG(opfdata,bundles)
          bundles[ncuts+1]=mpsoln
        else
	  println("Solver returned: $status")
        end
      end
      updateCenter(opfdata,mpsoln,ctr,bundles)
      @show "Traj",kk,mpsoln.linobjval,mpsoln.eta,ncuts,agg_norm,params
      if agg_norm < 1e-2
        println("Convergence to within tolerance: obj, feas ",mpsoln.linobjval," ",mpsoln.eta)
	break
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

function testProxPt0(opfdata,params,K,HEUR,node_data)
    println("Applying the algorithm of Delfino and de Oliveira 2018....")
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    trl_bundles=Dict()
    ctr_bundles=Dict()
    agg_bundles=Dict()
    ncuts=0
    params.rho,params.rhoUB = 0,0
    mpsoln=create_bundle(opfdata)
    ctr=create_bundle(opfdata)
    status = solveNodeProxPt(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
    while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
      params.tVal /= 2
      println("Resolving with reduced prox parameter value: ",params.tVal)
      status = solveNodeProxPt(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
    end	
    updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
    ctr_bundles[1]=mpsoln

    agg_bundle=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)

    v_est,ssc_cntr = 1e20,0
    tL,tU=params.tMin,params.tMax
    TOL = 1e-5
  # MAIN LOOP
    for kk=1:params.maxNSG
     # STEP 1
      mpsoln=create_bundle(opfdata)
      status = solveNodeProxPt(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
      while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
	params.tVal /= 2
        println("Resolving with reduced prox parameter value: ",params.tVal)
        status = solveNodeProxPt(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
      end	
      if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
       # STEP 2
        ncuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
	 # UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
         # COMPUTING AGGREGATION INFORMATION
          agg_norm = update_agg(opfdata,params,ctr,mpsoln,sg_agg)
          update_rho(params,trl_bundles,ctr_bundles,agg_bundles)
          epshat = compute_epshat(opfdata,params,mpsoln,ctr,sg_agg)
	  if params.rhoUB < params.rho
	    params.rhoUB = params.rho + 1
	  end
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
          if ctr.eta < TOL && agg_norm < 1e-3 && epshat < 1e-3 
	    println("Convergence to within tolerance: ")
	    @show kk,ncuts,ssc_cntr,params.tVal,params.rho
	    @show mpsoln.linobjval,mpsoln.eta,agg_norm,epshat,mpsoln.linerr
	    break
          end
         # STEP 3
	  sscval = ((mpsoln.linobjval - params.rhoUB*mpsoln.eta)-(ctr.linobjval - params.rhoUB*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - params.rhoUB*ctr.eta)) 
	  vval = (mpsoln.linobjval-(ctr.linobjval - params.rhoUB*ctr.eta)) 
          #params.tVal,v_est,ssc_cntr=KiwielRhoUpdate(opfdata,params,ctr,mpsoln,sscval,vval,agg_norm,epshat,agg_bundles[1],v_est,ssc_cntr)
          if sscval >= params.ssc 
         # UPDATE CENTER VALUES
	    if testSchrammZoweSSII(opfdata,params,ctr,mpsoln)
	      nctrcuts=purgeSG(opfdata,ctr_bundles)
	      ctr_bundles[nctrcuts+1]=mpsoln
	      updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
    	      tL,tU=params.tMin,params.tMax
	      @show kk,ncuts,ssc_cntr,params.tVal,params.rho
	      @show mpsoln.linobjval,mpsoln.eta,agg_norm,epshat
	    else
	      tU = params.tVal    
	      params.tVal = (tU+tL)/2.0
	      #@show params.tVal
	    end
	  else
	    if testSchrammZoweNSII(opfdata,params,ctr,mpsoln,agg_bundles)
	      ncuts=purgeSG(opfdata,trl_bundles)
              trl_bundles[ncuts+1]=mpsoln
    	      tL,tU=params.tMin,params.tMax
	    else
	      tL = params.tVal    
	      params.tVal = (tU+tL)/2.0
	      #@show params.tVal
	    end
          end
      else
	println("Solver returned: $status")
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

# Implements the approach of Sagastizabal and Solodov 2005
function testProxPt(opfdata,params,K,HEUR,node_data)
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    bundles=Dict()
    mpsoln=create_bundle(opfdata)
    mpsoln.soln.x[L] = x_val[L]
    sg_agg=create_soln(opfdata)
    ctr=mpsoln

    v_est,ssc_cntr = 1e20,0
  # MAIN LOOP
    for kk=1:params.maxNSG
     # STEP 1
      mpsoln=create_bundle(opfdata)
      status = solveNodeProxPt(opfdata,node_data,params,bundles,K,HEUR,ctr,PROX,mpsoln)
      while status != :Optimal
	params.tVal /= 2
        println("Resolving with reduced prox parameter value: ",params.tVal)
        status = solveNodeProxPt(opfdata,node_data,params,bundles,K,HEUR,ctr,PROX,mpsoln)
      end	
      if status == :Optimal 
	comp_agg(opfdata,ctr.soln,mpsoln.soln,sg_agg)
	agg_norm=comp_norm(opfdata,sg_agg)
	epshat=max(0,ctr.eta)-mpsoln.psival-(1.0/params.tVal)*agg_norm^2
	#del=epshat+0.5*(1.0/params.tVal)*agg_norm^2
        del=ctr.eta-mpsoln.objval
       # STEP 2
        if del < 1e-6 
	  println("Convergence to within tolerance: del = ",del," val: ",mpsoln.linobjval," and feas: ",mpsoln.eta)
	  break
        end
       # STEP 3
	hk = max(ctr.linobjval-mpsoln.linobjval,mpsoln.eta)
	sscval = max(0,ctr.eta) - hk
        if hk <= max(0,ctr.eta) - params.ssc*del || kk==1
          # UPDATE CENTER VALUES
	    updateCenter(opfdata,mpsoln,ctr,bundles)
	    #purgeSG(opfdata,bundles)
	    oldncuts = length(bundles)
	    ncuts = oldncuts
	    for n=oldncuts:-1:1
	      if ctr.linobjval-bundles[n].linobjval > bundles[n].psival
		bundles[n]=bundles[ncuts]
	        delete!(bundles,ncuts)
		ncuts -= 1
	      end
	    end
	    if ctr.linobjval-mpsoln.linobjval > mpsoln.psival
	      ncuts = length(bundles)
              bundles[ncuts+1]=mpsoln
	    end
            @show "Sagadov",kk,ctr.linobjval,ctr.eta,epshat,del,hk,length(bundles),params
	else
	  ncuts = length(bundles)
          bundles[ncuts+1]=mpsoln
        end
       # STEP 4
        #tVal,v_est,ssc_cntr = KiwielRhoUpdate(opfdata,params,mpsoln,sscval,del,agg_norm,epshat,v_est,ssc_cntr)
      else
	#println("Solve status: $status")
	#break
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

#Implements the approach of De Oliveira 2016
function testLevelBM(opfdata,params,K,HEUR)
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    optUB=1e20
    # INITIAL SOLN (ALL ZEROS)
      mpsoln=create_bundle(opfdata)
    bestsoln = mpsoln
    ctr=mpsoln 
    bundles = Dict()
    prox_term=LVLINF

    cpsoln=create_bundle(opfdata)
    solveNodeProxPt(opfdata,node_data,params,bundles,K,HEUR,ctr,CP,cpsoln)
    if cpsoln.status == :Optimal 
	optUB = cpsoln.objval
        hkctr = max(optUB,0)
    end

  # MAIN LOOP
    hkctr_updated = false
    optub_updated = false
    for kk=1:params.maxNSG
     # STEP 1
      hkval = max(optUB - ctr.linobjval,ctr.eta)
      bestsoln = ctr
      bestIdx = 0
      ncuts=length(bundles)
      for pp=ncuts:-1:1
	if max(optUB - bundles[pp].linobjval,bundles[pp].eta) < hkval
	  hkval = max(optUB - bundles[pp].linobjval,bundles[pp].eta) 
          bestsoln = bundles[pp]
	  bestIdx=pp
	end
      end
      node_data.nodeBd = optUB - params.ssc*hkval
      if hkval < 1e-6
	println("Convergence to within tolerance.")
	break
      end

     # STEP 2
        if hkval <= (1-params.ssc)*hkctr 
          hkctr = hkval
          hkctr_updated = true
          # UPDATE CENTER VALUES
	    if bestIdx > 0
	      updateCenter(opfdata,bestsoln,ctr,bundles)
	    end
	  @show "hk  update",optUB,ctr.linobjval,ctr.eta,hkval,ncuts
        end

     # STEP 3
      cpstatus = solveNodeProxPt(opfdata,node_data,params,bundles,K,HEUR,ctr,CP,cpsoln)
      if cpsoln.status != :Optimal
	println("FLAG! cpsoln.status = ",cpsoln.status)
      end
      if node_data.nodeBd-cpsoln.linobjval <= 0.0
	params.tVal = 0.1
        mpsoln=create_bundle(opfdata)
        status=solveNodeProxPt(opfdata,node_data,params,bundles,K,HEUR,ctr,prox_term,mpsoln)
        while status != :Optimal
	  params.tVal /= 2
          println("Resolving with reduced prox parameter value: ",params.tVal)
          status=solveNodeProxPt(opfdata,node_data,params,bundles,K,HEUR,ctr,prox_term,mpsoln)
        end	
	if mpsoln.status != :Optimal
	  println("Taking recourse since mpsoln does not have an optimal solution for MP. Feas: ",node_data.nodeBd - mpsoln.linobjval," val: ",node_data.nodeBd - cpsoln.linobjval)
	  cpy_soln(opfdata,cpsoln,mpsoln)
	else
	  if hkctr_updated
	    oldncuts = ncuts
	    ncuts=purgeSG(opfdata,bundles)
	    #println("There were ",oldncuts," cuts; after purging, there are ",ncuts)
	    hkctr_updated = false
	  end
        end

        if mpsoln.eta > 0
          bundles[ncuts+1] = mpsoln
	  ncuts += 1
        else
            println("Tolerance met for not generating a new lazy cut, obj,eta=",mpsoln.linobjval,mpsoln.eta,".")
	    if mpsoln.linobjval < optUB
	      optUB = mpsoln.linobjval
              optub_updated = true
	    end
	    @show "bnd opt update",optUB,mpsoln.linobjval,-mpsoln.eta,hkval,ncuts
	end
      else
	  optUB = node_data.nodeBd 
          optub_updated = true
          hkctr = max(optUB - bestsoln.linobjval,bestsoln.eta)
          hkctr_updated = true
	  if bestIdx > 0
	    updateCenter(opfdata,bestsoln,ctr,bundles)
	  end
	  @show "bnd update",kk,optUB,bestsoln.linobjval,bestsoln.eta,hkval,ncuts
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end


