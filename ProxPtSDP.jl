#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")

maxNSG = 100000
CP,PROX0,PROX,LVL1,LVL2,LVLINF,FEAS=0,1,2,3,4,5,6
tMin = 0.01
tMax = 20
tVal = 0.5
ssc=0.5

mutable struct NodeInfo
  x_lbs::Array{Float64}
  x_ubs::Array{Float64}
  nodeBd::Float64
end
mutable struct Soln
  α::Array{Float64}
  β::Array{Float64}
  γ::Array{Float64}
  δ::Array{Float64}
  ζpLB::Array{Float64}
  ζpUB::Array{Float64}
  ζqLB::Array{Float64}
  ζqUB::Array{Float64}
  x::Array{Float64}
  λF::Array{Float64}
  λT::Array{Float64}
  μF::Array{Float64}
  μT::Array{Float64}
end
function create_soln(opfdata)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  return Soln(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
    zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
    zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines))
end
function cpy_soln(opfdata,fromSoln,toSoln)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  toSoln.α[N] = fromSoln.α[N]
  toSoln.β[N] = fromSoln.β[N]
  toSoln.γ[N] = fromSoln.γ[N]
  toSoln.δ[N] = fromSoln.δ[N]
  toSoln.ζpLB[G] = fromSoln.ζpLB[G]
  toSoln.ζpUB[G] = fromSoln.ζpUB[G]
  toSoln.ζqLB[G] = fromSoln.ζqLB[G]
  toSoln.ζqUB[G] = fromSoln.ζqUB[G]
  toSoln.x[L] = fromSoln.x[L]
  toSoln.λF[L] = fromSoln.λF[L]
  toSoln.λT[L] = fromSoln.λT[L]
  toSoln.μF[L] = fromSoln.μF[L]
  toSoln.μT[L] = fromSoln.μT[L]
end
function comp_agg(opfdata,ctr,trl,agg)
  global tVal
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  agg.α[N] = tVal*(ctr.α[N]-trl.α[N])
  agg.β[N] = tVal*(ctr.β[N]-trl.β[N])
  agg.γ[N] = tVal*(ctr.γ[N]-trl.γ[N])
  agg.δ[N] = tVal*(ctr.δ[N]-trl.δ[N])
  agg.ζpLB[G] = tVal*(ctr.ζpLB[G]-trl.ζpLB[G])
  agg.ζpUB[G] = tVal*(ctr.ζpUB[G]-trl.ζpUB[G])
  agg.ζqLB[G] = tVal*(ctr.ζqLB[G]-trl.ζqLB[G])
  agg.ζqUB[G] = tVal*(ctr.ζqUB[G]-trl.ζqUB[G])
  agg.x[L] = tVal*(ctr.x[L]-trl.x[L])
  agg.λF[L] = tVal*(ctr.λF[L]-trl.λF[L])
  agg.λT[L] = tVal*(ctr.λT[L]-trl.λT[L])
  agg.μF[L] = tVal*(ctr.μF[L]-trl.μF[L])
  agg.μT[L] = tVal*(ctr.μT[L]-trl.μT[L])
end
function comp_norm(opfdata,soln)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  return sqrt(norm(soln.α[N])^2 + norm(soln.β[N])^2 + norm(soln.γ[N])^2 + norm(soln.δ[N])^2 
	+ norm(soln.ζpLB[G])^2 + norm(soln.ζpUB[G])^2 + norm(soln.ζqLB[G])^2 + norm(soln.ζqUB[G])^2
  	+ norm(soln.x[L])^2 
	+ norm(soln.λF[L])^2 + norm(soln.λT[L])^2 + norm(soln.μF[L])^2 + norm(soln.μT[L])^2)
end

mutable struct Bundle
  soln::Soln
  eta_sg::Soln
  objval::Float64
  linobjval::Float64
  penval::Float64
  psival::Float64
  eta::Float64
  etahat::Float64
  linerr::Float64
  cut_dual::Float64
  lvl_dual::Float64
  age::Float64
  solvetime::Float64
  status
end
function create_bundle(opfdata)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  return Bundle(create_soln(opfdata),create_soln(opfdata),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0)
end
function cpy_bundle(opfdata,fromBundle,toBundle)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  cpy_soln(opfdata,fromBundle.soln,toBundle.soln)
  cpy_soln(opfdata,fromBundle.eta_sg,toBundle.eta_sg)
  toBundle.objval = fromBundle.objval
  toBundle.linobjval = fromBundle.linobjval
  toBundle.penval = fromBundle.penval
  toBundle.psival = fromBundle.psival
  toBundle.eta = fromBundle.eta
  toBundle.etahat = fromBundle.etahat
  toBundle.linerr = fromBundle.linerr
  toBundle.cut_dual = fromBundle.cut_dual
  toBundle.lvl_dual = fromBundle.lvl_dual
  toBundle.solvetime = fromBundle.solvetime
  toBundle.status = fromBundle.status
end

mutable struct ConstrDuals
  var_bds::Soln
  zeta_eq::Array{Float64,2}
  budget::Float64
  AMcf1::Array{Float64} 
  AMcf2::Array{Float64} 
  AMcf3::Array{Float64} 
  AMcf4::Array{Float64} 
  AMct1::Array{Float64} 
  AMct2::Array{Float64} 
  AMct3::Array{Float64} 
  AMct4::Array{Float64} 
  BMcf1::Array{Float64} 
  BMcf2::Array{Float64} 
  BMcf3::Array{Float64} 
  BMcf4::Array{Float64} 
  BMct1::Array{Float64} 
  BMct2::Array{Float64} 
  BMct3::Array{Float64} 
  BMct4::Array{Float64} 
end
function create_constr_duals(opfdata)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  return ConstrDuals(create_soln(opfdata),zeros(nbuses,ngens),0,
    zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),
    zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines))
end

function solveNodeProxPt(opfdata,nodeinfo,bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,CTR_PARAM,
			mpsoln)
  global tVal
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
      @objective(mMP, Max, linobj - 0.5*tVal*(
	  sum( (ctr.soln.α[i] - α[i])^2 + (ctr.soln.β[i] - β[i])^2 + (ctr.soln.γ[i] - γ[i])^2 + (ctr.soln.δ[i] - δ[i])^2 for i in N)
	  +sum( (ctr.soln.ζpLB[g] - ζpLB[g])^2 + (ctr.soln.ζpUB[g] - ζpUB[g])^2 + (ctr.soln.ζqLB[g] - ζqLB[g])^2 + (ctr.soln.ζqUB[g] - ζqUB[g])^2 for g in G)
	  +sum( (ctr.soln.x[l]-x[l])^2 + (ctr.soln.λF[l] - λF[l])^2 + (ctr.soln.λT[l] - λT[l])^2 + (ctr.soln.μF[l] - μF[l])^2 + (ctr.soln.μT[l] - μT[l])^2 for l in L)
	)
      )
      if length(agg_bundles) > 0
        @constraint(mMP, CutPlanesAgg[n=1:length(agg_bundles)], 0 <= -agg_bundles[n].etahat
	  + sum( agg_bundles[n].eta_sg.α[i]*(α[i]-agg_bundles[n].soln.α[i]) + agg_bundles[n].eta_sg.β[i]*(β[i]-agg_bundles[n].soln.β[i])
	    + agg_bundles[n].eta_sg.γ[i]*(γ[i]-agg_bundles[n].soln.γ[i]) + agg_bundles[n].eta_sg.δ[i]*(δ[i]-agg_bundles[n].soln.δ[i]) for i in N)
          + sum( agg_bundles[n].eta_sg.λF[l]*(λF[l]-agg_bundles[n].soln.λF[l]) + agg_bundles[n].eta_sg.λT[l]*(λT[l]-agg_bundles[n].soln.λT[l]) 
	    + agg_bundles[n].eta_sg.μF[l]*(μF[l]-agg_bundles[n].soln.μF[l]) + agg_bundles[n].eta_sg.μT[l]*(μT[l]-agg_bundles[n].soln.μT[l]) for l in L)
        )
      end
      if length(bundles) > 0
        @constraint(mMP, CutPlanes[n=1:length(bundles)], 0 <= -bundles[n].eta 
	  + sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i])
	    + bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
          + sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
	    + bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L)
        )
      end
      if length(ctr_bundles) > 0
        @constraint(mMP, CutPlanesCtr[n=1:length(ctr_bundles)], 0 <= -ctr_bundles[n].eta 
	  + sum( ctr_bundles[n].eta_sg.α[i]*(α[i]-ctr_bundles[n].soln.α[i]) + ctr_bundles[n].eta_sg.β[i]*(β[i]-ctr_bundles[n].soln.β[i])
	    + ctr_bundles[n].eta_sg.γ[i]*(γ[i]-ctr_bundles[n].soln.γ[i]) + ctr_bundles[n].eta_sg.δ[i]*(δ[i]-ctr_bundles[n].soln.δ[i]) for i in N)
          + sum( ctr_bundles[n].eta_sg.λF[l]*(λF[l]-ctr_bundles[n].soln.λF[l]) + ctr_bundles[n].eta_sg.λT[l]*(λT[l]-ctr_bundles[n].soln.λT[l]) 
	    + ctr_bundles[n].eta_sg.μF[l]*(μF[l]-ctr_bundles[n].soln.μF[l]) + ctr_bundles[n].eta_sg.μT[l]*(μT[l]-ctr_bundles[n].soln.μT[l]) for l in L)
        )
      end
    elseif CTR_PARAM == PROX 
	@variable(mMP, psi)
        @objective(mMP, Min, psi + 0.5*tVal*(
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
        @objective(mMP, Max, linobj - tVal*sum(alphaSlack[i] + betaSlack[i] + gammaSlack[i] + deltaSlack[i] for i in N))
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
        @objective(mMP, Max, linobj - tVal*Slack)
      else
        @objective(mMP, Max, linobj - 0.5*tVal*sum( ( ctr.soln.α[i] - α[i])^2 + (ctr.soln.β[i] - β[i])^2 + (ctr.soln.γ[i] - γ[i])^2 + (ctr.soln.δ[i] - δ[i])^2 for i in N))
      end
      if length(bundles) > 0
        @constraint(mMP, CutPlanes[n=1:length(bundles)], 0 <= -bundles[n].eta 
	  + sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i])
	    + bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
          + sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
	    + bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L)
        )
      end
    elseif CTR_PARAM == CP
      @objective(mMP, Max, linobj)
      if length(bundles) > 0
        @constraint(mMP, CutPlanes[n=1:length(bundles)], 0 <= -bundles[n].eta 
	  + sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i])
	    + bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
          + sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
	    + bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L)
        )
      end
    elseif CTR_PARAM == FEAS
      @variable(mMP, sLev >= 0)
      @constraint(mMP, LVLConstr, linobj + sLev >= nodeinfo.nodeBd)
      @objective(mMP, Min, sLev)
      if length(bundles) > 0
        @constraint(mMP, CutPlanes[n=1:length(bundles)], 0 <= -bundles[n].eta 
	  + sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i])
	    + bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
          + sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
	    + bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L)
        )
      end
    else
      @objective(mMP, Min, 0)
    end

  # Adding the extra cuts
      if CTR_PARAM == LVL2 || CTR_PARAM == LVL2 || CTR_PARAM == LVLINF
      elseif CTR_PARAM == PROX
      end
  ### END DEFINING THE LaGRANGIAN DUAL PROBLEM

  #status=solve(mMP)
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
      mpsoln.linerr = -ctr.eta + mpsoln.eta 
	  - dot( mpsoln.eta_sg.α[N], (ctr.soln.α[N]-mpsoln.soln.α[N]) )
	  - dot( mpsoln.eta_sg.β[N], (ctr.soln.β[N]-mpsoln.soln.β[N]) )
	  - dot( mpsoln.eta_sg.γ[N], (ctr.soln.γ[N]-mpsoln.soln.γ[N]) )
	  - dot( mpsoln.eta_sg.δ[N], (ctr.soln.δ[N]-mpsoln.soln.δ[N]) )
          - dot( mpsoln.eta_sg.λF[L],(ctr.soln.λF[L]-mpsoln.soln.λF[L]) )
          - dot( mpsoln.eta_sg.λT[L],(ctr.soln.λT[L]-mpsoln.soln.λT[L]) )
          - dot( mpsoln.eta_sg.μF[L],(ctr.soln.μF[L]-mpsoln.soln.μF[L]) )
          - dot( mpsoln.eta_sg.μT[L],(ctr.soln.μT[L]-mpsoln.soln.μT[L]) )

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

function testProxTraj(opfdata,K,HEUR)
    global tVal, tMin, tMax
    TOL = 1e-5
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    x_val=zeros(opfdata.nlines)
    x_val[41],x_val[80]=1,1
    #x_val[8],x_val[9],x_val[10],x_val[40]=1,1,1,1
    fixedNode=NodeInfo(x_val,x_val,1e20)
    bundles=Dict()
    ncuts=0
    rho = 0
    mpsoln=create_bundle(opfdata)
    ctr=create_bundle(opfdata)
    agg_bundles=Dict()
    agg_bundle=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)


    v_est,ssc_cntr,agg_norm = 1e20,0,0.0
  # MAIN LOOP
    for kk=1:maxNSG
     # PHASE I
      tMax = max(tMin,0.5*(comp_norm(opfdata,ctr.eta_sg)^2)/(1+abs(ctr.eta)))
      println("Penalty will be in [",tMin,",",tMax,"].")
      tVal = (tMax+tMin)/2
     # PHASE II
      while true
     # STEP 1
        mpsoln=create_bundle(opfdata)
        while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
	  tVal /= 2
          println("Resolving with reduced prox parameter value: ",tVal)
          status = solveNodeProxPt(opfdata,fixedNode,bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
        end	
        if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
	  ncuts=length(bundles)
	  naggcuts=length(bundles)
	  if ncuts > 0
	    rho = sum(bundles[n].cut_dual for n in 1:ncuts) + sum(agg_bundles[n].cut_dual for n in 1:naggcuts) 
           # COMPUTING LINEARIZATION ERRORS
	    agg_norm,epshat,rho=update_agg(opfdata,bundles,ctr,mpsoln,sg_agg,ctr_bundles,agg_bundles,agg_bundle)
            if rho*mpsoln.eta <= ssc*(1/tVal)*agg_norm^2 
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
      @show "Traj",kk,mpsoln.linobjval,mpsoln.eta,ncuts,agg_norm,tVal
      if agg_norm < 1e-2
        println("Convergence to within tolerance: obj, feas ",mpsoln.linobjval," ",mpsoln.eta)
	break
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

function testProxPt0(opfdata,K,HEUR)
    global tVal
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    x_val=zeros(opfdata.nlines)
    #x_val[41],x_val[80]=1,1
    x_val[8],x_val[9],x_val[10],x_val[40]=1,1,1,1
    x_lbs=zeros(opfdata.nlines)
    x_ubs=ones(opfdata.nlines)
    x_lbs[L]=x_val[L]
    x_ubs[L]=x_val[L]
    fixedNode=NodeInfo(x_lbs,x_ubs,1e20)
    trl_bundles=Dict()
    ctr_bundles=Dict()
    agg_bundles=Dict()
    ncuts=0
    rho,rhoUB = 0,0
    mpsoln=create_bundle(opfdata)
    ctr=create_bundle(opfdata)
    status = solveNodeProxPt(opfdata,fixedNode,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
    while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
      tVal /= 2
      println("Resolving with reduced prox parameter value: ",tVal)
      status = solveNodeProxPt(opfdata,fixedNode,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
    end	
    updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
    ctr_bundles[1]=mpsoln

    agg_bundle=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)

    v_est,ssc_cntr = 1e20,0
    TOL = 1e-5
  # MAIN LOOP
    for kk=1:maxNSG
     # STEP 1
      #tVal = 10/(10+length(trl_bundles))
      mpsoln=create_bundle(opfdata)
      status = solveNodeProxPt(opfdata,fixedNode,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
      while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
	tVal /= 2
        println("Resolving with reduced prox parameter value: ",tVal)
        status = solveNodeProxPt(opfdata,fixedNode,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
      end	
      if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
       # STEP 2
        ncuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
	 # UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
         # COMPUTING AGGREGATION INFORMATION
	  agg_norm,epshat,rho=update_agg(opfdata,trl_bundles,ctr,mpsoln,sg_agg,ctr_bundles,agg_bundles,agg_bundle)
	  if rhoUB < rho
	    rhoUB = rho + 1
	  end
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
          if ctr.eta < TOL && agg_norm < 1e-3 && epshat < 1e-3 
	    println("Convergence to within tolerance: obj, feas ",mpsoln.linobjval," ",mpsoln.eta)
	    @show "Basic",kk,mpsoln.linobjval,mpsoln.eta,ncuts,ctr.linobjval,ctr.eta,agg_norm,epshat,rho,rhoUB
	    break
          end
         # STEP 3
	  sscval = ((mpsoln.linobjval - rhoUB*mpsoln.eta)-(ctr.linobjval - rhoUB*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - rhoUB*ctr.eta)) 
	  vval = (mpsoln.linobjval-(ctr.linobjval - rhoUB*ctr.eta)) 
          if sscval >= ssc 
         # UPDATE CENTER VALUES
	    updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
	    nctrcuts=purgeSG(opfdata,ctr_bundles)
	    ctr_bundles[nctrcuts+1]=mpsoln
	    @show "Basic",kk,mpsoln.linobjval,mpsoln.eta,ncuts,agg_norm,epshat,ssc_cntr,tVal,rho,rhoUB
	  else
	    ncuts=purgeSG(opfdata,trl_bundles)
            trl_bundles[ncuts+1]=mpsoln
          end
	  #@show "Basic",kk,mpsoln.linobjval,mpsoln.eta,ncuts,agg_norm,epshat,ssc_cntr,tVal
          #tVal,v_est,ssc_cntr=KiwielRhoUpdate(opfdata,mpsoln,sscval,vval,agg_norm,epshat,rho,tVal,v_est,ssc_cntr)
	  #tVal = min(tVal,tMax)
      else
	println("Solver returned: $status")
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

function KiwielRhoUpdate(opfdata,mpsoln,sscval,vval,agg_norm,epshat,rho,tVal,v_est,ssc_cntr)
  global ssc
  nextTVal=tVal
  oldTVal=tVal
  if sscval >= ssc
    if sscval >= 0.5 && ssc_cntr > 0
      nextTVal = 2*tVal*(1-sscval)
    elseif ssc_cntr > 3
      nextTVal = tVal/2
    end
    tVal = max(nextTVal, tVal/10, tMin)
    v_est = max(v_est,2*vval)
    ssc_cntr = max(ssc_cntr+1,1)
    if oldTVal != tVal
      ssc_cntr = 1
    end
  else
    v_est=min(v_est,agg_norm+epshat)
#@show rho*mpsoln.linerr,max(v_est,10*vval),ssc_cntr
    if rho*mpsoln.linerr > max(v_est,10*vval) && ssc_cntr < -3
      nextTVal = 2*tVal*(1-sscval)
    end
    tVal = min(nextTVal,10*tVal)
    ssc_cntr = min(ssc_cntr-1,-1)
    if oldTVal != tVal
      ssc_cntr = -1	
    end
  end
  return tVal,v_est,ssc_cntr
end

function update_agg(opfdata,bundles,ctr,mpsoln,sg_agg,ctr_bundles,agg_bundles,agg_bundle)
  nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
  N, L, G = opfdata.N, opfdata.L, opfdata.G 
  comp_agg(opfdata,ctr.soln,mpsoln.soln,sg_agg)
  agg_norm=comp_norm(opfdata,sg_agg)
  rho = 0 #ctr.cut_dual
  ncuts,naggcuts=length(bundles),length(agg_bundles)
  if ncuts > 0
    #@show sum(bundles[n].cut_dual for n in 1:ncuts) 
    rho += sum(bundles[n].cut_dual for n in 1:ncuts) 
  end
  if naggcuts > 0
    #@show sum(agg_bundles[n].cut_dual for n in 1:naggcuts)
    rho += sum(agg_bundles[n].cut_dual for n in 1:naggcuts)
  end
  epshat = mpsoln.linobjval - (ctr.linobjval - rho*ctr.eta)
  - (dot(sg_agg.α[N],(mpsoln.soln.α[N]-ctr.soln.α[N])) - dot(sg_agg.β[N],(mpsoln.soln.β[N]-ctr.soln.β[N])) 
  + dot(sg_agg.γ[N],(mpsoln.soln.γ[N]-ctr.soln.γ[N])) - dot(sg_agg.δ[N],(mpsoln.soln.δ[N]-ctr.soln.δ[N])) 
  + dot(sg_agg.λF[L],(mpsoln.soln.λF[L] - ctr.soln.λF[L])) - dot(sg_agg.λT[L],(mpsoln.soln.λT[L] - ctr.soln.λT[L]))
  + dot(sg_agg.x[L], (mpsoln.soln.x[L] - ctr.soln.x[L]))
  + dot(sg_agg.μF[L],(mpsoln.soln.μF[L] - ctr.soln.μF[L])) - dot(sg_agg.μT[L],(mpsoln.soln.μT[L] - ctr.soln.μT[L])) )
  return agg_norm,epshat,rho
end

# Implements the approach of Sagastizabal and Solodov 2005
function testProxPt(opfdata,K,HEUR)
    global tVal
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    x_val=zeros(opfdata.nlines)
    #x_val[41],x_val[80]=1,1
    x_val[8],x_val[9],x_val[10],x_val[40]=1,1,1,1
    fixedNode=NodeInfo(x_val,x_val,1e20)
    bundles=Dict()
    mpsoln=create_bundle(opfdata)
    mpsoln.soln.x[L] = x_val[L]
    sg_agg=create_soln(opfdata)
    ctr=mpsoln

    v_est,ssc_cntr = 1e20,0
  # MAIN LOOP
    for kk=1:maxNSG
     # STEP 1
      mpsoln=create_bundle(opfdata)
      status = solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,PROX,mpsoln)
      while status != :Optimal
	tVal /= 2
        println("Resolving with reduced prox parameter value: ",tVal)
        status = solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,PROX,mpsoln)
      end	
      if status == :Optimal 
	comp_agg(opfdata,ctr.soln,mpsoln.soln,sg_agg)
	agg_norm=comp_norm(opfdata,sg_agg)
	epshat=max(0,ctr.eta)-mpsoln.psival-(1.0/tVal)*agg_norm^2
	#del=epshat+0.5*(1.0/tVal)*agg_norm^2
        del=ctr.eta-mpsoln.objval
       # STEP 2
        if del < 1e-6 
	  println("Convergence to within tolerance: del = ",del," val: ",mpsoln.linobjval," and feas: ",mpsoln.eta)
	  break
        end
       # STEP 3
	hk = max(ctr.linobjval-mpsoln.linobjval,mpsoln.eta)
	sscval = max(0,ctr.eta) - hk
        if hk <= max(0,ctr.eta) - ssc*del || kk==1
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
            @show "Sagadov",kk,ctr.linobjval,ctr.eta,epshat,del,hk,length(bundles),tVal
	else
	  ncuts = length(bundles)
          bundles[ncuts+1]=mpsoln
        end
       # STEP 4
        #tVal,v_est,ssc_cntr = KiwielRhoUpdate(opfdata,mpsoln,sscval,del,agg_norm,epshat,tVal,v_est,ssc_cntr)
      else
	#println("Solve status: $status")
	#break
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

#Implements the approach of De Oliveira 2016
function testLevelBM(opfdata,K,HEUR)
  global tVal
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    optUB=1e20
    x_val=zeros(opfdata.nlines)
    x_val[41],x_val[80]=1,1
    #x_val[8],x_val[9],x_val[10],x_val[40]=1,1,1,1
    lbs = zeros(opfdata.nlines)
    ubs = ones(opfdata.nlines)
    fixedNode=NodeInfo(x_val,x_val,optUB)
    #fixedNode=NodeInfo(lbs,ubs,optUB)
    # INITIAL SOLN (ALL ZEROS)
      mpsoln=create_bundle(opfdata)
    bestsoln = mpsoln
    ctr=mpsoln 
    bundles = Dict()
    prox_term=LVLINF

    cpsoln=create_bundle(opfdata)
    solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,CP,cpsoln)
    if cpsoln.status == :Optimal 
	optUB = cpsoln.objval
        hkctr = max(optUB,0)
    end

  # MAIN LOOP
    hkctr_updated = false
    optub_updated = false
    for kk=1:maxNSG
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
      fixedNode.nodeBd = optUB - ssc*hkval
      if hkval < 1e-6
	println("Convergence to within tolerance.")
	break
      end

     # STEP 2
        if hkval <= (1-ssc)*hkctr 
          hkctr = hkval
          hkctr_updated = true
          # UPDATE CENTER VALUES
	    if bestIdx > 0
	      updateCenter(opfdata,bestsoln,ctr,bundles)
	    end
	  @show "hk  update",optUB,ctr.linobjval,ctr.eta,hkval,ncuts
        end

     # STEP 3
      cpstatus = solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,CP,cpsoln)
      if cpsoln.status != :Optimal
	println("FLAG! cpsoln.status = ",cpsoln.status)
      end
      if fixedNode.nodeBd-cpsoln.linobjval <= 0.0
	tVal = 0.1
        mpsoln=create_bundle(opfdata)
        status=solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,prox_term,mpsoln)
        while status != :Optimal
	  tVal /= 2
          println("Resolving with reduced prox parameter value: ",tVal)
          status=solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,prox_term,mpsoln)
        end	
	if mpsoln.status != :Optimal
	  println("Taking recourse since mpsoln does not have an optimal solution for MP. Feas: ",fixedNode.nodeBd - mpsoln.linobjval," val: ",fixedNode.nodeBd - cpsoln.linobjval)
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
	  optUB = fixedNode.nodeBd 
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


#USEFUL SUBROUTINES
function updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
  cpy_bundle(opfdata,mpsoln,ctr)
  for n=1:length(trl_bundles)
    updateLinErr(opfdata,ctr,trl_bundles[n])
  end
  for n=1:length(ctr_bundles)
    updateLinErr(opfdata,ctr,ctr_bundles[n])
  end
  for n=1:length(agg_bundles)
    updateLinErr(opfdata,ctr,agg_bundles[n])
  end
end

function updateLinErr(opfdata,ctr,bundle)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  bundle.linerr = ctr.eta - bundle.eta 
    - dot(bundle.eta_sg.α[N],(bundle.soln.α[N] - ctr.soln.α[N]))
    - dot(bundle.eta_sg.β[N],(bundle.soln.β[N] - ctr.soln.β[N]))
    - dot(bundle.eta_sg.γ[N],(bundle.soln.γ[N] - ctr.soln.γ[N]))
    - dot(bundle.eta_sg.δ[N],(bundle.soln.δ[N] - ctr.soln.δ[N]))
    - dot(bundle.eta_sg.λF[L],(bundle.soln.λF[L] - ctr.soln.λF[L]))
    - dot(bundle.eta_sg.λT[L],(bundle.soln.λT[L] - ctr.soln.λT[L]))
    - dot(bundle.eta_sg.μF[L],(bundle.soln.μF[L] - ctr.soln.μF[L]))
    - dot(bundle.eta_sg.μT[L],(bundle.soln.μT[L] - ctr.soln.μT[L]))
end

function computeSG(opfdata,mpsoln)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      vR = zeros(nbuses)
      vI = zeros(nbuses)
      try
        mpsoln.eta = -solveEta0Eigs(opfdata,mpsoln.soln,vR,vI)
      catch exc
        println("Exception caught with eigs(), solving η0Val subproblem with Ipopt as recourse.")
        println(exc)
        mpsoln.eta = -solveEta0SDP(opfdata,mpsoln.soln,vR,vI)
      end
      for i in N
        W_val = vR[i]^2 + vI[i]^2
        mpsoln.eta_sg.α[i] = Y["shR"][i] * W_val
	mpsoln.eta_sg.β[i] = -Y["shI"][i] * W_val
        mpsoln.eta_sg.δ[i] = W_val
	mpsoln.eta_sg.γ[i] = -W_val
      end
      for l in L
        from = fromBus[l]; to = toBus[l]
        e_valF = vR[from]; f_valF = vI[from]; W_valF = e_valF^2 + f_valF^2
        e_valT = vR[to]; f_valT = vI[to]; W_valT = e_valT^2 + f_valT^2
        Wr_val = e_valF*e_valT + f_valF*f_valT; Wi_val = e_valT*f_valF - e_valF*f_valT
        mpsoln.eta_sg.λF[l] = (Y["ffR"][l] * W_valF + Y["ftR"][l] * Wr_val + Y["ftI"][l] * Wi_val)
        mpsoln.eta_sg.λT[l] = (Y["ttR"][l] * W_valT + Y["tfR"][l] * Wr_val - Y["tfI"][l] * Wi_val)
        mpsoln.eta_sg.μF[l] = (-Y["ffI"][l] * W_valF - Y["ftI"][l] * Wr_val + Y["ftR"][l] * Wi_val)
        mpsoln.eta_sg.μT[l] = (-Y["ttI"][l] * W_valT - Y["tfI"][l] * Wr_val - Y["tfR"][l] * Wi_val)
      end
      return mpsoln.eta
end
# SUBROUTINE FOR COMPUTING THE MINIMUM EIGENVALUE OF H WITH A CORRESPONDING EIGENVECTOR
function solveEta0Eigs(opfdata,soln,vR,vI)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus
      H=spzeros(2*nbuses,2*nbuses)
      updateHess(opfdata,soln,H)
      E=eigs(H,nev=1,which=:SR, maxiter=100000, tol=1e-6)
      η0Val = E[1][1]
      for i in N
        vR[i] = E[2][i,1]; vI[i] = E[2][nbuses+i,1]
      end
      return η0Val
end
# Update Hessian
function updateHess(opfdata,pi_val,H)
      #lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus  
      Y = opfdata.Y_AC
      for i in N
        H[i,i] +=  pi_val.α[i] * Y["shR"][i] - pi_val.β[i] * Y["shI"][i]  + pi_val.δ[i] - pi_val.γ[i]
        H[nbuses+i,nbuses+i] += pi_val.α[i] * Y["shR"][i] - pi_val.β[i] * Y["shI"][i] + pi_val.δ[i] - pi_val.γ[i]
      end
      for l in L
        from = fromBus[l]; to = toBus[l]
        H[from,from] += pi_val.λF[l] * Y["ffR"][l] - pi_val.μF[l] * Y["ffI"][l]
        H[nbuses+from,nbuses+from] += pi_val.λF[l] * Y["ffR"][l] - pi_val.μF[l] * Y["ffI"][l]
        H[to,to] += pi_val.λT[l] * Y["ttR"][l] - pi_val.μT[l] * Y["ttI"][l]
        H[nbuses+to,nbuses+to] += pi_val.λT[l] * Y["ttR"][l] - pi_val.μT[l] * Y["ttI"][l]
        H[from,to] += 0.5*( pi_val.λF[l] * Y["ftR"][l] - pi_val.μF[l] * Y["ftI"][l] + pi_val.λT[l] * Y["tfR"][l] - pi_val.μT[l] * Y["tfI"][l] )
        H[to,from] += 0.5*( pi_val.λF[l] * Y["ftR"][l] - pi_val.μF[l] * Y["ftI"][l] + pi_val.λT[l] * Y["tfR"][l] - pi_val.μT[l] * Y["tfI"][l] )
        H[nbuses+from, nbuses+to] += 0.5*( pi_val.λF[l] * Y["ftR"][l] - pi_val.μF[l] * Y["ftI"][l] + pi_val.λT[l] * Y["tfR"][l] - pi_val.μT[l] * Y["tfI"][l] )
        H[nbuses+to, nbuses+from] += 0.5*( pi_val.λF[l] * Y["ftR"][l] - pi_val.μF[l] * Y["ftI"][l] + pi_val.λT[l] * Y["tfR"][l] - pi_val.μT[l] * Y["tfI"][l] )
        H[to, nbuses+from] += 0.5*( pi_val.λF[l] * Y["ftI"][l] - pi_val.λT[l] * Y["tfI"][l] + pi_val.μF[l] * Y["ftR"][l] - pi_val.μT[l] * Y["tfR"][l] )
        H[nbuses+from, to] += 0.5*( pi_val.λF[l] * Y["ftI"][l] - pi_val.λT[l] * Y["tfI"][l] + pi_val.μF[l] * Y["ftR"][l] - pi_val.μT[l] * Y["tfR"][l] )
        H[from,nbuses+to] -= 0.5*( pi_val.λF[l] * Y["ftI"][l] - pi_val.λT[l] * Y["tfI"][l] + pi_val.μF[l] * Y["ftR"][l] - pi_val.μT[l] * Y["tfR"][l] )
        H[nbuses+to,from] -= 0.5*( pi_val.λF[l] * Y["ftI"][l] - pi_val.λT[l] * Y["tfI"][l] + pi_val.μF[l] * Y["ftR"][l] - pi_val.μT[l] * Y["tfR"][l] )
      end
end


# SUBROUTINE FOR COMPUTING THE MINIMUM EIGENVALUE OF H WITH A CORRESPONDING EIGENVECTOR
  # VIA AN OPTIMIZATION PROBLEM
function solveEta0SDP(opfdata,soln,vR,vI)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus

      #The QP subproblem
      mSDP = Model(solver=IpoptSolver())
      @variable(mSDP, e[i=N], start=0); @variable(mSDP, f[i=N], start=0)
      η0Val = 0

      for i in N
        setvalue(e[i], 1); setvalue(f[i], 0)
      end

      H=spzeros(2*nbuses,2*nbuses)
      updateHess(opfdata,soln,H)

      # Adjust QP subproblem
      @NLobjective(mSDP, Min, sum( H[i,i]*(e[i]^2+f[i]^2) for i in N)
        + 2*sum( H[fromBus[l],toBus[l]]*(e[fromBus[l]]*e[toBus[l]]+f[fromBus[l]]*f[toBus[l]])   for l in L)
        - 2*sum( H[fromBus[l],nbuses+toBus[l]]*(f[fromBus[l]]*e[toBus[l]]-e[fromBus[l]]*f[toBus[l]])   for l in L)
      )
      status = solve(mSDP)
      if status == :Optimal || status == :UserLimit
        η0Val = getobjectivevalue(mSDP)
        for i in N
          vR[i],vI[i]=getvalue(e[i]),getvalue(f[i])
        end
        if(status == :UserLimit)
          println("solveEta0SDP solve status $status")
        end
      else
        println("solveEta0SDP: solve status $status")
        η0Val = 0
      end
      return η0Val
end

# SUBROUTINE FOR COMPUTING A SUBGRADIENT OF ETA(PI), WHICH IS THE FUNCTION TAKING THE VALUE OF THE MINIMUM EIGENVALUE OF H(PI)
function purgeSG(opfdata,bundle,minAge=20,maxAge=80)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC

      orig_ncuts = length(bundle)
      ncuts = orig_ncuts
      for n=orig_ncuts:-1:1
        if (abs(bundle[n].cut_dual) < 1e-8 && bundle[n].age > minAge) || bundle[n].age > maxAge
	  bundle[n]=bundle[ncuts]
	  delete!(bundle,ncuts)
	  ncuts -= 1
	end
      end
      for n=1:length(bundle)
	bundle[n].age += 1
      end
      return length(bundle)
end

function aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      agg_bundle=create_bundle(opfdata)

      cpy_soln(opfdata,mpsoln.soln,agg_bundle.soln)
      sumDuals,ncuts,nctrcuts,naggcuts = 0.0,length(trl_bundles),length(ctr_bundles),length(agg_bundles)
      if ncuts > 0
        sumDuals += sum(trl_bundles[n].cut_dual for n in 1:ncuts) 
      end
      if nctrcuts > 0
        sumDuals += sum(ctr_bundles[n].cut_dual for n in 1:nctrcuts) 
      end
      if naggcuts > 0
        sumDuals += sum(agg_bundles[n].cut_dual for n in 1:naggcuts) 
      end

      agg_bundle.eta_sg.α[N] .= 0
      agg_bundle.eta_sg.β[N] .= 0
      agg_bundle.eta_sg.γ[N] .= 0
      agg_bundle.eta_sg.δ[N] .= 0
      agg_bundle.eta_sg.ζpLB[G] .= 0
      agg_bundle.eta_sg.ζpUB[G] .= 0
      agg_bundle.eta_sg.ζqLB[G] .= 0
      agg_bundle.eta_sg.ζqUB[G] .= 0
      agg_bundle.eta_sg.x[L] .= 0
      agg_bundle.eta_sg.λF[L] .= 0
      agg_bundle.eta_sg.λT[L] .= 0
      agg_bundle.eta_sg.μF[L] .= 0
      agg_bundle.eta_sg.μT[L] .= 0

      agg_bundle.objval = 0
      agg_bundle.linobjval = 0
      agg_bundle.penval = 0
      agg_bundle.psival = 0

      if ncuts > 0
        agg_bundle.eta_sg.α[N] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.α[N] for n in 1:ncuts) )
        agg_bundle.eta_sg.β[N] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.β[N] for n in 1:ncuts) )
        agg_bundle.eta_sg.γ[N] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.γ[N] for n in 1:ncuts) )
        agg_bundle.eta_sg.δ[N] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.δ[N] for n in 1:ncuts) )
        agg_bundle.eta_sg.ζpLB[G] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.ζpLB[G] for n in 1:ncuts) )
        agg_bundle.eta_sg.ζpUB[G] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.ζpUB[G] for n in 1:ncuts) )
        agg_bundle.eta_sg.ζqLB[G] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.ζqLB[G] for n in 1:ncuts) )
        agg_bundle.eta_sg.ζqUB[G] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.ζqUB[G] for n in 1:ncuts) )
        agg_bundle.eta_sg.x[L]  += (1/sumDuals)* ( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.x[L] for n in 1:ncuts) )
        agg_bundle.eta_sg.λF[L] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.λF[L] for n in 1:ncuts) )
        agg_bundle.eta_sg.λT[L] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.λT[L] for n in 1:ncuts) )
        agg_bundle.eta_sg.μF[L] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.μF[L] for n in 1:ncuts) )
        agg_bundle.eta_sg.μT[L] += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].eta_sg.μT[L] for n in 1:ncuts) )

        agg_bundle.objval += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].objval for n in 1:ncuts) )
        agg_bundle.linobjval += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].linobjval for n in 1:ncuts) )
        agg_bundle.penval += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].penval for n in 1:ncuts) )
        agg_bundle.psival += (1/sumDuals)*( sum( trl_bundles[n].cut_dual*trl_bundles[n].psival for n in 1:ncuts) )

        agg_bundle.etahat -= (1/sumDuals)*sum( trl_bundles[n].cut_dual*(trl_bundles[n].eta 
	+ dot(trl_bundles[n].eta_sg.α[N],trl_bundles[n].soln.α[N]) + dot(trl_bundles[n].eta_sg.β[N],trl_bundles[n].soln.β[N])
	  + dot(trl_bundles[n].eta_sg.γ[N],trl_bundles[n].soln.γ[N]) + dot(trl_bundles[n].eta_sg.δ[N],trl_bundles[n].soln.δ[N]) 
        + dot(trl_bundles[n].eta_sg.λF[L],trl_bundles[n].soln.λF[L]) + dot(trl_bundles[n].eta_sg.λT[L],trl_bundles[n].soln.λT[L]) 
	  + dot(trl_bundles[n].eta_sg.μF[L],trl_bundles[n].soln.μF[L]) + dot(trl_bundles[n].eta_sg.μT[L],trl_bundles[n].soln.μT[L]) ) for n in 1:ncuts)
      end
      if nctrcuts > 0
        agg_bundle.eta_sg.α[N] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.α[N] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.β[N] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.β[N] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.γ[N] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.γ[N] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.δ[N] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.δ[N] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.ζpLB[G] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.ζpLB[G] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.ζpUB[G] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.ζpUB[G] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.ζqLB[G] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.ζqLB[G] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.ζqUB[G] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.ζqUB[G] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.x[L]  += (1/sumDuals)* ( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.x[L] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.λF[L] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.λF[L] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.λT[L] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.λT[L] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.μF[L] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.μF[L] for n in 1:nctrcuts) )
        agg_bundle.eta_sg.μT[L] += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].eta_sg.μT[L] for n in 1:nctrcuts) )

        agg_bundle.objval += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].objval for n in 1:nctrcuts) )
        agg_bundle.linobjval += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].linobjval for n in 1:nctrcuts) )
        agg_bundle.penval += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].penval for n in 1:nctrcuts) )
        agg_bundle.psival += (1/sumDuals)*( sum( ctr_bundles[n].cut_dual*ctr_bundles[n].psival for n in 1:nctrcuts) )

        agg_bundle.etahat -= (1/sumDuals)*sum( ctr_bundles[n].cut_dual*(ctr_bundles[n].eta 
	+ dot(ctr_bundles[n].eta_sg.α[N],ctr_bundles[n].soln.α[N]) + dot(ctr_bundles[n].eta_sg.β[N],ctr_bundles[n].soln.β[N])
	  + dot(ctr_bundles[n].eta_sg.γ[N],ctr_bundles[n].soln.γ[N]) + dot(ctr_bundles[n].eta_sg.δ[N],ctr_bundles[n].soln.δ[N]) 
        + dot(ctr_bundles[n].eta_sg.λF[L],ctr_bundles[n].soln.λF[L]) + dot(ctr_bundles[n].eta_sg.λT[L],ctr_bundles[n].soln.λT[L]) 
	  + dot(ctr_bundles[n].eta_sg.μF[L],ctr_bundles[n].soln.μF[L]) + dot(ctr_bundles[n].eta_sg.μT[L],ctr_bundles[n].soln.μT[L]) ) for n in 1:nctrcuts)
      end
      if naggcuts > 0
        agg_bundle.eta_sg.α[N] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.α[N] for n in 1:naggcuts) )
        agg_bundle.eta_sg.β[N] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.β[N] for n in 1:naggcuts) )
        agg_bundle.eta_sg.γ[N] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.γ[N] for n in 1:naggcuts) )
        agg_bundle.eta_sg.δ[N] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.δ[N] for n in 1:naggcuts) )
        agg_bundle.eta_sg.ζpLB[G] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.ζpLB[G] for n in 1:naggcuts) )
        agg_bundle.eta_sg.ζpUB[G] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.ζpUB[G] for n in 1:naggcuts) )
        agg_bundle.eta_sg.ζqLB[G] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.ζqLB[G] for n in 1:naggcuts) )
        agg_bundle.eta_sg.ζqUB[G] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.ζqUB[G] for n in 1:naggcuts) )
        agg_bundle.eta_sg.x[L]  += (1/sumDuals)* ( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.x[L] for n in 1:naggcuts) )
        agg_bundle.eta_sg.λF[L] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.λF[L] for n in 1:naggcuts) )
        agg_bundle.eta_sg.λT[L] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.λT[L] for n in 1:naggcuts) )
        agg_bundle.eta_sg.μF[L] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.μF[L] for n in 1:naggcuts) )
        agg_bundle.eta_sg.μT[L] += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].eta_sg.μT[L] for n in 1:naggcuts) )

        agg_bundle.objval += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].objval for n in 1:naggcuts) )
        agg_bundle.linobjval += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].linobjval for n in 1:naggcuts) )
        agg_bundle.penval += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].penval for n in 1:naggcuts) )
        agg_bundle.psival += (1/sumDuals)*( sum( agg_bundles[n].cut_dual*agg_bundles[n].psival for n in 1:naggcuts) )

        agg_bundle.etahat -= (1/sumDuals)*sum( agg_bundles[n].cut_dual*(agg_bundles[n].eta 
	+ dot(agg_bundles[n].eta_sg.α[N],agg_bundles[n].soln.α[N]) + dot(agg_bundles[n].eta_sg.β[N],agg_bundles[n].soln.β[N])
	  + dot(agg_bundles[n].eta_sg.γ[N],agg_bundles[n].soln.γ[N]) + dot(agg_bundles[n].eta_sg.δ[N],agg_bundles[n].soln.δ[N]) 
        + dot(agg_bundles[n].eta_sg.λF[L],agg_bundles[n].soln.λF[L]) + dot(agg_bundles[n].eta_sg.λT[L],agg_bundles[n].soln.λT[L]) 
	  + dot(agg_bundles[n].eta_sg.μF[L],agg_bundles[n].soln.μF[L]) + dot(agg_bundles[n].eta_sg.μT[L],agg_bundles[n].soln.μT[L]) ) for n in 1:naggcuts)
      end
      agg_bundle.age = 1
      agg_bundle.eta = mpsoln.eta
      return agg_bundle
end

