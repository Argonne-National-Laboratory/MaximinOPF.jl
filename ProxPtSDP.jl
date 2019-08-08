#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")

maxNSG = 20000
CP,PROX0,PROX,LVL1,LVL2,LVLINF,FEAS=0,1,2,3,4,5,6
tVal = 1
ssc=0.5

type NodeInfo
  x_lbs::Array{Float64}
  x_ubs::Array{Float64}
  nodeBd::Float64
end
type Soln
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

type Bundle
  soln::Soln
  eta_sg::Soln
  objval::Float64
  linobjval::Float64
  penval::Float64
  psival::Float64
  eta::Float64
  cut_dual::Float64
  lvl_dual::Float64
  solvetime::Float64
  status
end
function create_bundle(opfdata)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  return Bundle(create_soln(opfdata),create_soln(opfdata),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0)
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
  toBundle.cut_dual = fromBundle.cut_dual
  toBundle.lvl_dual = fromBundle.lvl_dual
  toBundle.solvetime = fromBundle.solvetime
  toBundle.status = fromBundle.status
end

#=
function evalConstraints(opfdata,soln,K)
  nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
  fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
  BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC

  maxViolation = soln.eta
  xUBV = soln.x - 1
  xLBV = 0 - soln.x
  αUBV = soln.α - 1 
  αLBV = -1 - soln.α 
  βUBV = soln.β - 1 
  βLBV = -1 - soln.β 
  γLBV = 0 - soln.γ 
  δLBV = 0 - soln.δ 
  γδUBV = soln.γ+soln.δ - 1 
  ζpUBV = soln.ζpUB - 1 
  ζpLBV = 1 - soln.ζpLB  
  ζqUBV = soln.ζqUB - 1 
  ζqLBV = 1 - soln.ζqLB
  maxViolation = max(maximum(xUBV),maximum(xLBV),maximum(αUBV),maximum(αLBV),maximum(βUBV),maximum(βLBV),
	maximum(γLBV),maximum(δLBV),maximum(γδUBV),maximum(ζpUBV),maximum(ζpLBV),maximum(ζqUBV),maximum(ζqLBV),maxViolation)

  for i in N
    for g in BusGeners[i]
      maxViolation = max(abs( -soln.α[i] + soln.ζpUB[g] - soln.ζpLB[g] ),maxViolation)
      maxViolation = max(abs( -soln.β[i] + soln.ζqUB[g] - soln.ζqLB[g] ),maxViolation)
    end
  end
  maxViolation = max(sum(soln.x) - K, maxViolation)
  
end

function solveNodeProxPtSagadov(opfdata,nodeinfo,bundles,K,HEUR,ctr,CTR_PARAM,
			mpsoln)
  global tVal
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata


  # The master problem MP
    #if CTR_PARAM == CP || CTR_PARAM == FEAS
      mMP = Model(solver=CplexSolver(
        CPX_PARAM_SCRIND=0,
        CPX_PARAM_TILIM=MAX_TIME,
        CPX_PARAM_THREADS=1
	)
      )
    @variable(mMP, psi)
    # each x[l] is either 0 (line l active) or 1 (line l is cut)
      @variable(mMP, 0 <= x[l=L] <= 1, start=0)
      @constraint(mMP, [l=L], x[l] - nodeinfo.x_ubs[l] <= psi)
      @constraint(mMP, [l=L], nodeinfo.x_lbs[l] -x[l] <= psi)
    # dual multipliers associated with active power flow balance constraints
      @variable(mMP, -1 <= α[i=N] <= 1, start=0)
      @constraint(mMP, [i=N], α[i] - 1 <= psi)
      @constraint(mMP, [i=N], -1 - α[i] <= psi)
    # dual multipliers associated with active power flow balance constraints
      @variable(mMP, -1 <= β[i=N] <= 1, start=0)
      @constraint(mMP, [i=N], β[i] - 1 <= psi)
      @constraint(mMP, [i=N], -1 - β[i] <= psi)
    # dual multipliers associated with the voltage magnitude bounds
      @variable(mMP, δ[i=N] >= 0) 
      @constraint(mMP, [i=N], -δ[i] <= psi)
      @variable(mMP, γ[i=N] >= 0) 
      @constraint(mMP, [i=N], -γ[i] <= psi)
      @constraint(mMP, [i=N], δ[i]+γ[i] -1 <= psi)
    # dual multipliers associated with the power generation bounds
      @variable(mMP, 0 <= ζpUB[g=G] <= 1); @variable(mMP, 0 <= ζpLB[g=G] <= 1); @variable(mMP, 0 <= ζqUB[g=G] <= 1); @variable(mMP, 0 <= ζqLB[g=G] <= 1)
      @constraint(mMP, [g=G], ζpUB[g] - 1 <= psi)
      @constraint(mMP, [g=G], -ζpUB[g] <= psi)
      @constraint(mMP, [g=G], ζpLB[g] - 1 <= psi)
      @constraint(mMP, [g=G], -ζpLB[g] <= psi)
      @constraint(mMP, [g=G], ζqUB[g] - 1 <= psi)
      @constraint(mMP, [g=G], -ζqUB[g] <= psi)
      @constraint(mMP, [g=G], ζqLB[g] - 1 <= psi)
      @constraint(mMP, [g=G], -ζqLB[g] <= psi)
    # constraints needed for terms in pG[g] and qG[g] in the Lagrangian to vanish (needed for dual feasibility)
      for i in N
        for g in BusGeners[i]
          @constraint(mMP, -α[i] + ζpUB[g] - ζpLB[g] == 0 )
          @constraint(mMP, -β[i] + ζqUB[g] - ζqLB[g] == 0 )
        end
      end
  @constraint(mMP, sum(x[l] for l in L) - K <= psi)
 # McCormick inequalities enforcing bilinear equalities
    # auxiliary dual variables due to McCormick reformulation of cross terms appearing in the Lagrangian
      @variable(mMP, λF[l=L], start=0); @variable(mMP, λT[l=L], start=0); @variable(mMP, μF[l=L], start=0); @variable(mMP, μT[l=L], start=0)
    @constraint(mMP, AMcf1[l in L], α[fromBus[l]] - x[l] - λF[l] <= psi); @constraint(mMP, AMcf2[l in L], α[fromBus[l]] + x[l] - λF[l] >= -psi)
    @constraint(mMP, AMcf3[l in L], -(1 - x[l]) - λF[l] <= psi); @constraint(mMP, AMcf4[l in L], (1 - x[l]) - λF[l] >= -psi)
    @constraint(mMP, AMct1[l in L], α[toBus[l]] - x[l] - λT[l] <= psi); @constraint(mMP, AMct2[l in L], α[toBus[l]] + x[l] - λT[l] >= -psi)
    @constraint(mMP, AMct3[l in L], -(1 - x[l]) - λT[l] <= psi); @constraint(mMP, AMct4[l in L], (1 - x[l]) - λT[l] >= -psi)

    @constraint(mMP, BMcf1[l in L], β[fromBus[l]] - x[l] - μF[l] <= psi); @constraint(mMP, BMcf2[l in L], β[fromBus[l]] + x[l] - μF[l] >= -psi)
    @constraint(mMP, BMcf3[l in L], -(1 - x[l]) - μF[l] <= psi); @constraint(mMP, BMcf4[l in L], (1 - x[l]) - μF[l] >= -psi)
    @constraint(mMP, BMct1[l in L], β[toBus[l]] - x[l] - μT[l] <= psi); @constraint(mMP, BMct2[l in L], β[toBus[l]] + x[l] - μT[l] >= -psi)
    @constraint(mMP, BMct3[l in L], -(1 - x[l]) - μT[l] <= psi); @constraint(mMP, BMct4[l in L], (1 - x[l]) - μT[l] >= -psi)

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
    #@constraint(mMP, ctr.linobjval - linobj <= psi)
    @constraint(mMP, 0.937 - linobj <= psi)
   
    @objective(mMP, Min, psi + 0.5*tVal*(
      sum( (ctr.soln.α[i] - α[i])^2 + (ctr.soln.β[i] - β[i])^2 + (ctr.soln.γ[i] - γ[i])^2 + (ctr.soln.δ[i] - δ[i])^2 for i in N)
	+sum( (ctr.soln.ζpLB[g] - ζpLB[g])^2 + (ctr.soln.ζpUB[g] - ζpUB[g])^2 + (ctr.soln.ζqLB[g] - ζqLB[g])^2 + (ctr.soln.ζqUB[g] - ζqUB[g])^2 for g in G)
	+sum( (ctr.soln.x[l]-x[l])^2 + (ctr.soln.λF[l] - λF[l])^2 + (ctr.soln.λT[l] - λT[l])^2 + (ctr.soln.μF[l] - μF[l])^2 + (ctr.soln.μT[l] - μT[l])^2 for l in L)
      )
    )
    for n=1:length(bundles)
      @constraint(mMP, bundles[n].eta - sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i]) 
	+ bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
        - sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
	+ bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L) - psi <= 0.0
      )
    end
    status=solve(mMP)
    mpsoln.status=status
    if status == :Optimal 
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
      mpsoln.psival = getvalue(psi)

      computeSG(opfdata,mpsoln) #This computes mpsoln.eta 

    else
      #println("solveNodeAC: Return status $status")
    end
    return status
end
=#

function solveNodeProxPt(opfdata,nodeinfo,bundles,K,HEUR,ctr,CTR_PARAM,
			mpsoln)
  global tVal
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata


  # The master problem MP
    #if CTR_PARAM == CP || CTR_PARAM == FEAS
      mMP = Model(solver=CplexSolver(
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
    #else
    #mMP = Model(solver=IpoptSolver())
    #end
    #mMP = Model(solver=CplexSolver(CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_LPMETHOD=4,CPX_PARAM_SOLUTIONTYPE=2,CPX_PARAM_STARTALG=4))
    # each x[l] is either 0 (line l active) or 1 (line l is cut)
      @variable(mMP, 0 <= x[l=L] <= 1, start=0)
      for l in L
        setlowerbound(x[l],nodeinfo.x_lbs[l])
        setupperbound(x[l],nodeinfo.x_ubs[l])
      end
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
    @constraint(mMP, AMcf1[l in L], α[fromBus[l]] - x[l] <= λF[l]); @constraint(mMP, AMcf2[l in L], α[fromBus[l]] + x[l] >= λF[l])
    @constraint(mMP, AMcf3[l in L], -(1 - x[l]) <= λF[l]); @constraint(mMP, AMcf4[l in L], (1 - x[l]) >= λF[l])
    @constraint(mMP, AMct1[l in L], α[toBus[l]] - x[l] <= λT[l]); @constraint(mMP, AMct2[l in L], α[toBus[l]] + x[l] >= λT[l])
    @constraint(mMP, AMct3[l in L], -(1 - x[l]) <= λT[l]); @constraint(mMP, AMct4[l in L], (1 - x[l]) >= λT[l])

    @constraint(mMP, BMcf1[l in L], β[fromBus[l]] - x[l] <= μF[l]); @constraint(mMP, BMcf2[l in L], β[fromBus[l]] + x[l] >= μF[l])
    @constraint(mMP, BMcf3[l in L], -(1 - x[l]) <= μF[l]); @constraint(mMP, BMcf4[l in L], (1 - x[l]) >= μF[l])
    @constraint(mMP, BMct1[l in L], β[toBus[l]] - x[l] <= μT[l]); @constraint(mMP, BMct2[l in L], β[toBus[l]] + x[l] >= μT[l])
    @constraint(mMP, BMct3[l in L], -(1 - x[l]) <= μT[l]); @constraint(mMP, BMct4[l in L], (1 - x[l]) >= μT[l])

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
	+sum((ctr.soln.λF[l] - λF[l])^2 + (ctr.soln.λT[l] - λT[l])^2 + (ctr.soln.μF[l] - μF[l])^2 + (ctr.soln.μT[l] - μT[l])^2 for l in L)
	)
      )
      if length(bundles) > 0
        @constraint(mMP, CutPlanes[n=1:length(bundles)], 0 <= -bundles[n].eta 
	  + sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i])
	    + bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
          + sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
	    + bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L)
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
	addnncon=true
	addlincon=true
	#if ctr.linobjval >= 0
	   @constraint(mMP, ctr.linobjval - linobj - psi <= 0)
	   #addlincon = false
	#end
	for n=1:length(bundles)
	  #if ctr.linobjval-bundles[n].linobjval < bundles[n].eta
	  #@show ctr.linobjval-bundles[n].linobjval,bundles[n].eta
	    #if bundles[n].eta == 0.0
	      #@constraint(mMP, psi >= 0)
	    #  addnncon=false
	    #if bundles[n].eta > 0.0
              @constraint(mMP, bundles[n].eta - sum( bundles[n].eta_sg.α[i]*(α[i]-bundles[n].soln.α[i]) + bundles[n].eta_sg.β[i]*(β[i]-bundles[n].soln.β[i]) 
	        + bundles[n].eta_sg.γ[i]*(γ[i]-bundles[n].soln.γ[i]) + bundles[n].eta_sg.δ[i]*(δ[i]-bundles[n].soln.δ[i]) for i in N)
                - sum( bundles[n].eta_sg.λF[l]*(λF[l]-bundles[n].soln.λF[l]) + bundles[n].eta_sg.λT[l]*(λT[l]-bundles[n].soln.λT[l]) 
		+ bundles[n].eta_sg.μF[l]*(μF[l]-bundles[n].soln.μF[l]) + bundles[n].eta_sg.μT[l]*(μT[l]-bundles[n].soln.μT[l]) for l in L) - psi <= 0.0
              )
	    #end
	  #else
#=
	    if addlincon
	      @constraint(mMP, ctr.linobjval - linobj - psi <= 0)
	      addlincon = false
	    end
	    #println("Not adding cut")
=#
	  #end
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

  status=solve(mMP)
  mpsoln.status=status
  if status == :Optimal || status == :CPX_STAT_NUM_BEST
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

      if status == :Optimal && CTR_PARAM!=PROX
	for n=1:length(bundles)
	  bundles[n].cut_dual = getdual(CutPlanes[n])
	end
      end
      if status == :Optimal && (CTR_PARAM==LVL1 || CTR_PARAM == LVL2 || CTR_PARAM == LVLINF )
	mpsoln.lvl_dual = -getdual(LVLConstr)
      end
  else
    #println("solveNodeAC: Return status $status")
  end
  return status
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
    fixedNode=NodeInfo(x_val,x_val,1e20)
    bundles=Dict()
    ncuts=0
    mpsoln=create_bundle(opfdata)
    ctr=create_bundle(opfdata)
    status = solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,PROX0,mpsoln)
    if status == :Optimal
      if mpsoln.eta < 1e-6 
	  println("Convergence to within tolerance: obj, feas ",mpsoln.linobjval,mpsoln.eta)
          time_End = (time_ns()-time_Start)/1e9
          println("Done after ",time_End," seconds.")
	  return
      end
      updateCenter(opfdata,mpsoln,ctr)
      bundles[1]=mpsoln
      ncuts=1
    else
	println("Solver returned: $status")
    end
    @show 0,mpsoln.linobjval,mpsoln.eta

  # MAIN LOOP
    for kk=1:maxNSG
     # STEP 1
      mpsoln=create_bundle(opfdata)
      status = solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,PROX0,mpsoln)
      if status == :Optimal
       # STEP 2
        if mpsoln.eta < 1e-6 
	  println("Convergence to within tolerance: obj, feas ",mpsoln.linobjval,mpsoln.eta)
          time_End = (time_ns()-time_Start)/1e9
          println("Done after ",time_End," seconds.")
	  return
        end
       # STEP 3
        ncuts = length(bundles)
        if -(mpsoln.eta-ctr.eta)/ctr.eta >= 0.05 
          # UPDATE CENTER VALUES
	    updateCenter(opfdata,mpsoln,ctr)
	    ncuts=purgeSG(opfdata,bundles)
        end
        bundles[ncuts+1]=mpsoln
      else
	println("Solver returned: $status")
      end
@show kk,mpsoln.linobjval,mpsoln.eta
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end
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

  # MAIN LOOP
    for kk=1:maxNSG
     # STEP 1
      mpsoln=create_bundle(opfdata)
      status = solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,PROX,mpsoln)
      if status == :Optimal || status == :CPX_STAT_NUM_BEST
	comp_agg(opfdata,ctr.soln,mpsoln.soln,sg_agg)
	agg_norm=comp_norm(opfdata,sg_agg)
	epshat=max(0,ctr.eta)-mpsoln.psival-(1.0/tVal)*agg_norm^2
	#del=epshat+0.5*(1.0/tVal)*agg_norm^2
        del=ctr.eta-mpsoln.objval
#@show del,ctr.linobjval,mpsoln.linobjval,mpsoln.eta,mpsoln.psival,agg_norm^2,epshat,mpsoln.status
       # STEP 2
        if del < 1e-6 
	  println("Convergence to within tolerance: del = ",del)
	  break
        end
       # STEP 3
	hk = max(ctr.linobjval-mpsoln.linobjval,mpsoln.eta)
      @show kk,ctr.linobjval,ctr.eta,mpsoln.linobjval,mpsoln.eta,epshat,del,hk
        #if (η_val-ctr.eta) >= -ssc*ctr.eta 
        #if true
#@show hk, mpsoln.psival,ctr.eta, -ssc*del
        if hk <= max(0,ctr.eta) - ssc*del || kk==1
        #if -(mpsoln.psival-ctr.psival)/ctr.psival >= 0.05 
println("Updating center")
	  #purgeSG(opfdata,sg,ctr)
          # UPDATE CENTER VALUES
	    updateCenter(opfdata,mpsoln,ctr)
	    empty!(bundles)
	else
	  ncuts = length(bundles)
          bundles[ncuts+1]=mpsoln
        end
       # STEP 4
        #updateSG(opfdata,sg)
      else
	#println("Solve status: $status")
	#break
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

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
    #x_val[41],x_val[80]=1,1
    x_val[8],x_val[9],x_val[10],x_val[40]=1,1,1,1
    lbs = zeros(opfdata.nlines)
    ubs = ones(opfdata.nlines)
    fixedNode=NodeInfo(x_val,x_val,optUB)
    #fixedNode=NodeInfo(lbs,ubs,optUB)
    # INITIAL SOLN (ALL ZEROS)
      mpsoln=create_bundle(opfdata)
    bestsoln = mpsoln
    ctr=mpsoln 
    bundles = Dict()

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
	      updateCenter(opfdata,bestsoln,ctr)
	    end
	  @show "hk  update",optUB,ctr.linobjval,ctr.eta,hkval,ncuts
        end

     # STEP 3
      cpstatus = solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,CP,cpsoln)
      if cpsoln.status != :Optimal
	println("FLAG! cpsoln.status = ",cpsoln.status)
      end
      if fixedNode.nodeBd-cpsoln.linobjval <= 0.0
        mpsoln=create_bundle(opfdata)
        solveNodeProxPt(opfdata,fixedNode,bundles,K,HEUR,ctr,LVLINF,mpsoln)
	if mpsoln.status != :Optimal
	  println("Taking recourse since mpsoln does not have an optimal solution for MP. Feas: ",fixedNode.nodeBd - mpsoln.linobjval," val: ",fixedNode.nodeBd - cpsoln.linobjval)
	  cpy_soln(opfdata,cpsoln,mpsoln)
	else
	  if hkctr_updated
	    oldncuts = ncuts
	    ncuts=purgeSG(opfdata,bundles)
	    println("There were ",oldncuts," cuts; after purging, there are ",ncuts)
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
	    updateCenter(opfdata,bestsoln,ctr)
	  end
	  @show "bnd update",kk,optUB,bestsoln.linobjval,bestsoln.eta,hkval,ncuts
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end


  #USEFUL SUBROUTINES
    function updateCenter(opfdata,mpsoln,ctr)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      cpy_bundle(opfdata,mpsoln,ctr)
    end
    function computeSG(opfdata,mpsoln)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      vR = zeros(nbuses)
      vI = zeros(nbuses)
      mpsoln.eta = -solveEta0Eigs(opfdata,mpsoln.soln,vR,vI)
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
      E=eigs(H,nev=6,which=:SR, maxiter=100000, tol=1e-8)
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
#=
    function solveEta0SDP(H,opfdata,v)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus
    # DATA RELATED TO THE STORAGE OF VOLTAGE VARIABLE VALUES
      vR = zeros(nbuses)
      vI = zeros(nbuses)

      #The QP subproblem
      mSDP = Model(solver=IpoptSolver())
      @variable(mSDP, e[i=N], start=0); @variable(mSDP, f[i=N], start=0)
      η0Val = 0

      for i in N
        setvalue(e[i], 1); setvalue(f[i], 0)
      end

      # Adjust QP subproblem
      @NLobjective(mSDP, Min, sum( H[i,i]*(e[i]^2+f[i]^2) for i in N)
        + 2*sum( H[fromBus[l],toBus[l]]*(e[fromBus[l]]*e[toBus[l]]+f[fromBus[l]]*f[toBus[l]])   for l in L)
        - 2*sum( H[fromBus[l],nbuses+toBus[l]]*(f[fromBus[l]]*e[toBus[l]]-e[fromBus[l]]*f[toBus[l]])   for l in L)
      )
      status = solve(mSDP)
      if status == :Optimal || status == :UserLimit
        η0Val = getobjectivevalue(mSDP)
        for i in N
          vR[i]=getvalue(e[i]); vI[i]=getvalue(f[i])
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
=#

  # SUBROUTINE FOR COMPUTING A SUBGRADIENT OF ETA(PI), WHICH IS THE FUNCTION TAKING THE VALUE OF THE MINIMUM EIGENVALUE OF H(PI)
    function purgeSG(opfdata,bundle)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      orig_ncuts = length(bundle)
      ncuts = orig_ncuts
      for n=orig_ncuts:-1:1
        if abs(bundle[n].cut_dual) < 1e-8 
	  bundle[n] = bundle[ncuts]
	  delete!(bundle,ncuts)
	  ncuts -= 1
	end
      end
      return ncuts
    end
    function updateSG(opfdata,sg)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      sg.nSGs += 1
      newcut = sg.nSGs
      for i in N
        sg.α[i,newcut] = sg.α_new[i]
	sg.β[i,newcut] = sg.β_new[i]
	sg.γ[i,newcut] = sg.γ_new[i]
        sg.δ[i,newcut] = sg.δ_new[i]
      end
      for l in L
        sg.λF[l,newcut] = sg.λF_new[l]
        sg.λT[l,newcut] = sg.λT_new[l]
        sg.μF[l,newcut] = sg.μF_new[l]
        sg.μT[l,newcut] = sg.μT_new[l]
      end
      sg.trl_soln[newcut]=sg.trl_soln[0]
    end

