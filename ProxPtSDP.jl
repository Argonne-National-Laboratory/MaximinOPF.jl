#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")

type NodeInfo
  x_lbs::Array{Float64}
  x_ubs::Array{Float64}
  nodeBd::Float64
end
type SolnInfo
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
  objval::Float64
  linobjval::Float64
  eta::Float64
  solvetime::Float64
end
type EtaSGs
  nSGs::Int
  α::Array{Float64,2}
  β::Array{Float64,2}
  γ::Array{Float64,2}
  δ::Array{Float64,2}
  λF::Array{Float64,2}
  λT::Array{Float64,2}
  μF::Array{Float64,2}
  μT::Array{Float64,2}
  α_new::Array{Float64}
  β_new::Array{Float64}
  γ_new::Array{Float64}
  δ_new::Array{Float64}
  λF_new::Array{Float64}
  λT_new::Array{Float64}
  μF_new::Array{Float64}
  μT_new::Array{Float64}
  trl_soln::Dict
  cut_duals::Array{Float64}
end

maxNSG = 100000
CP,PROX,LVL=0,1,2
tVal = 1
ssc=0.5

function solveNodeProxPt(opfdata,nodeinfo,sg,K,HEUR,ctr,CTR_PARAM,
			mpsoln)
  global tVal
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata


  # The master problem MP
    mMP = Model(solver=CplexSolver(
      CPX_PARAM_SCRIND=0,
      CPX_PARAM_TILIM=MAX_TIME,
      # CPX_PARAM_MIPDISPLAY=4,
      # CPX_PARAM_MIPINTERVAL=1,
      # CPX_PARAM_NODELIM=1,
      # CPX_PARAM_HEURFREQ=-1,
      CPX_PARAM_THREADS=1,
      #CPX_PARAM_ADVIND=0
	)
      )
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
      @variable(mMP, ζpUB[g=G] >=0); @variable(mMP, ζpLB[g=G] >=0); @variable(mMP, ζqUB[g=G] >=0); @variable(mMP, ζqLB[g=G] >=0)

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
    @expression(mMP, linobj, sum(ζpLB[g]*opfdata.Pmin[g] - ζpUB[g]*opfdata.Pmax[g] + ζqLB[g]*opfdata.Qmin[g] - ζqUB[g]*opfdata.Qmax[g]  for g in G)
    			+ sum( γ[i]*opfdata.Wmin[i]-δ[i]*opfdata.Wmax[i] + α[i]*opfdata.PD[i] + β[i]*opfdata.QD[i] for i in N) 
    )
    if CTR_PARAM == PROX || CTR_PARAM == LVL
      @constraint(mMP, LVLConstr, linobj >= nodeinfo.nodeBd)
#=
      @variable(mMP, alphaSlack[i=N] >= 0)
      @variable(mMP, betaSlack[i=N] >= 0)
      @variable(mMP, gammaSlack[i=N] >= 0)
      @variable(mMP, deltaSlack[i=N] >= 0)
      @constraint(mMP, alphaUBSlack[i=N], ctr.α[i] - α[i] <= alphaSlack[i])
      @constraint(mMP, alphaLBSlack[i=N], ctr.α[i] - α[i] >= -alphaSlack[i])
      @constraint(mMP, betaUBSlack[i=N], ctr.β[i] - β[i] <= betaSlack[i])
      @constraint(mMP, betaLBSlack[i=N], ctr.β[i] - β[i] >= -betaSlack[i])
      @constraint(mMP, gammaUBSlack[i=N], ctr.γ[i] - γ[i] <= gammaSlack[i])
      @constraint(mMP, gammaLBSlack[i=N], ctr.γ[i] - γ[i] >= -gammaSlack[i])
      @constraint(mMP, deltaUBSlack[i=N], ctr.δ[i] - δ[i] <= deltaSlack[i])
      @constraint(mMP, deltaLBSlack[i=N], ctr.δ[i] - δ[i] >= -deltaSlack[i])
      @objective(mMP, Max, linobj - tVal*sum(alphaSlack[i] + betaSlack[i] + gammaSlack[i] + deltaSlack[i] for i in N))
      @variable(mMP, Slack >= 0)
      @constraint(mMP, alphaUBSlack[i=N], ctr.α[i] - α[i] <= Slack)
      @constraint(mMP, alphaLBSlack[i=N], ctr.α[i] - α[i] >= -Slack)
      @constraint(mMP, betaUBSlack[i=N], ctr.β[i] - β[i] <= Slack)
      @constraint(mMP, betaLBSlack[i=N], ctr.β[i] - β[i] >= -Slack)
      @constraint(mMP, gammaUBSlack[i=N], ctr.γ[i] - γ[i] <= Slack)
      @constraint(mMP, gammaLBSlack[i=N], ctr.γ[i] - γ[i] >= -Slack)
      @constraint(mMP, deltaUBSlack[i=N], ctr.δ[i] - δ[i] <= Slack)
      @constraint(mMP, deltaLBSlack[i=N], ctr.δ[i] - δ[i] >= -Slack)
      @objective(mMP, Max, linobj - tVal*Slack)
=#
      @objective(mMP, Max, linobj - 0.5*tVal*sum( ( ctr.α[i] - α[i])^2 + (ctr.β[i] - β[i])^2 + (ctr.γ[i] - γ[i])^2 + (ctr.δ[i] - δ[i])^2 for i in N))
    else
      @objective(mMP, Max, linobj)
    end

  # Adding the extra cuts
    if sg.nSGs > 0
      @constraint(mMP, CP[n=1:sg.nSGs], 0 <= sum( sg.α[i,n]*α[i] + sg.β[i,n]*β[i] +  sg.γ[i,n]*γ[i] + sg.δ[i,n]*δ[i] for i in N)
            + sum( sg.λF[l,n]*λF[l] + sg.λT[l,n]*λT[l] + sg.μF[l,n]*μF[l] + sg.μT[l,n]*μT[l] for l in L)
      )
    end
  ### END DEFINING THE LaGRANGIAN DUAL PROBLEM

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

  status=solve(mMP)
  if status == :Optimal || status == :CPX_STAT_NUM_BEST
      for i in N
          mpsoln.α[i],mpsoln.β[i],mpsoln.γ[i],mpsoln.δ[i] = getvalue(α[i]), getvalue(β[i]), getvalue(γ[i]), getvalue(δ[i])
      end
      for g in G
          mpsoln.ζpLB[g],mpsoln.ζpUB[g],mpsoln.ζqLB[g],mpsoln.ζqUB[g] = getvalue(ζpLB[g]), getvalue(ζpUB[g]), getvalue(ζqLB[g]), getvalue(ζqUB[g])
      end
      for l in L
        mpsoln.x[l] = getvalue(x[l])
        mpsoln.λF[l], mpsoln.λT[l], mpsoln.μF[l], mpsoln.μT[l] = getvalue(λF[l]), getvalue(λT[l]), getvalue(μF[l]), getvalue(μT[l])
      end
      mpsoln.objval,mpsoln.solvetime = getobjectivevalue(mMP), getsolvetime(mMP)
      mpsoln.linobjval = getvalue(linobj)
      mpsoln.eta=0
      if status == :Optimal
        for n=1:sg.nSGs
	  sg.cut_duals[n] = getdual(CP[n])
        end
      end
#=
      if CTR_PARAM==LVL || CTR_PARAM==PROX
        @show getdual(LVLConstr)
      end
=#
  else
    #println("solveNodeAC: Return status $status")
  end
  return status
end


function testProxPt(opfdata,K,HEUR)
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # DATA RELATED TO SUBGRADIENT INFORMATION
    sg = EtaSGs(0,
      zeros(nbuses,maxNSG),zeros(nbuses,maxNSG),zeros(nbuses,maxNSG),zeros(nbuses,maxNSG),
      zeros(nlines,maxNSG),zeros(nlines,maxNSG),zeros(nlines,maxNSG),zeros(nlines,maxNSG), 
      zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
      zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines), 
      Dict(),zeros(maxNSG)
    )
  # INITIAL ITERATION
    ctr=SolnInfo(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
			zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
			zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),0.0,0.0,0.0,0.0)
    x_val=zeros(opfdata.nlines)
    #x_val[41],x_val[80]=1,1
    x_val[8],x_val[9],x_val[10],x_val[40]=1,1,1,1
    fixedNode=NodeInfo(x_val,x_val,1e20)
    optUB = 1e20
    mpsoln=SolnInfo(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
	zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
	zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),0.0,0.0,0.0,0.0)
    bestsoln = mpsoln
    status = solveNodeProxPt(opfdata,fixedNode,sg,K,HEUR,ctr,CP,mpsoln)
    if status == :Optimal || status == :CPX_STAT_NUM_BEST
	updateCenter(opfdata,mpsoln,ctr)
        η_val = computeSG(opfdata,mpsoln,sg)
	ctr.eta = η_val
        if η_val < 0
          updateSG(opfdata,sg)
	end
	optUB = mpsoln.objval
	@show mpsoln.objval,-ctr.eta,sg.nSGs
    end
  # MAIN LOOP
    for kk=1:maxNSG
     # STEP 1
      mpsoln=SolnInfo(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
			zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
			zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),0.0,0.0,0.0,0.0)
      fixedNode.nodeBd = optUB + ctr.eta
      status = solveNodeProxPt(opfdata,fixedNode,sg,K,HEUR,ctr,PROX,mpsoln)
      if status == :Optimal
        η_val = computeSG(opfdata,mpsoln,sg)
       # STEP 2
        if -ctr.eta < 1e-6
	  println("Convergence to within tolerance.")
	  break
        end
       # STEP 3
        if (η_val-ctr.eta) >= -ssc*ctr.eta
	  purgeSG(opfdata,sg,ctr)
          # UPDATE CENTER VALUES
	    updateCenter(opfdata,mpsoln,ctr)
	    ctr.eta = η_val
	  @show optUB,mpsoln.linobjval,-ctr.eta,sg.nSGs
        end
       # STEP 4
        updateSG(opfdata,sg)
      else
	optUB = fixedNode.nodeBd 
	@show kk,optUB,-ctr.eta,sg.nSGs
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

  # DATA RELATED TO SUBGRADIENT INFORMATION
    sg = EtaSGs(0,
      zeros(nbuses,maxNSG),zeros(nbuses,maxNSG),zeros(nbuses,maxNSG),zeros(nbuses,maxNSG),
      zeros(nlines,maxNSG),zeros(nlines,maxNSG),zeros(nlines,maxNSG),zeros(nlines,maxNSG), 
      zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
      zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines), 
      Dict(),zeros(maxNSG)
    )
  # INITIAL ITERATION
    optUB=1e20
    ctr=SolnInfo(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
			zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
			zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),0.0,0.0,0.0,0.0)
    x_val=zeros(opfdata.nlines)
    x_val[41],x_val[80]=1,1
    #x_val[8],x_val[9],x_val[10],x_val[40]=1,1,1,1
    lbs = zeros(opfdata.nlines)
    ubs = ones(opfdata.nlines)
    fixedNode=NodeInfo(lbs,ubs,optUB)
    mpsoln=SolnInfo(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
	zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
	zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),0.0,0.0,0.0,0.0)
    bestsoln = mpsoln
    status = solveNodeProxPt(opfdata,fixedNode,sg,K,HEUR,ctr,CP,mpsoln)
    hkval = 0
    hkctr = hkval
    if status == :Optimal 
	updateCenter(opfdata,mpsoln,ctr)
        η_val = computeSG(opfdata,mpsoln,sg)
	ctr.eta = η_val
	hkval = -η_val
        hkctr = hkval
	optUB = mpsoln.objval
        if η_val < 0
          updateSG(opfdata,sg)
	end
	@show mpsoln.objval,optUB,-ctr.eta,hkval,sg.nSGs
        fixedNode.nodeBd=optUB - ssc*hkval
    end
  # MAIN LOOP
    center_updated = false
    for kk=1:maxNSG
     # STEP 1
      hkval = max(optUB - sg.trl_soln[0].linobjval,-sg.trl_soln[0].eta)
      bestsoln = sg.trl_soln[0]
      objval = sg.trl_soln[0].objval
      linobjval = sg.trl_soln[0].linobjval
      etaval = -sg.trl_soln[0].eta
      bestIdx = 0
      for pp=sg.nSGs:-1:1
	if max(optUB - sg.trl_soln[pp].linobjval,-sg.trl_soln[pp].eta) < hkval
	  hkval = max(optUB - sg.trl_soln[pp].linobjval,-sg.trl_soln[pp].eta) 
          bestsoln = sg.trl_soln[pp]
          objval = sg.trl_soln[pp].objval
          linobjval = sg.trl_soln[pp].linobjval
          etaval = -sg.trl_soln[pp].eta
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
      #if -η_val <= -ctr.eta - (1-ssc)*hkval 
        hkctr = hkval
        # UPDATE CENTER VALUES
	  updateCenter(opfdata,bestsoln,ctr)
	  ctr.eta = η_val
	  #purgeSG(opfdata,sg,ctr)
	  if hkval > 1e-4
	    #purgeSG2(opfdata,sg,bestIdx)
	  end
	  @show optUB,linobjval,etaval,hkval,sg.nSGs
          center_updated = true
      end
     # STEP 3
      mpsoln=SolnInfo(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
			zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
			zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),0.0,0.0,0.0,0.0)
      status = solveNodeProxPt(opfdata,fixedNode,sg,K,HEUR,ctr,LVL,mpsoln)
      if status == :Optimal 
        η_val = computeSG(opfdata,mpsoln,sg)
	#purgeSG(opfdata,sg,ctr)
        if η_val < 0
          updateSG(opfdata,sg)
        else
          println("Tolerance met for not generating a new lazy cut, eta=",η_val,".")
	  #@show kk,mpsoln.objval,mpsoln.linobjval,optUB,η_val,hkctr
	  break
	end
	#@show optUB,linobjval,etaval,hkval,η_val,sg.nSGs
	#@show kk,mpsoln.objval,mpsoln.linobjval,optUB,η_val,hkctr,sg.nSGs
        if center_updated
	  purgeSG(opfdata,sg,ctr)
	  center_updated = false
	end
      else
        hkctr = hkval
	optUB = fixedNode.nodeBd 
	#@show kk,optUB,linobjval,etaval,hkval,sg.nSGs
	updateCenter(opfdata,bestsoln,ctr)
	ctr.eta = η_val
	#purgeSG(opfdata,sg,ctr)
	if hkval > 1e-4
	  #purgeSG2(opfdata,sg,bestIdx)
	end
	@show optUB,linobjval,etaval,hkval,sg.nSGs
        center_updated = true
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end



  #USEFUL SUBROUTINES
    function updateCenter(opfdata,mpsoln,ctr)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
        for i in N
          ctr.α[i],ctr.β[i],ctr.γ[i],ctr.δ[i] = mpsoln.α[i],mpsoln.β[i],mpsoln.γ[i],mpsoln.δ[i]
        end 
        for g in G
          ctr.ζpLB[g],ctr.ζpUB[g],ctr.ζqLB[g],ctr.ζqUB[g] = mpsoln.ζpLB[g],mpsoln.ζpUB[g],mpsoln.ζqLB[g],mpsoln.ζqUB[g]
        end
        for l in L
          ctr.x[l] =  mpsoln.x[l] 
          ctr.λF[l], ctr.λT[l], ctr.μF[l], ctr.μT[l] = mpsoln.λF[l], mpsoln.λT[l], mpsoln.μF[l], mpsoln.μT[l] 
        end
        ctr.objval = mpsoln.objval
    end
    function computeSG(opfdata,mpsoln,sg)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      vR = zeros(nbuses)
      vI = zeros(nbuses)
      eta_val = solveEta0Eigs(opfdata,mpsoln,vR,vI)
      for i in N
        W_val = vR[i]^2 + vI[i]^2
        sg.α_new[i] = Y["shR"][i] * W_val
	sg.β_new[i] = -Y["shI"][i] * W_val
        sg.δ_new[i] = W_val
	sg.γ_new[i] = -W_val
      end
      for l in L
        from = fromBus[l]; to = toBus[l]
        e_valF = vR[from]; f_valF = vI[from]; W_valF = e_valF^2 + f_valF^2
        e_valT = vR[to]; f_valT = vI[to]; W_valT = e_valT^2 + f_valT^2
        Wr_val = e_valF*e_valT + f_valF*f_valT; Wi_val = e_valT*f_valF - e_valF*f_valT
        sg.λF_new[l] = (Y["ffR"][l] * W_valF + Y["ftR"][l] * Wr_val + Y["ftI"][l] * Wi_val)
        sg.λT_new[l] = (Y["ttR"][l] * W_valT + Y["tfR"][l] * Wr_val - Y["tfI"][l] * Wi_val)
        sg.μF_new[l] = (-Y["ffI"][l] * W_valF - Y["ftI"][l] * Wr_val + Y["ftR"][l] * Wi_val)
        sg.μT_new[l] = (-Y["ttI"][l] * W_valT - Y["tfI"][l] * Wr_val - Y["tfR"][l] * Wi_val)
      end
      sg.trl_soln[0]=mpsoln
      sg.trl_soln[0].eta = eta_val
      return eta_val
    end
  # SUBROUTINE FOR COMPUTING THE MINIMUM EIGENVALUE OF H WITH A CORRESPONDING EIGENVECTOR
    function solveEta0Eigs(opfdata,mpsoln,vR,vI)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus
      H=spzeros(2*nbuses,2*nbuses)
      updateHess(opfdata,mpsoln,H)
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

  # SUBROUTINE FOR COMPUTING A SUBGRADIENT OF ETA(PI), WHICH IS THE FUNCTION TAKING THE VALUE OF THE MINIMUM EIGENVALUE OF H(PI)
    function purgeSG(opfdata,sg,ctr)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      ncuts = sg.nSGs
      #@show maximum(abs.(sg.cut_duals))
      #sg.cut_duals /= scal
      for n=ncuts:-1:1
	if abs(sg.cut_duals[n]) < 1e-8 
	#if abs(sg.cut_duals[n]) < 1e-6 && (sg.trl_soln[n].eta > ctr.eta)
	#if abs(sg.cut_duals[n]) < 1e-6 && (sg.trl_soln[n].eta > sg.trl_soln[0].eta)
	#if (sg.trl_soln[n].eta > sg.trl_soln[0].eta)
	#if (sg.trl_soln[n].eta > ctr.eta)
	  sg.α[N,n].=sg.α[N,sg.nSGs]
	  sg.β[N,n].=sg.β[N,sg.nSGs]
	  sg.γ[N,n].=sg.γ[N,sg.nSGs]
	  sg.δ[N,n].=sg.δ[N,sg.nSGs]
	  sg.λF[L,n].=sg.λF[L,sg.nSGs]
	  sg.λT[L,n].=sg.λT[L,sg.nSGs]
	  sg.μF[L,n].=sg.μF[L,sg.nSGs]
	  sg.μT[L,n].=sg.μT[L,sg.nSGs]
          sg.trl_soln[n]=sg.trl_soln[sg.nSGs]
	  sg.nSGs -= 1
	end
      end
    end
    function purgeSG2(opfdata,sg,bestIdx)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      if bestIdx == 0
        sg.α[N,1].=sg.α_new[N]
        sg.β[N,1].=sg.β_new[N]
        sg.γ[N,1].=sg.γ_new[N]
        sg.δ[N,1].=sg.δ_new[N]
        sg.λF[L,1].=sg.λF_new[L]
        sg.λT[L,1].=sg.λT_new[L]
        sg.μF[L,1].=sg.μF_new[L]
        sg.μT[L,1].=sg.μT_new[L]
        sg.trl_soln[1]=sg.trl_soln[0]
      else
        sg.α[N,1].=sg.α[N,bestIdx]
        sg.β[N,1].=sg.β[N,bestIdx]
        sg.γ[N,1].=sg.γ[N,bestIdx]
        sg.δ[N,1].=sg.δ[N,bestIdx]
        sg.λF[L,1].=sg.λF[L,bestIdx]
        sg.λT[L,1].=sg.λT[L,bestIdx]
        sg.μF[L,1].=sg.μF[L,bestIdx]
        sg.μT[L,1].=sg.μT[L,bestIdx]
        sg.trl_soln[1]=sg.trl_soln[bestIdx]
      end
      sg.nSGs = 1
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
