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
  α_val::Array{Float64}
  β_val::Array{Float64}
  γ_val::Array{Float64}
  δ_val::Array{Float64}
  ζpLB_val::Array{Float64}
  ζpUB_val::Array{Float64}
  ζqLB_val::Array{Float64}
  ζqUB_val::Array{Float64}
  x_val::Array{Float64}
  λF_val::Array{Float64}
  λT_val::Array{Float64}
  μF_val::Array{Float64}
  μT_val::Array{Float64}
  objval::Float64
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
  cut_duals::Array{Float64}
end

maxNSG = 3000

function solveNodeProxPt(opfdata,nodeinfo,sg,K,HEUR,mpsoln, ctr)

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
    @objective(mMP, Max, sum(ζpLB[g]*opfdata.Pmin[g] - ζpUB[g]*opfdata.Pmax[g] + ζqLB[g]*opfdata.Qmin[g] - ζqUB[g]*opfdata.Qmax[g]  for g in G)
    	+ sum( γ[i]*opfdata.Wmin[i]-δ[i]*opfdata.Wmax[i] + α[i]*opfdata.PD[i] + β[i]*opfdata.QD[i] for i in N) 
	- 0.5*sum( ( ctr.α_val[i] - α[i])^2 + (ctr.β_val[i] - β[i])^2 + (ctr.γ_val[i] - γ[i])^2 + (ctr.δ_val[i] - δ[i])^2 for i in N)
    )

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
  if status == :Optimal 
      for i in N
          mpsoln.α_val[i],mpsoln.β_val[i],mpsoln.γ_val[i],mpsoln.δ_val[i] = getvalue(α[i]), getvalue(β[i]), getvalue(γ[i]), getvalue(δ[i])
      end
      for g in G
          mpsoln.ζpLB_val[g],mpsoln.ζpUB_val[g],mpsoln.ζqLB_val[g],mpsoln.ζqUB_val[g] = getvalue(ζpLB[g]), getvalue(ζpUB[g]), getvalue(ζqLB[g]), getvalue(ζqUB[g])
      end
      for l in L
        mpsoln.x_val[l] = getvalue(x[l])
        mpsoln.λF_val[l], mpsoln.λT_val[l], mpsoln.μF_val[l], mpsoln.μT_val[l] = getvalue(λF[l]), getvalue(λT[l]), getvalue(μF[l]), getvalue(μT[l])
      end
      mpsoln.objval,mpsoln.solvetime = getobjectivevalue(mMP), getsolvetime(mMP)
      for n=1:sg.nSGs
	sg.cut_duals[n] = getdual(CP[n])
      end
  else
    println("solveNodeAC: Return status $status")
  end
end



function testProxPt(opfdata,K,HEUR)
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # DATA RELATED TO SUBGRADIENT INFORMATION
    sg = EtaSGs(0,
      zeros(nbuses,maxNSG),zeros(nbuses,maxNSG),zeros(nbuses,maxNSG),zeros(nbuses,maxNSG),
      zeros(nlines,maxNSG),zeros(nlines,maxNSG),zeros(nlines,maxNSG),zeros(nlines,maxNSG), zeros(maxNSG)
    )
  # DATA RELATED TO THE STORAGE OF DUAL VARIABLE VALUES FOR THE MP
    mpsoln=SolnInfo(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
			zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
			zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),0.0,0.0)
    ctr=SolnInfo(zeros(nbuses),zeros(nbuses),zeros(nbuses),zeros(nbuses),
			zeros(ngens),zeros(ngens),zeros(ngens),zeros(ngens),
			zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),zeros(nlines),0.0,0.0)
    x_val=zeros(opfdata.nlines)
    x_val[41],x_val[80]=1,1
    fixedNode=NodeInfo(x_val,x_val,1e20)
  # DATA RELATED TO THE STORAGE OF VOLTAGE VARIABLE VALUES
    v = Dict()
    v["R"] = zeros(nbuses)
    v["I"] = zeros(nbuses)
    for kk=1:maxNSG
      solveNodeProxPt(opfdata,fixedNode,sg,K,HEUR,mpsoln,ctr)
@show mpsoln.objval
      H=spzeros(2*nbuses,2*nbuses)
      updateHess(opfdata,mpsoln,H)
      η_val = solveEta0Eigs(H,opfdata,v)
      if η_val < -TOL
        updateSG(opfdata,v,sg)
      else
        println("Tolerance met for not generating a new lazy cut.")
      end

  # UPDATE CENTER VALUES
      if kk==1 || (ctr.objval-η_val)/(ctr.objval) >= 0.1
        for i in N
          ctr.α_val[i],ctr.β_val[i],ctr.γ_val[i],ctr.δ_val[i] = mpsoln.α_val[i],mpsoln.β_val[i],mpsoln.γ_val[i],mpsoln.δ_val[i]
        end 
        for g in G
          ctr.ζpLB_val[g],ctr.ζpUB_val[g],ctr.ζqLB_val[g],ctr.ζqUB_val[g] = mpsoln.ζpLB_val[g],mpsoln.ζpUB_val[g],mpsoln.ζqLB_val[g],mpsoln.ζqUB_val[g]
        end
        for l in L
          ctr.x_val[l] =  mpsoln.x_val[l] 
          ctr.λF_val[l], ctr.λT_val[l], ctr.μF_val[l], ctr.μT_val[l] = mpsoln.λF_val[l], mpsoln.λT_val[l], mpsoln.μF_val[l], mpsoln.μT_val[l] 
        end
        ctr.objval = η_val
      end
    end
end



  #USEFUL SUBROUTINES
  # Update Hessian
    function updateHess(opfdata,pi_val,H)
      #lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus  
      Y = opfdata.Y_AC
      for i in N
        H[i,i] +=  pi_val.α_val[i] * Y["shR"][i] - pi_val.β_val[i] * Y["shI"][i]  + pi_val.δ_val[i] - pi_val.γ_val[i]
        H[nbuses+i,nbuses+i] += pi_val.α_val[i] * Y["shR"][i] - pi_val.β_val[i] * Y["shI"][i] + pi_val.δ_val[i] - pi_val.γ_val[i]
      end
      for l in L
        from = fromBus[l]; to = toBus[l]
        H[from,from] += pi_val.λF_val[l] * Y["ffR"][l] - pi_val.μF_val[l] * Y["ffI"][l]
        H[nbuses+from,nbuses+from] += pi_val.λF_val[l] * Y["ffR"][l] - pi_val.μF_val[l] * Y["ffI"][l]
        H[to,to] += pi_val.λT_val[l] * Y["ttR"][l] - pi_val.μT_val[l] * Y["ttI"][l]
        H[nbuses+to,nbuses+to] += pi_val.λT_val[l] * Y["ttR"][l] - pi_val.μT_val[l] * Y["ttI"][l]
        H[from,to] += 0.5*( pi_val.λF_val[l] * Y["ftR"][l] - pi_val.μF_val[l] * Y["ftI"][l] + pi_val.λT_val[l] * Y["tfR"][l] - pi_val.μT_val[l] * Y["tfI"][l] )
        H[to,from] += 0.5*( pi_val.λF_val[l] * Y["ftR"][l] - pi_val.μF_val[l] * Y["ftI"][l] + pi_val.λT_val[l] * Y["tfR"][l] - pi_val.μT_val[l] * Y["tfI"][l] )
        H[nbuses+from, nbuses+to] += 0.5*( pi_val.λF_val[l] * Y["ftR"][l] - pi_val.μF_val[l] * Y["ftI"][l] + pi_val.λT_val[l] * Y["tfR"][l] - pi_val.μT_val[l] * Y["tfI"][l] )
        H[nbuses+to, nbuses+from] += 0.5*( pi_val.λF_val[l] * Y["ftR"][l] - pi_val.μF_val[l] * Y["ftI"][l] + pi_val.λT_val[l] * Y["tfR"][l] - pi_val.μT_val[l] * Y["tfI"][l] )
        H[to, nbuses+from] += 0.5*( pi_val.λF_val[l] * Y["ftI"][l] - pi_val.λT_val[l] * Y["tfI"][l] + pi_val.μF_val[l] * Y["ftR"][l] - pi_val.μT_val[l] * Y["tfR"][l] )
        H[nbuses+from, to] += 0.5*( pi_val.λF_val[l] * Y["ftI"][l] - pi_val.λT_val[l] * Y["tfI"][l] + pi_val.μF_val[l] * Y["ftR"][l] - pi_val.μT_val[l] * Y["tfR"][l] )
        H[from,nbuses+to] -= 0.5*( pi_val.λF_val[l] * Y["ftI"][l] - pi_val.λT_val[l] * Y["tfI"][l] + pi_val.μF_val[l] * Y["ftR"][l] - pi_val.μT_val[l] * Y["tfR"][l] )
        H[nbuses+to,from] -= 0.5*( pi_val.λF_val[l] * Y["ftI"][l] - pi_val.λT_val[l] * Y["tfI"][l] + pi_val.μF_val[l] * Y["ftR"][l] - pi_val.μT_val[l] * Y["tfR"][l] )
      end
    end

  # SUBROUTINE FOR COMPUTING THE MINIMUM EIGENVALUE OF H WITH A CORRESPONDING EIGENVECTOR
    function solveEta0Eigs(H,optdata,v)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus

      E=eigs(H,nev=6,which=:SR, maxiter=100000, tol=1e-8)
      η0Val = E[1][1]
      for i in N
        v["R"][i] = E[2][i,1]; v["I"][i] = E[2][nbuses+i,1]
      end
      return η0Val
    end

  # SUBROUTINE FOR COMPUTING THE MINIMUM EIGENVALUE OF H WITH A CORRESPONDING EIGENVECTOR
    # VIA AN OPTIMIZATION PROBLEM
    function solveEta0SDP(H,opfdata,v)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus

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
          v["R"][i]=getvalue(e[i]); v["I"][i]=getvalue(f[i])
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
    function updateSG(opfdata,v,sg)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, opfdata.N, opfdata.L, opfdata.G 
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC
      ncuts = sg.nSGs
      for n=ncuts:-1:1
	if sg.cut_duals[n] < 1e-6
	  sg.α[N,n].=sg.α[N,sg.nSGs]
	  sg.β[N,n].=sg.β[N,sg.nSGs]
	  sg.γ[N,n].=sg.γ[N,sg.nSGs]
	  sg.δ[N,n].=sg.δ[N,sg.nSGs]
	  sg.λF[L,n].=sg.λF[L,sg.nSGs]
	  sg.λT[L,n].=sg.λT[L,sg.nSGs]
	  sg.μF[L,n].=sg.μF[L,sg.nSGs]
	  sg.μT[L,n].=sg.μT[L,sg.nSGs]
	  sg.nSGs -= 1
	end
      end
      sg.nSGs += 1
      newcut = sg.nSGs
      for i in N
        W_val = v["R"][i]^2 + v["I"][i]^2
        sg.α[i,newcut] = Y["shR"][i] * W_val; sg.β[i,newcut] = -Y["shI"][i] * W_val
        sg.δ[i,newcut] = W_val; sg.γ[i,newcut] = -W_val
      end
      for l in L
        from = fromBus[l]; to = toBus[l]
        e_valF = v["R"][from]; f_valF = v["I"][from]; W_valF = e_valF^2 + f_valF^2
        e_valT = v["R"][to]; f_valT = v["I"][to]; W_valT = e_valT^2 + f_valT^2
        Wr_val = e_valF*e_valT + f_valF*f_valT; Wi_val = e_valT*f_valF - e_valF*f_valT
        sg.λF[l,newcut] = (Y["ffR"][l] * W_valF + Y["ftR"][l] * Wr_val + Y["ftI"][l] * Wi_val)
        sg.λT[l,newcut] = (Y["ttR"][l] * W_valT + Y["tfR"][l] * Wr_val - Y["tfI"][l] * Wi_val)
        sg.μF[l,newcut] = (-Y["ffI"][l] * W_valF - Y["ftI"][l] * Wr_val + Y["ftR"][l] * Wi_val)
        sg.μT[l,newcut] = (-Y["ttI"][l] * W_valT - Y["tfI"][l] * Wr_val - Y["tfR"][l] * Wi_val)
      end
    end
