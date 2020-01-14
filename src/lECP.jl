#=
Template for branch-and-cut method

July 5, 2018
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")
using CPLEX 
using GLPK 
TOL=1e-4

function solveLECP(opfdata,K,HEUR)

  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = 1:nbuses, 1:nlines, 1:ngens
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # indicate to enable chordal decomposition
    chordal_decomposition = false

    if chordal_decomposition
        # Constrcut a chordal extension of the network topology
        chordal = get_chordal_extension_complex(opfdata)

        # Get the maximal cliques from the chordal extension
        max_cliques = maximal_cliques(get_graph(chordal))
        @show length(max_cliques)

        # Precalculate the number of nodes in each clique
        num_nodes = Dict{Int64,Int64}()
        for k in 1:length(max_cliques)
            clique = sort(max_cliques[k])
            num_nodes[k] = length(clique)
        end
    end


  # The master problem MP
    #mMP = Model(with_optimizer(CPLEX.Optimizer))
    #mMP = direct_model(CPLEX.Optimizer())
    mMP = direct_model(GLPK.Optimizer())
    #MOI.set(mMP, MOI.RawParameter("CPX_PARAM_THREADS"), 1)
#=,
      CPX_PARAM_SCRIND=1,
      CPX_PARAM_TILIM=MAX_TIME,
      CPX_PARAM_MIPDISPLAY=4,
      CPX_PARAM_MIPINTERVAL=1,
      # CPX_PARAM_NODELIM=1,
      # CPX_PARAM_HEURFREQ=-1,
      CPX_PARAM_THREADS=1,
      CPX_PARAM_ADVIND=0))
=#

    #mMP = Model(with_optimizer(Cplex.Optimizer,CPX_PARAM_SCRIND=1,CPX_PARAM_TILIM=MAX_TIME,CPX_PARAM_MIPINTERVAL=50,CPX_PARAM_LPMETHOD=4,CPX_PARAM_SOLUTIONTYPE=2,CPX_PARAM_STARTALG=4))
    # each x[l] is either 0 (line l active) or 1 (line l is cut)
      @variable(mMP, x[l=L], Bin)
    # dual multipliers associated with active power flow balance constraints
      @variable(mMP, -1 <= α[i=N] <= 1)
    # dual multipliers associated with active power flow balance constraints
      @variable(mMP, -1 <= β[i=N] <= 1)
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
      @variable(mMP, λF[l=L]); @variable(mMP, λT[l=L]); @variable(mMP, μF[l=L]); @variable(mMP, μT[l=L])
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
    	+ sum( γ[i]*opfdata.Wmin[i]-δ[i]*opfdata.Wmax[i] + α[i]*opfdata.PD[i] + β[i]*opfdata.QD[i] for i in N) )
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


  # CONTINUING WITH CHORDAL DECOMPOSITION
    @expression(mMP, C[i=1:(2*nbuses),j=i:(2*nbuses)],AffExpr(0.0))
    for i in N
        C[i,i] += δ[i] - γ[i] + α[i]*Y["shR"][i] - β[i]*Y["shI"][i]
        C[nbuses+i,nbuses+i] += δ[i] - γ[i] + α[i]*Y["shR"][i] - β[i]*Y["shI"][i]
        for l in fromLines[i]
            C[i,i] += λF[l]*Y["ffR"][l] - μF[l]*Y["ffI"][l]
            C[nbuses+i,nbuses+i] += λF[l]*Y["ffR"][l] - μF[l]*Y["ffI"][l]
        end
        for l in toLines[i]
            C[i,i] += λT[l]*Y["ttR"][l] - μT[l]*Y["ttI"][l]
            C[nbuses+i,nbuses+i] += λT[l]*Y["ttR"][l] - μT[l]*Y["ttI"][l]
        end
    end
    for l in L
        from=fromBus[l]; to=toBus[l]
        C[from,nbuses+to] -= 0.5*(  λF[l]*Y["ftI"][l] + μF[l]*Y["ftR"][l] - λT[l]*Y["tfI"][l] - μT[l]*Y["tfR"][l]  )
        C[to,nbuses+from] += 0.5*(  λF[l]*Y["ftI"][l] + μF[l]*Y["ftR"][l] - λT[l]*Y["tfI"][l] - μT[l]*Y["tfR"][l]  )
        if from < to
            C[from,to] += 0.5*(λF[l]*Y["ftR"][l]  - μF[l]*Y["ftI"][l] + λT[l]*Y["tfR"][l] - μT[l]*Y["tfI"][l])
            C[nbuses+from,nbuses+to] += 0.5*(λF[l]*Y["ftR"][l]  - μF[l]*Y["ftI"][l] + λT[l]*Y["tfR"][l] - μT[l]*Y["tfI"][l])
        else
            C[to,from] += 0.5*(λF[l]*Y["ftR"][l]  - μF[l]*Y["ftI"][l] + λT[l]*Y["tfR"][l] - μT[l]*Y["tfI"][l])
            C[nbuses+to,nbuses+from] += 0.5*(λF[l]*Y["ftR"][l]  - μF[l]*Y["ftI"][l] + λT[l]*Y["tfR"][l] - μT[l]*Y["tfI"][l])
        end
    end

    if chordal_decomposition
        # Smaller PSD matrices to apply Agler's theorem
        # Each matrix will be approximated by linear inequalities.
        @variable(mMP, subH[k=1:length(max_cliques),i=1:num_nodes[k],j=i:num_nodes[k]])
        @expression(mMP, sumH[i=1:(2*nbuses),j=i:(2*nbuses)], AffExpr(0.0))
        for k in 1:length(max_cliques)
            clique = sort(max_cliques[k])
            # @show (k,clique)
            for i=1:num_nodes[k], j=i:num_nodes[k]
                sumH[clique[i],clique[j]] += subH[k,i,j]
            end
        end
        # Map each smaller PSD matrix to the corresponding element
        @constraint(mMP, [i=1:(2*nbuses),j=i:(2*nbuses)], C[i,j] - sumH[i,j] == 0)
    end
  # END OF CHORDAL DECOMPOSITION

  # DATA RELATED TO SUBGRADIENT INFORMATION
    maxNSG = 1
    sg = Dict()
    sg["α"] = zeros(nbuses,maxNSG)
    sg["β"] = zeros(nbuses,maxNSG)
    sg["γ"] = zeros(nbuses,maxNSG)
    sg["δ"] = zeros(nbuses,maxNSG)
    sg["λF"] = zeros(nlines,maxNSG)
    sg["μF"] = zeros(nlines,maxNSG)
    sg["λT"] = zeros(nlines,maxNSG)
    sg["μT"] = zeros(nlines,maxNSG)

    η0Val = 0

  # DATA RELATED TO THE STORAGE OF X VARIABLE VALUES
    x_val = zeros(nlines)

  # DATA RELATED TO THE STORAGE OF DUAL VARIABLE VALUES
    pi_val = Dict()
    pi_val["α"] = zeros(nbuses)
    pi_val["β"] = zeros(nbuses)
    pi_val["γ"] = zeros(nbuses)
    pi_val["δ"] = zeros(nbuses)
    pi_val["λF"] = zeros(nlines)
    pi_val["λT"] = zeros(nlines)
    pi_val["μF"] = zeros(nlines)
    pi_val["μT"] = zeros(nlines)

  # DATA RELATED TO THE STORAGE OF VOLTAGE VARIABLE VALUES
    v = Dict()
    v["R"] = zeros(nbuses)
    v["I"] = zeros(nbuses)

  # CALLBACK FUNCTION FOR ADDING CUTS
    function addMPCutsLazyCB(cb)
      for l in L
        x_val[l] = round(callback_value(cb,x[l]))
      end
      time_Eta0Start = time_ns()

      if chordal_decomposition
#=
        for k in 1:length(max_cliques)
          clique = sort(max_cliques[k])

          Is = Int64[]; Js = Int64[]; Vs = Float64[]
          for i=1:num_nodes[k], j=i:num_nodes[k]
            # @show (k,i,j,callback_value(cb,subH[k,i,j]))
            if abs(callback_value(cb,subH[k,i,j])) > 1e-10
              push!(Is, i)
              push!(Js, j)
              push!(Vs, callback_value(cb,subH[k,i,j]))
              # Make it symmetric for off-diagonal elements
              if i != j
                push!(Is, j)
                push!(Js, i)
                push!(Vs, callback_value(cb,subH[k,i,j]))
              end
            end
          end
          Hk = sparse(Is,Js,Vs,num_nodes[k],num_nodes[k])
          # @show (typeof(Hk),num_nodes[k],Hk)

          E = eigs(Hk, which=:SR, nev=1)
          # @show E
          η0Val = E[1][1]
          # @show η0Val
          if η0Val < -TOL
            @lazyconstraint(cb,
              sum(E[2][i,1]*E[2][i,1]*subH[k,i,i] for i=1:num_nodes[k])
              + sum(2*E[2][i,1]*E[2][j,1]*subH[k,i,j] for i=1:(num_nodes[k]-1), j=(i+1):num_nodes[k]) >= 0,
                localcut=useLocalCuts
	    )
            ncuts += 1
          end
        end # for loop
=#
    else  # NOT USING CHORDAL DECOMPOSITION
      # generate cut(s)
        for i in N
          pi_val["α"][i] = callback_value(cb,α[i]); pi_val["β"][i] = callback_value(cb,β[i]); pi_val["γ"][i] = callback_value(cb,γ[i]); pi_val["δ"][i] = callback_value(cb,δ[i])
        end
        for l in L
          pi_val["λF"][l] = callback_value(cb,λF[l])
          pi_val["λT"][l] = callback_value(cb,λT[l])
          pi_val["μF"][l] = callback_value(cb,μF[l])
          pi_val["μT"][l] = callback_value(cb,μT[l])
        end
        H=spzeros(2*nbuses,2*nbuses)
        updateHess(opfdata,pi_val,H)

          η0Val = solveEta0Eigs(H,opfdata,v)
        try
        catch exc
          println("Exception caught with eigs(), solving η0Val subproblem with Ipopt as recourse.")
          println(exc)
          η0Val = solveEta0SDP(H,opfdata,v)
        end
        if η0Val < -TOL
          computeSG(opfdata,v,sg)
          con = @build_constraint( sum( (sg["α"][i,1])* (α[i]-pi_val["α"][i]) for i in N) + sum( (sg["β"][i,1])*(β[i]-pi_val["β"][i]) for i in N )
            + sum( (sg["γ"][i,1])*(γ[i]-pi_val["γ"][i]) for i in N)  + sum( (sg["δ"][i,1])*(δ[i]-pi_val["δ"][i]) for i in N)
            + sum( (sg["λF"][l,1])*(λF[l]-pi_val["λF"][l]) + (sg["λT"][l,1])*(λT[l]-pi_val["λT"][l]) for l in L) 
	    + sum( (sg["μF"][l,1])*(μF[l]-pi_val["μF"][l]) + (sg["μT"][l,1])*(μT[l]-pi_val["μT"][l]) for l in L) >= -η0Val
	  )
#println(con)
#println("ncuts: ",ncuts," etaVal: ",η0Val)
#printX(opfdata,x_val)
	  MOI.submit(mMP, MOI.LazyConstraint(cb), con)
#println("submit successful")
          ncuts += 1
        else
          println("Tolerance met for not generating a new lazy cut.")
printX(opfdata,x_val)
        end
    end
    time_Eta0SP += (time_ns()-time_Eta0Start)/1e9
    numcalls_Eta0SP += 1
  end

  MOI.set(mMP, MOI.LazyConstraintCallback(), addMPCutsLazyCB)
  ncuts = 0   # collect number of cuts applied
  mpRootStartTime = time_ns() #For solving root node
  # Collect the timing results
  time_Eta0SP = 0.0
  # some counters
  numcalls_Eta0SP = 0
  nMPUpdates=0
  wtime = @elapsed JuMP.optimize!(mMP)
  status=JuMP.termination_status(mMP)
  @show wtime
  @printf("LAST ATTACK COMPUTED&& ")
  for l in L
    if round(JuMP.value(x[l])) == 1
        println("Cut line $l")
    end
  end
  @printf("\n")
  @printf("\n********************FINAL RESULT FOR CASE %s WITH %d LINES CUT H%dR%d*******************\n",CASE_NUM,K, HEUR,FORM)


  ###Printing summary as one line to aid in latex writeup
  @printf("Summary\t")
  if(HEUR == 0 )
    @printf("%d  &  ",K)
  else
    @printf("    &   ")
  end
  @printf("\$\\HSet^{%d}\$  &  ",HEUR)

  @printf("ECP &    ")


  @printf("%d  & ", ncuts)
  #@printf("%.2f  & ", 100*time_Eta0SP/JuMP.solve_time(mMP))
  @printf("%.2f  & ", 100*time_Eta0SP/JuMP.solve_time(mMP))
  #@printf("%d   &  ", MOI.get(mMP,MOI.NodeCount()))
  #@printf("%.2f   &  ", JuMP.solve_time(mMP)/MOI.get(mMP,MOI.NodeCount()))
  #@printf("%d   &  ", JuMP.solve_time(mMP))
  @printf("%d   &  ", wtime)
 #@printf("%.3f   &  ", JuMP.objective_bound(mMP))

 #@printf("%d    &   ", round(time_root-tShift))
 #@printf("%d   &  ", round(time_Eta0SP))
 #@printf("%d   &  ", round(total_TimeMP))
  for l in L
   if JuMP.value(x[l]) > 0.5
        @printf(" %d",l)
   end
  end
  #@printf("\t&\t%.3f \t&\t%.3f\t & \t%.3f", bestAttack[nlines+1+SDP], bestAttack[nlines+1+SOCP],bestAttack[nlines+1+DC])
 @printf("  & %.3f ", JuMP.objective_value(mMP))
 if abs(JuMP.objective_bound(mMP)-JuMP.objective_value(mMP)) < 1e-3
  @printf("&  0.0\\%% ")
 elseif JuMP.objective_value(mMP) > 1e-3
  percentGap= 100*(JuMP.objective_bound(mMP)-JuMP.objective_value(mMP))/JuMP.objective_value(mMP)
  @printf("&  %.1f\\%% ",percentGap)
 else
  print("& ---     ")
 end
 @printf(" \\\\ \n")
   return mMP
 end

  #USEFUL SUBROUTINES
  # Update Hessian
    function updateHess(opfdata,pi_val,H)
      #lines, buses, generators, baseMVA = opfdata.lines, opfdata.buses, opfdata.generators, opfdata.baseMVA
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, 1:opfdata.nbuses, 1:opfdata.nlines, 1:opfdata.ngens 
      fromBus,toBus = opfdata.fromBus, opfdata.toBus  
      Y = opfdata.Y_AC
      for i in N
        H[i,i] +=  pi_val["α"][i] * Y["shR"][i] - pi_val["β"][i] * Y["shI"][i]  + pi_val["δ"][i] - pi_val["γ"][i]
        H[nbuses+i,nbuses+i] += pi_val["α"][i] * Y["shR"][i] - pi_val["β"][i] * Y["shI"][i] + pi_val["δ"][i] - pi_val["γ"][i]
      end
      for l in L
        from = fromBus[l]; to = toBus[l]
        H[from,from] += pi_val["λF"][l] * Y["ffR"][l] - pi_val["μF"][l] * Y["ffI"][l]
        H[nbuses+from,nbuses+from] += pi_val["λF"][l] * Y["ffR"][l] - pi_val["μF"][l] * Y["ffI"][l]
        H[to,to] += pi_val["λT"][l] * Y["ttR"][l] - pi_val["μT"][l] * Y["ttI"][l]
        H[nbuses+to,nbuses+to] += pi_val["λT"][l] * Y["ttR"][l] - pi_val["μT"][l] * Y["ttI"][l]
        H[from,to] += 0.5*( pi_val["λF"][l] * Y["ftR"][l] - pi_val["μF"][l] * Y["ftI"][l] + pi_val["λT"][l] * Y["tfR"][l] - pi_val["μT"][l] * Y["tfI"][l] )
        H[to,from] += 0.5*( pi_val["λF"][l] * Y["ftR"][l] - pi_val["μF"][l] * Y["ftI"][l] + pi_val["λT"][l] * Y["tfR"][l] - pi_val["μT"][l] * Y["tfI"][l] )
        H[nbuses+from, nbuses+to] += 0.5*( pi_val["λF"][l] * Y["ftR"][l] - pi_val["μF"][l] * Y["ftI"][l] + pi_val["λT"][l] * Y["tfR"][l] - pi_val["μT"][l] * Y["tfI"][l] )
        H[nbuses+to, nbuses+from] += 0.5*( pi_val["λF"][l] * Y["ftR"][l] - pi_val["μF"][l] * Y["ftI"][l] + pi_val["λT"][l] * Y["tfR"][l] - pi_val["μT"][l] * Y["tfI"][l] )
        H[to, nbuses+from] += 0.5*( pi_val["λF"][l] * Y["ftI"][l] - pi_val["λT"][l] * Y["tfI"][l] + pi_val["μF"][l] * Y["ftR"][l] - pi_val["μT"][l] * Y["tfR"][l] )
        H[nbuses+from, to] += 0.5*( pi_val["λF"][l] * Y["ftI"][l] - pi_val["λT"][l] * Y["tfI"][l] + pi_val["μF"][l] * Y["ftR"][l] - pi_val["μT"][l] * Y["tfR"][l] )
        H[from,nbuses+to] -= 0.5*( pi_val["λF"][l] * Y["ftI"][l] - pi_val["λT"][l] * Y["tfI"][l] + pi_val["μF"][l] * Y["ftR"][l] - pi_val["μT"][l] * Y["tfR"][l] )
        H[nbuses+to,from] -= 0.5*( pi_val["λF"][l] * Y["ftI"][l] - pi_val["λT"][l] * Y["tfI"][l] + pi_val["μF"][l] * Y["ftR"][l] - pi_val["μT"][l] * Y["tfR"][l] )
      end
    end

  # SUBROUTINE FOR COMPUTING THE MINIMUM EIGENVALUE OF H WITH A CORRESPONDING EIGENVECTOR
    function solveEta0Eigs(H,optdata,v)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, 1:opfdata.nbuses, 1:opfdata.nlines, 1:opfdata.ngens
      fromBus,toBus = opfdata.fromBus, opfdata.toBus

      E=eigs(H,nev=1,which=:SR, maxiter=100000, tol=1e-4)
      η0Val = E[1][1]
      for i in N
        v["R"][i] = E[2][i,1]; v["I"][i] = E[2][nbuses+i,1]
      end
      return η0Val
    end

  # SUBROUTINE FOR COMPUTING THE MINIMUM EIGENVALUE OF H WITH A CORRESPONDING EIGENVECTOR
    # VIA AN OPTIMIZATION PROBLEM
    function solveEta0SDP(H,opfdata,v)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, 1:opfdata.nbuses, 1:opfdata.nlines, 1:opfdata.ngens
      fromBus,toBus = opfdata.fromBus, opfdata.toBus

      #The QP subproblem
      mSDP = Model(with_optimizer(Ipopt.Optimizer))
      @variable(mSDP, e[i=N], start=1); @variable(mSDP, f[i=N], start=0)
      η0Val = 0


      # Adjust QP subproblem
      @NLobjective(mSDP, Min, sum( H[i,i]*(e[i]^2+f[i]^2) for i in N)
        + 2*sum( H[fromBus[l],toBus[l]]*(e[fromBus[l]]*e[toBus[l]]+f[fromBus[l]]*f[toBus[l]])   for l in L)
        - 2*sum( H[fromBus[l],nbuses+toBus[l]]*(f[fromBus[l]]*e[toBus[l]]-e[fromBus[l]]*f[toBus[l]])   for l in L)
      )
      JuMP.optimize!(mSDP)
      status=JuMP.termination_status(mSDP)
      #if status == :Optimal || status == :UserLimit
      if true
        η0Val = getobjectivevalue(mSDP)
        for i in N
          v["R"][i]=JuMP.value(e[i]); v["I"][i]=JuMP.value(f[i])
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
    function computeSG(opfdata,v,sg)
      nbuses, nlines, ngens, N, L, G = opfdata.nbuses, opfdata.nlines, opfdata.ngens, 1:opfdata.nbuses, 1:opfdata.nlines, 1:opfdata.ngens
      fromBus,toBus,Y = opfdata.fromBus, opfdata.toBus, opfdata.Y_AC

      for i in N
        W_val = v["R"][i]^2 + v["I"][i]^2
        sg["α"][i,1] = Y["shR"][i] * W_val; sg["β"][i,1] = -Y["shI"][i] * W_val
        sg["δ"][i,1] = W_val; sg["γ"][i,1] = -W_val
      end
      for l in L
        from = fromBus[l]; to = toBus[l]
        e_valF = v["R"][from]; f_valF = v["I"][from]; W_valF = e_valF^2 + f_valF^2
        e_valT = v["R"][to]; f_valT = v["I"][to]; W_valT = e_valT^2 + f_valT^2
        Wr_val = e_valF*e_valT + f_valF*f_valT; Wi_val = e_valT*f_valF - e_valF*f_valT
        sg["λF"][l,1] = (Y["ffR"][l] * W_valF + Y["ftR"][l] * Wr_val + Y["ftI"][l] * Wi_val)
        sg["λT"][l,1] = (Y["ttR"][l] * W_valT + Y["tfR"][l] * Wr_val - Y["tfI"][l] * Wi_val)
        sg["μF"][l,1] = (-Y["ffI"][l] * W_valF - Y["ftI"][l] * Wr_val + Y["ftR"][l] * Wi_val)
        sg["μT"][l,1] = (-Y["ttI"][l] * W_valT - Y["tfI"][l] * Wr_val - Y["tfR"][l] * Wi_val)
      end
    end
