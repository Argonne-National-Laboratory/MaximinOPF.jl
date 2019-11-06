#=
Creating model for scip-sdp

10 Oct 2019
Brian Dandurand
=#

include("utils.jl")

using MathOptFormat

### THIS FUNCTION CREATES A JuMP MODEL, AND CONVERTS IT TO A CBF FORMATTED FILE. 
### NOTE: WHILE THIS FUNCTION GENERATES THE FILE, THE SCIPSDP SOLVER IS PICKY ABOUT THE VERSION AND THE FORMAT OF THE CBF FILE
###    AND SO UNFORTUNATELY, THE FILES GENERATED BY THIS FUNCTION ARE NOT USABLE BY THE SCIPSDP SOLVER.
###    FURTHERMORE, THE IMPLEMENTATION OF MathOptFormat AS OF 11 OCT 2019 FOR CONVERTING JuMP MODELS IS SPOTTY, SEE BELOW NOTES FOR EXAMPLE.
function createSCIPModel(opfdata)
  # OBTAIN PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
    # indicate to enable chordal decomposition
    chordal_decomposition = false

  # Instantiating the model and solver for the dual problem
  # NOTE: UNDER THE CURRENT STATE OF DEVELOPMENT, MathOptFormat REQUIRES MODELS CREATED THROUGH JUMP TO ADD BOUNDS AS CONSTRAINTS
    nThreads=1
    #mMP = Model(solver=MosekSolver(MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=nThreads))
    #mMP = Model(solver=SCSSolver(verbose=1,max_iters=1000000))
    #mMP = Model(with_optimizer(SCS.Optimizer,verbose=1,max_iters=100000,rho_x=1.0))
    mMP = Model(with_optimizer(Mosek.Optimizer,MSK_IPAR_LOG=0,MSK_IPAR_NUM_THREADS=1))
    #mMP = Model(with_optimizer(SCS.Optimizer,verbose=1,max_iters=100000))
    @variable(mMP, α[i=N])
    @constraint(mMP, aLB[i=N], α[i] >= -1)
    @constraint(mMP, aUB[i=N], α[i] <= 1)
    @variable(mMP, β[i=N])
    @constraint(mMP, bLB[i=N], β[i] >= -1)
    @constraint(mMP, bUB[i=N], β[i] <= 1)
    @variable(mMP, γp[i=N] )
    @variable(mMP, γm[i=N] )
    @constraint(mMP, gpLB[i=N], γp[i] >= 0)
    @constraint(mMP, gmLB[i=N], γm[i] >= 0)
    @constraint(mMP, [i=N], γp[i]+γm[i] <= 1)
    @variable(mMP, ζpUB[g=G])
    @constraint(mMP, zpu_nneg[g=G], ζpUB[g] >= 0)
    @variable(mMP, ζpLB[g=G])
    @constraint(mMP, zpl_nneg[g=G], ζpLB[g] >= 0)
    @variable(mMP, ζqUB[g=G])
    @constraint(mMP, zqu_nneg[g=G], ζqUB[g] >= 0)
    @variable(mMP, ζqLB[g=G])
    @constraint(mMP, zql_nneg[g=G], ζqLB[g] >= 0)
  # NOTE: MathOptFormal DOES NOT SEEM TO ACKNOWLEDGE THE Bin TYPE
    @variable(mMP, x[l=L], Int)
    @constraint(mMP, x_lb[l=L], x[l] >= 0)
    @constraint(mMP, x_ub[l=L], x[l] <= 1)
    @constraint(mMP, sum(x[l] for l in L) <= K)


  for i in N
   for g in BusGeners[i]
    @constraint(mMP,
    -α[i] + ζpUB[g] - ζpLB[g] == 0 )
    @constraint(mMP,
    -β[i] + ζqUB[g] - ζqLB[g] == 0 )
   end
  end

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

  @expression(mMP, C[i=1:(2*nbuses),j=i:(2*nbuses)], 0)
  for i in N
    C[i,i] += γp[i] - γm[i] + α[i]*Y["shR"][i] - β[i]*Y["shI"][i]
        C[nbuses+i,nbuses+i] += γp[i] - γm[i] + α[i]*Y["shR"][i] - β[i]*Y["shI"][i]
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
    # Compute new topology for each node subproblem

    # Get maximal cliques from the chordal extension of the network topology
      chordal = get_chordal_extension_complex(opfdata)
      max_cliques = maximal_cliques(get_graph(chordal))

    # Applying Agler's theorem (matrix decomposition)
      subH = Dict{Int64,Any}()
      @expression(mMP, sumH[i=1:(2*nbuses),j=i:(2*nbuses)], 0)
      @show length(max_cliques)
      for k in 1:length(max_cliques)
        clique = sort(max_cliques[k])
        num_nodes = length(clique)
        # @show (k,clique)
        subH[k] = @variable(mMP, [i=1:num_nodes,j=1:num_nodes], PSD)
        for i=1:num_nodes, j=i:num_nodes
          sumH[clique[i],clique[j]] += subH[k][i,j]
        end
      end
      @constraint(mMP, [i=1:(2*nbuses),j=i:(2*nbuses)], C[i,j] - sumH[i,j] == 0)
  else
    # SDP matrix
    @variable(mMP, H[i=1:(2*nbuses),j=1:(2*nbuses)], PSD)
    @constraint(mMP, SetH[i=1:(2*nbuses),j=i:(2*nbuses)], C[i,j] - H[i,j] == 0)
  end

  @objective(mMP, Max, sum(ζpLB[g]*opfdata.Pmin[g] - ζpUB[g]*opfdata.Pmax[g] + ζqLB[g]*opfdata.Qmin[g] - ζqUB[g]*opfdata.Qmax[g]  for g in G)
    + sum( γm[i]*opfdata.Wmin[i]-γp[i]*opfdata.Wmax[i] + α[i]*opfdata.PD[i] + β[i]*opfdata.QD[i] for i in N))



  #These constraints are not quite valid, but their inclusion often results in much faster time to near optimal solution.

  if HEUR == 1
      @constraint(mMP, LambdaMuConstr1[l in L], λF[l]*Y["ftI"][l] - λT[l]*Y["tfI"][l] + μF[l]*Y["ftR"][l] - μT[l]*Y["tfR"][l] == 0.0)
  elseif HEUR == 2
    @constraint(mMP, LambdaFequalsT[l in L], λF[l] - λT[l] == 0)
    @constraint(mMP, muFequalsT[l in L], μF[l] - μT[l] == 0)
  elseif HEUR == 3
    @constraint(mMP, LambdaMuConstr2[l in L], λF[l]*Y["tfR"][l] - λT[l]*Y["ftR"][l] - μF[l]*Y["tfI"][l] + μT[l]*Y["ftI"][l] == 0.0)
  end

  # Now write the model
  
  fn_base=string("CBF/","case",CASE_NUM,"-",K,"-",HEUR,".cbf")
  mathoptformat_model = MathOptFormat.CBF.Model()
  MOI.copy_to(MOI.Bridges.full_bridge_optimizer(mathoptformat_model,Float64), backend(mMP))
  MOI.write_to_file(mathoptformat_model, fn_base)
end #end of function

