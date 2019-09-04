#=
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")
include("MP.jl")

function PBM_DelfinoOliveira(opfdata,params,K,HEUR,node_data)
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
    status = solveNodeMP(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
    while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
      params.tVal /= 2
      println("Status was: ",mpsoln.status,". Resolving with reduced prox parameter value: ",params.tVal)
      status = solveNodeMP(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
    end	
    updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
    ctr_bundles[1]=mpsoln

    agg_bundle=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)

    v_est,ssc_cntr = 1e20,0
    tL,tU=params.tMin,params.tMax
  # MAIN LOOP
    for kk=1:params.maxNSG
     # STEP 1
      mpsoln=create_bundle(opfdata)
      status = solveNodeMP(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
      while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
	params.tVal /= 2
        println("Status was: ",mpsoln.status,". Resolving with reduced prox parameter value: ",params.tVal)
        status = solveNodeMP(opfdata,node_data,params,trl_bundles,ctr_bundles,agg_bundles,K,HEUR,ctr,PROX0,mpsoln)
      end	
      if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
       # STEP 2
        ntrlcuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
	 # UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
         # COMPUTING AGGREGATION INFORMATION
          agg_norm = update_agg(opfdata,params,ctr,mpsoln,sg_agg)
          update_rho(params,trl_bundles,ctr_bundles,agg_bundles)
          epshat = compute_epshat(opfdata,params,mpsoln,ctr,sg_agg)
	  params.rhoUB = params.rho
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
          if ctr.eta < params.tol1 && agg_norm < params.tol2 && epshat < params.tol3 
	    println("Convergence to within tolerance: ")
	    @show kk,ntrlcuts,ssc_cntr,params.tVal,params.rho
	    @show ctr.linobjval,ctr.eta,agg_norm,epshat
	    break
	  else
	    #@show ctr.eta,agg_norm,epshat,params.tVal
          end
         # STEP 3
	  sscval = ((mpsoln.linobjval - params.rhoUB*mpsoln.eta)-(ctr.linobjval - params.rhoUB*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - params.rhoUB*ctr.eta)) 
	  vval = (mpsoln.linobjval-(ctr.linobjval - params.rhoUB*ctr.eta)) 
	  ntrlcuts=purgeSG(opfdata,trl_bundles)
          if sscval >= params.ssc 
         # UPDATE CENTER VALUES
	    nctrcuts=purgeSG(opfdata,ctr_bundles,5,20)
	    ctr_bundles[nctrcuts+1]=mpsoln
	    updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
	    if ctr.eta < params.tol1
	      params.tVal /= 2.0
	    elseif epshat<params.tol3
	      params.tVal /= 1.05
	    end
	    @show kk,nctrcuts,ntrlcuts,ssc_cntr,params.tVal,params.rho
	    @show mpsoln.linobjval,mpsoln.eta,agg_norm,epshat
	  else
            trl_bundles[ntrlcuts+1]=mpsoln
	    if agg_norm-params.tol2 <= epshat-params.tol3
	      params.tVal *= 1.05 
	    end
          end
	  #@show kk,nctrcuts,ntrlcuts,ssc_cntr,params.tVal,params.rho
	  #@show mpsoln.linobjval,mpsoln.eta,agg_norm,epshat
      else
	println("Solver returned: $status")
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

