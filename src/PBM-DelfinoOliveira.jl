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
    TOL = 1e-5
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
	    @show ctr.linobjval,ctr.eta,agg_norm,epshat
	    break
          end
         # STEP 3
	  sscval = ((mpsoln.linobjval - params.rhoUB*mpsoln.eta)-(ctr.linobjval - params.rhoUB*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - params.rhoUB*ctr.eta)) 
	  vval = (mpsoln.linobjval-(ctr.linobjval - params.rhoUB*ctr.eta)) 
          if sscval >= params.ssc 
         # UPDATE CENTER VALUES
	    nctrcuts=purgeSG(opfdata,ctr_bundles)
	    ctr_bundles[nctrcuts+1]=mpsoln
	    updateCenter(opfdata,mpsoln,ctr,trl_bundles,ctr_bundles,agg_bundles)
	    if mpsoln.eta < 1e-5
	      params.tVal *= 0.8
	    end
	    @show kk,ncuts,ssc_cntr,params.tVal,params.rho
	    @show mpsoln.linobjval,mpsoln.eta,agg_norm,epshat
	  else
	    ncuts=purgeSG(opfdata,trl_bundles)
            trl_bundles[ncuts+1]=mpsoln
          end
      else
	println("Solver returned: $status")
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

