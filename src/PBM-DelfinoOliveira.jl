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
    mpsoln=create_bundle(opfdata)
    ctr=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)
    mMP = createBasicMP(opfdata,node_data,K,PROX0)
    setObjMP(opfdata,mMP,node_data,ctr,PROX0)
    status = solveNodeMP(opfdata,mMP,node_data,trl_bundles,ctr_bundles,agg_bundles,ctr,PROX0,mpsoln,sg_agg)
    while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
      node_data.tVal /= 2
      println("Status was: ",mpsoln.status,". Resolving with reduced prox parameter value: ",node_data.tVal)
      setObjMP(opfdata,mMP,node_data,ctr,PROX0)
      status = solveNodeMP(opfdata,mMP,node_data,trl_bundles,ctr_bundles,agg_bundles,ctr,PROX0,mpsoln,sg_agg)
    end	
    cpy_bundle(opfdata,mpsoln,ctr)
    ctr_bundles[1]=mpsoln

    agg_bundle=create_bundle(opfdata)

    tL,tU=params.tMin,params.tMax
    mMP=nothing
    GC.gc()
  # MAIN LOOP
    for kk=1:params.maxNSG
      node_data.iter=kk
     # STEP 1
      mpsoln=create_bundle(opfdata)
      mMP = createBasicMP(opfdata,node_data,K,PROX0)
      setObjMP(opfdata,mMP,node_data,ctr,PROX0)
      status = solveNodeMP(opfdata,mMP,node_data,trl_bundles,ctr_bundles,agg_bundles,ctr,PROX0,mpsoln,sg_agg)
      while mpsoln.status != MOI.OPTIMAL && mpsoln.status != MOI.LOCALLY_SOLVED
	node_data.tVal /= 2
        println("Status was: ",mpsoln.status,". Resolving with reduced prox parameter value: ",node_data.tVal)
        setObjMP(opfdata,mMP,node_data,ctr,PROX0)
        status = solveNodeMP(opfdata,mMP,node_data,trl_bundles,ctr_bundles,agg_bundles,ctr,PROX0,mpsoln,sg_agg)
      end	
      if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
       # STEP 2
        ntrlcuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
	node_data.ncuts = ntrlcuts
	 # UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
         # COMPUTING AGGREGATION INFORMATION
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
          if ctr.eta < params.tol1 && node_data.agg_sg_norm < params.tol2 && node_data.epshat < params.tol3 
	    println("Convergence to within tolerance: ")
	    @printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho)=(%.3f,%.3f)\t(||s||,epshat)=(%.2e,%.2e)\n",kk,ctr.linobjval,ctr.eta,node_data.tVal,node_data.rho,node_data.agg_sg_norm,node_data.epshat)
	    break
          end
         # STEP 3
	  node_data.sscval = ((mpsoln.linobjval - node_data.rhoUB*mpsoln.eta)-(ctr.linobjval - node_data.rhoUB*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - node_data.rhoUB*ctr.eta)) 
	  node_data.descent_est = mpsoln.linobjval-(ctr.linobjval - node_data.rhoUB*ctr.eta) 
          KiwielRhoUpdate(opfdata,params,node_data)
	  ntrlcuts=purgeSG(opfdata,trl_bundles,20,50)
          if node_data.sscval >= params.ssc 
         # UPDATE CENTER VALUES
	    trl_bundles[ntrlcuts+1]=ctr_bundles[1] #Move old ctr bundle to the collection of trial bundles
	    ctr_bundles[1]=mpsoln
            cpy_bundle(opfdata,mpsoln,ctr)
	    @printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho)=(%.3f,%.3f)\t(||s||,epshat)=(%.2e,%.2e)\n",kk,mpsoln.linobjval,mpsoln.eta,node_data.tVal,node_data.rho,node_data.agg_sg_norm,node_data.epshat)
	  else
            trl_bundles[ntrlcuts+1]=mpsoln
          end
      else
	println("Solver returned: $status")
      end
      mMP=nothing
      GC.gc()
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

