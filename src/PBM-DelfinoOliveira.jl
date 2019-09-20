#=
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")
include("MP.jl")

function PBM_DelfinoOliveira(opfdata,params,K,HEUR,node_data)
    println("Applying the algorithm of Delfino and de Oliveira 2018 with the proximal parameter update of Zowe and Schramm 1992....")
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # FOR STORING EXPERIMENTAL DATA
    last_kk=params.maxNIter
    plot_data = zeros(params.maxNIter,7)

  # INITIAL ITERATION
    trl_bundles=Dict()
    ctr_bundles=Dict()
    agg_bundles=Dict()
    ctr=create_bundle(opfdata)
    mpsoln=computeMPSoln(opfdata,node_data,K,PROX0,ctr,trl_bundles,ctr_bundles,agg_bundles)
    ctr=mpsoln
    ctr_bundles[1]=mpsoln


  # MAIN LOOP
    tLow,tHigh=params.tMin,params.tMax
    for kk=1:params.maxNIter
      node_data.iter=kk
     # STEP 1
      mpsoln=computeMPSoln(opfdata,node_data,K,PROX0,ctr,trl_bundles,ctr_bundles,agg_bundles)
      plot_data[kk,1],plot_data[kk,2],plot_data[kk,3],plot_data[kk,4],plot_data[kk,5],plot_data[kk,6],plot_data[kk,7]=mpsoln.linobjval,log(10,mpsoln.eta),mpsoln.init_time,mpsoln.solve_time,mpsoln.sg_time,mpsoln.pp_time,mpsoln.bundle_time

     # STEP 2
      ntrlcuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
      node_data.ncuts = ntrlcuts
	# UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
        # COMPUTING AGGREGATION INFORMATION
      if ctr.eta < params.tol1 && node_data.agg_sg_norm < params.tol2 && node_data.epshat < params.tol3 
	println("Convergence to within tolerance: ")
	@printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho)=(%.3f,%.3f)\t(err,||s||,epshat)=(%.2e,%.2e,%.2e)\n",
	  kk,ctr.linobjval,ctr.eta,node_data.tVal,node_data.rho,node_data.linerr,node_data.agg_sg_norm,node_data.epshat)
        last_kk=kk
        plot_data[kk,1],plot_data[kk,2]=ctr.linobjval,log(10,ctr.eta)
	break
      end
     # STEP 3
      node_data.sscval = ((mpsoln.linobjval - node_data.rhoUB*mpsoln.eta)-(ctr.linobjval - node_data.rhoUB*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - node_data.rhoUB*ctr.eta)) 
      node_data.descent_est = mpsoln.linobjval-(ctr.linobjval - node_data.rhoUB*ctr.eta) 
      if node_data.sscval >= params.ssc 
        # UPDATE CENTER VALUES
        if testSchrammZoweSSII(opfdata,params,node_data,mpsoln,ctr) 
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	  ntrlcuts=purgeSG(opfdata,trl_bundles,10,params.maxNSG)
	  #ntrlcuts=length(trl_bundles)
	  trl_bundles[ntrlcuts+1]=ctr_bundles[1] 	#Move old ctr bundle to the collection of trial bundles
	  ctr_bundles[1]=mpsoln
	  ctr=mpsoln
          tLow,tHigh=params.tMin,params.tMax
	  @printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho)=(%.3f,%.3f)\t(err,||s||,epshat)=(%.2e,%.2e,%.2e)\n",
	      kk,mpsoln.linobjval,mpsoln.eta,node_data.tVal,node_data.rho,node_data.linerr,node_data.agg_sg_norm,node_data.epshat)
          node_data.tVal = max(0.7*node_data.tVal,params.tMin)
	else
          tHigh=node_data.tVal
	  node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
	end
      else
	if testSchrammZoweNSII(opfdata,params,ctr,node_data,mpsoln) 
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	  ntrlcuts=purgeSG(opfdata,trl_bundles,10,params.maxNSG)
	  #ntrlcuts=length(trl_bundles)
          trl_bundles[ntrlcuts+1]=mpsoln
          tLow,tHigh=params.tMin,params.tMax
	else 
          tLow=node_data.tVal
	  node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
	end
      end
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
    return plot_data[1:last_kk,:]
end

