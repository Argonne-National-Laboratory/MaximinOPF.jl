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

  # INITIAL ITERATION
    trl_bundles=Dict()
    ctr_bundles=Dict()
    agg_bundles=Dict()
    #initialSG(opfdata,trl_bundles)
    ctr=create_bundle(opfdata)
    mpsoln=computeMPSoln(opfdata,node_data,K,PROX0,ctr,trl_bundles,ctr_bundles,agg_bundles)
    ctr_bundles[1]=mpsoln
    ctr=mpsoln


  # MAIN LOOP
    tLow,tHigh=params.tMin,params.tMax
    for kk=1:params.maxNSG
      node_data.iter=kk
     # STEP 1
      mpsoln=computeMPSoln(opfdata,node_data,K,PROX0,ctr,trl_bundles,ctr_bundles,agg_bundles)
      node_data.sscval = ((mpsoln.linobjval - node_data.rhoUB*mpsoln.eta)-(ctr.linobjval - node_data.rhoUB*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - node_data.rhoUB*ctr.eta)) 
      node_data.descent_est = mpsoln.linobjval-(ctr.linobjval - node_data.rhoUB*ctr.eta) 

     # STEP 2
      ntrlcuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
      node_data.ncuts = ntrlcuts
	# UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
        # COMPUTING AGGREGATION INFORMATION
      if mpsoln.eta < params.tol1 && node_data.agg_sg_norm < params.tol2 && node_data.epshat < params.tol3 
        ctr=mpsoln
	println("Convergence to within tolerance: ")
	@printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho)=(%.3f,%.3f)\t(err,||s||,epshat)=(%.2e,%.2e,%.2e)\n",
	  kk,ctr.linobjval,ctr.eta,node_data.tVal,node_data.rho,node_data.linerr,node_data.agg_sg_norm,node_data.epshat)
	break
      elseif mpsoln.eta < 1e-8 && abs(tHigh-tLow) > 1e-2
        tHigh=node_data.tVal
	node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
	@printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho)=(%.3f,%.3f)\t(err,||s||,epshat,desc_est)=(%.2e,%.2e,%.2e,%.5e)\n",
	  kk,mpsoln.linobjval,mpsoln.eta,node_data.tVal,node_data.rho,node_data.linerr,node_data.agg_sg_norm,node_data.epshat,node_data.descent_est)
	continue
      end
     # STEP 3
      if node_data.sscval >= params.ssc1
        # UPDATE CENTER VALUES
        if testSchrammZoweSSII(opfdata,params,node_data,mpsoln,ctr) || abs(tHigh-tLow) <= 1e-2
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	  ntrlcuts=purgeSG(opfdata,trl_bundles,10,100)
	  for n=1:length(ctr_bundles) 
	    trl_bundles[ntrlcuts+n]=ctr_bundles[n] 	#Move old ctr bundle to the collection of trial bundles
	    delete!(ctr_bundles,n)
	  end
	  ctr_bundles[1]=mpsoln
	  ctr=mpsoln
          tLow,tHigh=params.tMin,params.tMax
	  @printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho)=(%.3f,%.3f)\t(err,||s||,epshat,desc_est)=(%.2e,%.2e,%.2e,%.5e)\n",
	      kk,mpsoln.linobjval,mpsoln.eta,node_data.tVal,node_data.rho,node_data.linerr,node_data.agg_sg_norm,node_data.epshat,node_data.descent_est)
	  node_data.tVal = max(0.5*node_data.tVal,params.tMin)
	else
          tHigh=node_data.tVal
	  node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
	end
      else
	if testSchrammZoweNSII(opfdata,params,ctr,node_data,mpsoln) || abs(tHigh-tLow) <= 1e-2
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	  ntrlcuts=purgeSG(opfdata,trl_bundles,10,100)
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
end

