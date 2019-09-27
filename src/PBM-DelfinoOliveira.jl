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
    plot_params = zeros(params.maxNIter,3)
    plot_opt = zeros(params.maxNIter,5)
    plot_data = zeros(params.maxNIter,7)
    plot_params_ssteps = zeros(params.maxNIter,4)
    plot_opt_ssteps = zeros(params.maxNIter,6)
    plot_data_ssteps = zeros(params.maxNIter,8)

  # INITIAL ITERATION
    trl_bundles=Dict()
    ctr_bundles=Dict()
    agg_bundles=Dict()
    ctr=create_bundle(opfdata)
    sstep_no=1

  # MAIN LOOP
    tLow,tHigh=params.tMin,params.tMax
    for kk=1:params.maxNIter
      node_data.iter=kk
     # STEP 1
      mpsoln=computeMPSoln(opfdata,node_data,K,PROX0,ctr,trl_bundles,ctr_bundles,agg_bundles)
      rho_est=node_data.rho
      node_data.sscval = ((mpsoln.linobjval - rho_est*mpsoln.eta)-(ctr.linobjval - rho_est*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - rho_est*ctr.eta)) 
      node_data.descent_est = mpsoln.linobjval-(ctr.linobjval - node_data.rho*ctr.eta) 
      plot_data[kk,3],plot_data[kk,4],plot_data[kk,5],plot_data[kk,6],plot_data[kk,7]=mpsoln.init_time,mpsoln.solve_time,mpsoln.sg_time,mpsoln.pp_time,mpsoln.bundle_time
      plot_opt[kk,1],plot_opt[kk,2],plot_opt[kk,3],plot_opt[kk,4],plot_opt[kk,5]=mpsoln.linobjval,mpsoln.eta,node_data.linerr,node_data.agg_sg_norm,node_data.epshat
      plot_data_ssteps[sstep_no,3:7] += plot_data[kk,3:7]

     # STEP 2
      ntrlcuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
      node_data.ncuts = ntrlcuts+nctrcuts+naggcuts
      plot_params[kk,1],plot_params[kk,2],plot_params[kk,3]=node_data.tVal,node_data.rho,node_data.ncuts
	# UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
        # COMPUTING AGGREGATION INFORMATION
      if mpsoln.eta < params.tol1 && node_data.agg_sg_norm < params.tol2 && node_data.epshat < params.tol3 
        ctr=mpsoln
	println("Convergence to within tolerance: ")
        last_kk=kk
	break
      end
     # STEP 3
      if node_data.sscval >= params.ssc1 
        # UPDATE CENTER VALUES
        if testSchrammZoweSSII(opfdata,params,node_data,mpsoln,ctr) || tHigh-tLow < 1e-2 
          plot_opt_ssteps[sstep_no,1],plot_opt_ssteps[sstep_no,2],plot_opt_ssteps[sstep_no,3],plot_opt_ssteps[sstep_no,4],plot_opt_ssteps[sstep_no,5]=plot_opt[kk,1],plot_opt[kk,2],plot_opt[kk,3],plot_opt[kk,4],plot_opt[kk,5]
          plot_params_ssteps[sstep_no,1],plot_params_ssteps[sstep_no,2],plot_params_ssteps[sstep_no,3]=plot_params[kk,1],plot_params[kk,2],plot_params[kk,3]
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	  ntrlcuts=purgeSG(opfdata,trl_bundles,params.maxNSG)
	  for n=1:length(ctr_bundles)
	    trl_bundles[ntrlcuts+n]=ctr_bundles[n] 	#Move old ctr bundle to the collection of trial bundles
	  end
	  ctr_bundles[1]=mpsoln
	  ctr=mpsoln
          tLow,tHigh=params.tMin,params.tMax
	  @printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho,ncuts)=(%.3f,%.3f,%d)\t(err,||s||,epshat,desc)=(%.2e,%.2e,%.2e,%.5e)\n",
	      kk,mpsoln.linobjval,mpsoln.eta,node_data.tVal,node_data.rho,node_data.ncuts,node_data.linerr,node_data.agg_sg_norm,node_data.epshat,node_data.descent_est)
	  plot_params_ssteps[sstep_no,4] = kk
          plot_data_ssteps[sstep_no,8] = kk
          plot_opt_ssteps[sstep_no,6] = kk
	  sstep_no += 1
          node_data.tVal = max(0.95*node_data.tVal,params.tMin)
	else
          tHigh=node_data.tVal
	  node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
	end
      else
	if testSchrammZoweNSII(opfdata,params,ctr,node_data,mpsoln) || tHigh-tLow < 1e-2
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	  ntrlcuts=purgeSG(opfdata,trl_bundles,params.maxNSG)
          trl_bundles[ntrlcuts+1]=mpsoln
          tLow,tHigh=params.tMin,params.tMax
	else 
          tLow=node_data.tVal
	  node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
	end
      end
    end
    plot_opt[last_kk,1],plot_opt[last_kk,2],plot_opt[last_kk,3],plot_opt[last_kk,4],plot_opt[last_kk,5]=ctr.linobjval,ctr.eta,node_data.linerr,node_data.agg_sg_norm,node_data.epshat
    plot_params[last_kk,1],plot_params[last_kk,2],plot_params[last_kk,3]=node_data.tVal,node_data.rho,node_data.ncuts
    plot_params_ssteps[sstep_no,4] = last_kk
    plot_opt_ssteps[sstep_no,6] = last_kk
    plot_data_ssteps[sstep_no,8] = last_kk
    @printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho,ncuts)=(%.3f,%.3f,%d)\t(err,||s||,epshat,desc)=(%.2e,%.2e,%.2e,%.5e)\n",
      last_kk,ctr.linobjval,ctr.eta,node_data.tVal,node_data.rho,node_data.ncuts,node_data.linerr,node_data.agg_sg_norm,node_data.epshat,node_data.descent_est)
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
    printX2(opfdata,ctr.soln.x)
    return plot_params[1:last_kk,:],plot_opt[1:last_kk,:],plot_data[1:last_kk,:],plot_params_ssteps[1:sstep_no,:],plot_opt_ssteps[1:sstep_no,:],plot_data_ssteps[1:sstep_no,:]
end

