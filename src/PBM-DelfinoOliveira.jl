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
    plot_opt = zeros(params.maxNIter,5)
    plot_data = zeros(params.maxNIter,7)
    plot_opt_ssteps = zeros(params.maxNIter,6)
    plot_data_ssteps = zeros(params.maxNIter,8)

  # INITIAL ITERATION
    trl_bundles=Dict()
    ctr_bundles=Dict()
    agg_bundles=Dict()
    #initialSG(opfdata,trl_bundles)
    ctr=create_bundle(opfdata)
    mpsoln=computeMPSoln(opfdata,node_data,K,PROX0,ctr,trl_bundles,ctr_bundles,agg_bundles)
    ctr=mpsoln
    ctr_bundles[1]=mpsoln
    sstep_no=1
    #min_it=params.maxNSG
    min_it=50

  # MAIN LOOP
    tLow,tHigh=params.tMin,params.tMax
    for kk=1:params.maxNIter
      node_data.iter=kk
     # STEP 1
      mpsoln=computeMPSoln(opfdata,node_data,K,PROX0,ctr,trl_bundles,ctr_bundles,agg_bundles)
      plot_data[kk,3],plot_data[kk,4],plot_data[kk,5],plot_data[kk,6],plot_data[kk,7]=mpsoln.init_time,mpsoln.solve_time,mpsoln.sg_time,mpsoln.pp_time,mpsoln.bundle_time
      plot_opt[kk,1],plot_opt[kk,2],plot_opt[kk,3],plot_opt[kk,4],plot_opt[kk,5]=mpsoln.linobjval,mpsoln.eta,node_data.linerr,node_data.agg_sg_norm,node_data.epshat
      plot_data_ssteps[sstep_no,3:7] += plot_data[kk,3:7]

     # STEP 2
      ntrlcuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
      node_data.ncuts = ntrlcuts
	# UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
        # COMPUTING AGGREGATION INFORMATION
      if ctr.eta < params.tol1 && node_data.agg_sg_norm < params.tol2 && node_data.epshat < params.tol3 
	println("Convergence to within tolerance: ")
        last_kk=kk
	break
      end
     # STEP 3
      node_data.sscval = ((mpsoln.linobjval - node_data.rhoUB*mpsoln.eta)-(ctr.linobjval - node_data.rhoUB*ctr.eta))/(mpsoln.linobjval-(ctr.linobjval - node_data.rhoUB*ctr.eta)) 
      node_data.descent_est = mpsoln.linobjval-(ctr.linobjval - node_data.rhoUB*ctr.eta) 
      if node_data.sscval >= params.ssc 
        # UPDATE CENTER VALUES
        if testSchrammZoweSSII(opfdata,params,node_data,mpsoln,ctr) || kk < min_it 
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	  ntrlcuts=purgeSG(opfdata,trl_bundles,min(10,params.maxNSG),params.maxNSG)
	  trl_bundles[ntrlcuts+1]=ctr_bundles[1] 	#Move old ctr bundle to the collection of trial bundles
	  ctr_bundles[1]=mpsoln
	  ctr=mpsoln
          tLow,tHigh=params.tMin,params.tMax
	  @printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho)=(%.3f,%.3f)\t(err,||s||,epshat)=(%.2e,%.2e,%.2e)\n",
	      kk,mpsoln.linobjval,mpsoln.eta,node_data.tVal,node_data.rho,node_data.linerr,node_data.agg_sg_norm,node_data.epshat)
          plot_data_ssteps[sstep_no,8] = kk
          plot_opt_ssteps[sstep_no,6] = kk
          if sstep_no > 1
	    plot_data_ssteps[sstep_no,3:7] += plot_data_ssteps[sstep_no-1,3:7]  
          end
	  sstep_no += 1
	  if kk >= min_it
            node_data.tVal = max(0.9*node_data.tVal,params.tMin)
	  end
	else
          tHigh=node_data.tVal
print("Decreasing t from ",node_data.tVal)
	  node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
println("to ",node_data.tVal)
	  
	end
      else
	if testSchrammZoweNSII(opfdata,params,ctr,node_data,mpsoln) || kk < min_it
          agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	  ntrlcuts=purgeSG(opfdata,trl_bundles,min(10,params.maxNSG),params.maxNSG)
          trl_bundles[ntrlcuts+1]=mpsoln
          tLow,tHigh=params.tMin,params.tMax
	else 
          tLow=node_data.tVal
	  node_data.tVal=2*tLow*tHigh/(tLow+tHigh)
	end
      end
    end
    plot_opt[last_kk,1],plot_opt[last_kk,2],plot_opt[last_kk,3],plot_opt[last_kk,4],plot_opt[last_kk,5]=ctr.linobjval,ctr.eta,node_data.linerr,node_data.agg_sg_norm,node_data.epshat
    plot_opt_ssteps[sstep_no,6] = last_kk
    plot_data_ssteps[sstep_no,8] = last_kk
    if sstep_no > 1
      plot_data_ssteps[sstep_no,3:7] += plot_data_ssteps[sstep_no-1,3:7]  
    end
    @printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(t,rho)=(%.3f,%.3f)\t(err,||s||,epshat)=(%.2e,%.2e,%.2e)\n",
      last_kk,ctr.linobjval,ctr.eta,node_data.tVal,node_data.rho,node_data.linerr,node_data.agg_sg_norm,node_data.epshat)
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
    return plot_opt[1:last_kk,:],plot_opt_ssteps[1:sstep_no,:],plot_data[1:last_kk,:],plot_data_ssteps[1:sstep_no,:]
end

