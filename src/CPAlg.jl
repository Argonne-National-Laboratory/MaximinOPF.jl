#=
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")
include("MP.jl")

function CPAlg(opfdata,params,K,HEUR,node_data)
    println("Applying the algorithm of Delfino and de Oliveira 2018 with the proximal parameter update of Zowe and Schramm 1992....")
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # FOR STORING EXPERIMENTAL DATA
    last_kk=params.maxNSG
    plot_data = zeros(params.maxNSG,7)

  # INITIAL ITERATION
    trl_bundles=Dict()
    ctr_bundles=Dict()
    agg_bundles=Dict()
    node_data.tVal=0.0
    ctr=create_bundle(opfdata)
    mpsoln=computeMPSoln(opfdata,node_data,K,CP,ctr,trl_bundles,ctr_bundles,agg_bundles)
    ctr=mpsoln
    trl_bundles[1]=mpsoln

  # MAIN LOOP
    tLow,tHigh=params.tMin,params.tMax
    for kk=1:params.maxNSG
      node_data.iter=kk
     # STEP 1
      mpsoln=computeMPSoln(opfdata,node_data,K,CP,ctr,trl_bundles,ctr_bundles,agg_bundles)
      plot_data[kk,1],plot_data[kk,2],plot_data[kk,3],plot_data[kk,4],plot_data[kk,5],plot_data[kk,6],plot_data[kk,7]=mpsoln.linobjval,log(10,mpsoln.eta),mpsoln.init_time,mpsoln.solve_time,mpsoln.sg_time,mpsoln.pp_time,mpsoln.bundle_time

     # STEP 2
      ntrlcuts,nctrcuts,naggcuts = length(trl_bundles),length(ctr_bundles),length(agg_bundles)
      node_data.ncuts = ntrlcuts
	# UPDATING RHO AS NECESSARY TO CORRESPOND TO EXACT PENALTY
        # COMPUTING AGGREGATION INFORMATION
      if max(node_data.rho*mpsoln.eta,mpsoln.eta) < params.tol1
	println("Convergence to within tolerance: ")
	@printf("iter: %d\t(objval,eta)=(%.4f,%.2e)\t(rho)=(%.3f)\n",
	  kk,ctr.linobjval,ctr.eta,node_data.rho)
        last_kk=kk
	break
      end
     # STEP 3
      #ntrlcuts=purgeSG(opfdata,trl_bundles,10,100000)
      ntrlcuts=length(trl_bundles)
      trl_bundles[ntrlcuts+1]=mpsoln 
      @printf("iter: %d\t(objval,eta,rho,rho*eta)=(%.4f,%.2e,%.3f,%.3e)\n",
	      kk,mpsoln.linobjval,mpsoln.eta,node_data.rho,node_data.rho*mpsoln.eta)
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
    return plot_data[1:last_kk,:]
end

