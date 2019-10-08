#=
Kibaek Kim
Brian Dandurand
=#

include("utils.jl")
include("MP.jl")

# Implements the approach of Sagastizabal and Solodov 2005
function PBM_SagastizabalSolodov(opfdata,params,K,HEUR,node_data)
    time_Start = time_ns()
  # OBTAIN SHORTHAND PROBLEM INFORMATION FROM opfdata
    nbuses, nlines, ngens = opfdata.nbuses, opfdata.nlines, opfdata.ngens
    N, L, G = opfdata.N, opfdata.L, opfdata.G 
    fromLines,toLines,fromBus,toBus = opfdata.fromLines, opfdata.toLines, opfdata.fromBus, opfdata.toBus
    BusGeners, Y = opfdata.BusGeners, opfdata.Y_AC
  # DONE OBTAINING PROBLEM INFORMATION FROM opfdata

  # INITIAL ITERATION
    ctr_bundles=Dict()
    trl_bundles=Dict()
    agg_bundles=Dict()
    mpsoln=create_bundle(opfdata)
    sg_agg=create_soln(opfdata)
    ctr=mpsoln
    mpsoln=computeMPSoln(opfdata,node_data,K,PROX,ctr,trl_bundles,ctr_bundles,agg_bundles)
    ctr=mpsoln
    node_data.tVal=0.1
  # MAIN LOOP
    for kk=1:params.maxNIter
     # STEP 1
      mpsoln=computeMPSoln(opfdata,node_data,K,PROX,ctr,trl_bundles,ctr_bundles,agg_bundles)
       # STEP 2
        #if node_data.descent_est < 1e-4 
        #if mpsoln.eta < params.tol1 && node_data.agg_sg_norm < params.tol2 && node_data.epshat < params.tol3 
        if mpsoln.eta < params.tol1 && node_data.agg_sg_norm < 1e-5 && node_data.epshat < 1e-5
	  println("Convergence to within tolerance: del = ",node_data.descent_est," val: ",mpsoln.linobjval," and feas: ",mpsoln.eta)
          @show "Sagadov",kk,mpsoln.linobjval,mpsoln.eta,node_data.epshat,node_data.agg_sg_norm
	  break
        end
       # STEP 3
#@show node_data.descent,node_data.descent_est
        if node_data.descent >= params.ssc1*node_data.descent_est || kk==1
          # UPDATE CENTER VALUES
	    ctr=mpsoln
            #agg_bundles[1]=aggregateSG(opfdata,trl_bundles,mpsoln,ctr,ctr_bundles,agg_bundles)
	    #ntrlcuts,nnzcuts=purgeSG(opfdata,trl_bundles,params.maxNSG)
	    ntrlcuts=length(trl_bundles)
	    if ctr.linobjval-mpsoln.linobjval > mpsoln.psival
              trl_bundles[ntrlcuts+1]=mpsoln
	    end
            @show "Sagadov",kk,ctr.linobjval,ctr.eta,node_data.epshat,node_data.agg_sg_norm,node_data.linerr
	else
	  ntrlcuts=length(trl_bundles)
	  #ntrlcuts,nnzcuts=purgeSG(opfdata,trl_bundles,params.maxNSG)
          trl_bundles[ntrlcuts+1]=mpsoln
        end
       # STEP 4
    end
    time_End = (time_ns()-time_Start)/1e9
    println("Done after ",time_End," seconds.")
end

