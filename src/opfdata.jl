using DelimitedFiles
using Printf
using PowerModels

#=

mutable struct Gener
  # .gen fields
  bus::Int
  Pg::Float64
  Qg::Float64
  Qmax::Float64
  Qmin::Float64
  Vg::Float64
  mBase::Float64
  status::Int
  Pmax::Float64
  Pmin::Float64
  Pc1::Float64
  Pc2::Float64
  Qc1min::Float64
  Qc1max::Float64
  Qc2min::Float64
  Qc2max::Float64
  ramp_agc::Float64
  # .gencost fields
  gentype::Int
  startup::Float64
  shutdown::Float64
  n::Int
  coeff::Array
end
=#

mutable struct OPFData
  nbuses::Int
  PD::Array
  QD::Array
  Wmin::Array
  Wmax::Array
  nlines::Int
  fromLines::Array         #From lines for each bus (Array of Array)
  toLines::Array           #To lines for each bus (Array of Array)
  fromBus::Array
  toBus::Array
  ngens::Int
  Pmin::Array{Float64}
  Pmax::Array{Float64}
  Qmin::Array{Float64}
  Qmax::Array{Float64}
  BusGeners
  bus_ref::Int
  Y_AC
  Y_DC
end


function opf_loaddata(case_name)

  data_path = "data/case"
  fn_base=string(data_path,case_name)
  fn_mat=string(fn_base,".m")
  pm_data = PowerModels.parse_file(fn_mat)
  #@show pm_data

#=

  #
  # load generators
  #
  # gen_arr = readdlm("data/" * case_name * ".gen")
  gen_arr = readdlm(data_path * string(case_name) * ".gen")
  # costgen_arr = readdlm("data/" * case_name * ".gencost")
  costgen_arr = readdlm(data_path * string(case_name) * ".gencost")
  num_gens = size(gen_arr,1)


  @assert num_gens == size(costgen_arr,1)

  gens_on=findall(x->x!=0, gen_arr[:,8]); num_on=length(gens_on)
  if num_gens-num_on>0
    println("loaddata: ", num_gens-num_on, " generators are off and will be discarded (out of ", num_gens, ")")
  end

  generators = Array{Gener}(undef,num_on)
  i=0
  for git in gens_on
    i += 1
    generators[i] = Gener(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, Array{Int}(undef,0)) #gen_arr[i,1:end]...)

    generators[i].bus      = gen_arr[git,1]
    generators[i].Pg       = gen_arr[git,2] / pm_data["baseMVA"]
    generators[i].Qg       = gen_arr[git,3] / pm_data["baseMVA"]
    generators[i].Qmax     = gen_arr[git,4] / pm_data["baseMVA"]
    generators[i].Qmin     = gen_arr[git,5] / pm_data["baseMVA"]
    generators[i].Vg       = gen_arr[git,6]
    generators[i].mBase    = gen_arr[git,7]
    generators[i].status   = gen_arr[git,8]
    @assert generators[i].status==1
    generators[i].Pmax     = gen_arr[git,9]  / pm_data["baseMVA"]
    generators[i].Pmin     = gen_arr[git,10] / pm_data["baseMVA"]
    generators[i].Pc1      = gen_arr[git,11]
    generators[i].Pc2      = gen_arr[git,12]
    generators[i].Qc1min   = gen_arr[git,13]
    generators[i].Qc1max   = gen_arr[git,14]
    generators[i].Qc2min   = gen_arr[git,15]
    generators[i].Qc2max   = gen_arr[git,16]
    generators[i].gentype  = costgen_arr[git,1]
    generators[i].startup  = costgen_arr[git,2]
    generators[i].shutdown = costgen_arr[git,3]
    generators[i].n        = costgen_arr[git,4]
    if generators[i].gentype == 1
      generators[i].coeff = costgen_arr[git,5:end]
      error("Piecewise linear costs remains to be implemented.")
    else
      if generators[i].gentype == 2
        generators[i].coeff = costgen_arr[git,5:end]
        #println(generators[i].coeff, " ", length(generators[i].coeff), " ", generators[i].coeff[2])
      else
        error("Invalid generator cost model in the data.")
      end
    end
  end
=#

  #println(generators)
  #println(bus_ref)
  pm_buses = pm_data["bus"]
  pm_lines = pm_data["branch"]
  pm_gens = pm_data["gen"]
  nbuses, nlines, ngens = length(pm_buses), length(pm_lines), length(pm_gens)
 # compute reference bus
  bus_ref=-1
  for (k,bus) in pm_buses
    if bus["bus_type"]==3
      if bus_ref>0
        error("More than one reference bus present in the data")
      else
         bus_ref=bus["index"]
      end
    end
  end
  # set bus demands
      PD = zeros(nbuses); QD = zeros(nbuses)
      for (k,load_data) in pm_data["load"]
        PD[load_data["load_bus"]] = load_data["pd"]; QD[load_data["load_bus"]] = load_data["qd"]
      end
  # set squared bounds for bus voltage magnitude
      Wmin = zeros(nbuses); Wmax = zeros(nbuses)
      for (k,bus) in pm_buses
        Wmin[bus["index"]] = bus["vmin"]^2; Wmax[bus["index"]] = bus["vmax"]^2
      end
  # set bounds for generator power generation
      Pmin = zeros(ngens); Pmax = zeros(ngens)
      Qmin = zeros(ngens); Qmax = zeros(ngens)
      for (k,gen) in pm_gens
        Pmin[gen["index"]] = gen["pmin"]; Pmax[gen["index"]] = gen["pmax"]
        Qmin[gen["index"]] = gen["qmin"]; Qmax[gen["index"]] = gen["qmax"]
      end
  # set up the fromLines and toLines for each bus
      fromLines, toLines = mapLinesToBuses(pm_data)
      fromBus,toBus=zeros(Int,nlines),zeros(Int,nlines)
      for (k,line) in pm_lines 
        fromBus[line["index"]] = line["f_bus"]; toBus[line["index"]] = line["t_bus"]
      end
    # generators at each bus
      BusGeners = mapGenersToBuses(pm_data)
  # obtain entries of the admittance matrix
  Y_AC = computeAdmittances(pm_data, 0)
  Y_DC = computeAdmittances(pm_data, 2)
  return OPFData(nbuses, PD, QD, Wmin, Wmax, 
	nlines, fromLines, toLines, fromBus, toBus, 
	ngens, Pmin, Pmax, Qmin, Qmax, BusGeners,
	bus_ref, Y_AC, Y_DC)
end

function  computeAdmittances(pm_data, relaxType=0)
  Y=Dict()
  nlines = length(pm_data["branch"])
  Y["ffR"]=zeros(nlines)
  Y["ffI"]=zeros(nlines)
  Y["ttR"]=zeros(nlines)
  Y["ttI"]=zeros(nlines)
  Y["ftR"]=zeros(nlines)
  Y["ftI"]=zeros(nlines)
  Y["tfR"]=zeros(nlines)
  Y["tfI"]=zeros(nlines)

  for (k,line) in pm_data["branch"]
    @assert line["br_status"] == 1
    Ys = 1/(line["br_r"] + line["br_x"]*im)
    if relaxType==2  #DC relaxation
      Ys = 1/(0.0 + line["br_x"]*im)  #using DC, line resistances are set to zero
    end

    #assign nonzero tap ratio
    tap=line["tap"]

    #add phase shifters
    tap *= exp(line["shift"] * im)

    Yff = Ys + line["b_fr"]*im
    Ytt = Ys + line["b_to"]*im
    if relaxType==2
      Ytt = Ys #for DC, line susceptances are set to zero
    end
    Yff /= (tap*conj(tap))
    Yft = -Ys / conj(tap)
    Ytf = -Ys / tap

    #split into real and imag parts
    Y["ffR"][line["index"]] = real(Yff); Y["ffI"][line["index"]] = imag(Yff)
    Y["ttR"][line["index"]] = real(Ytt); Y["ttI"][line["index"]] = imag(Ytt)
    Y["tfR"][line["index"]] = real(Ytf); Y["tfI"][line["index"]] = imag(Ytf)
    Y["ftR"][line["index"]] = real(Yft); Y["ftI"][line["index"]] = imag(Yft)
  end

  nbuses = length(pm_data["bus"])
  Y["shR"] = zeros(nbuses)
  Y["shI"] = zeros(nbuses)
  for (k,shunt) in pm_data["shunt"]
    Y["shR"][shunt["shunt_bus"]] = shunt["gs"] 
    Y["shI"][shunt["shunt_bus"]] = shunt["bs"]
  end

  return Y
end

# Builds a map from lines to buses.
# For each line we store an array with zero or one element containing
# the  'From' and 'To'  bus number.
function mapLinesToBuses(pm_data)
  nbus = length(pm_data["bus"])
  FromLines = [Int[] for b in 1:nbus]
  ToLines   = [Int[] for b in 1:nbus]
  for (k,line) in pm_data["branch"]
      if(line["br_status"]>0)
          push!(FromLines[line["f_bus"]], line["index"])
          push!(ToLines[line["t_bus"]], line["index"])
      end
  end
  return FromLines,ToLines
end

# Builds a map between buses and generators.
# For each bus we keep an array of corresponding generators number (as array).
#
# (Can be more than one generator per bus, in future implementations at least)
function mapGenersToBuses(pm_data)
  gen2bus = [Int[] for b in 1:length(pm_data["bus"])]
  for (k,gen) in pm_data["gen"]
    push!(gen2bus[gen["gen_bus"]], gen["index"])
  end
  return gen2bus
end
