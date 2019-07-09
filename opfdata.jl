type Bus
  bus_i::Int
  bustype::Int
  Pd::Float64
  Qd::Float64
  Gs::Float64
  Bs::Float64
  area::Int
  Vm::Float64
  Va::Float64
  baseKV::Float64
  zone::Int
  Vmax::Float64
  Vmin::Float64
end

type Line
  from::Int
  to::Int
  r::Float64
  x::Float64
  b::Float64
  rateA::Float64
  rateB::Float64
  rateC::Float64
  ratio::Float64 #TAP
  angle::Float64 #SHIFT
  status::Int
  angmin::Float64
  angmax::Float64
end
Line() = Line(0,0,0.,0.,0.,0.,0.,0.,0.,0.,0,0.,0.)

type Gener
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

type Admittances
    YffR::Array{Float64}
    YffI::Array{Float64}
    YttR::Array{Float64}
    YttI::Array{Float64}
    YftR::Array{Float64}
    YftI::Array{Float64}
    YtfR::Array{Float64}
    YtfI::Array{Float64}
    YshR::Array{Float64}
    YshI::Array{Float64}
end

type OPFData
  buses::Array{Bus}
  nbuses::Int
  N
  PD::Array
  QD::Array
  Wmin::Array
  Wmax::Array
  lines::Array{Line}
  nlines::Int
  L
  fromLines::Array         #From lines for each bus (Array of Array)
  toLines::Array           #To lines for each bus (Array of Array)
  fromBus::Array
  toBus::Array
  generators::Array{Gener}
  ngens::Int
  G
  Pmin::Array{Float64}
  Pmax::Array{Float64}
  Qmin::Array{Float64}
  Qmax::Array{Float64}
  BusGeners
  bus_ref::Int
  Y
  # baseMVA::Float64
#  BusIdx::Dict{Int,Int}    #map from bus ID to bus index
#  BusGenerators::Array     #list of generators for each bus (Array of Array)
  #admittancesAC::Admittances   # admittances
  #admittancesDC::Admittances   # admittances
  
end


function opf_loaddata(case_name, lineOff=Line())

  data_path = "data/case"
  #
  # load buses
  #
  # bus_arr = readdlm("data/" * case_name * ".bus")
  bus_arr = readdlm(data_path * string(case_name) * ".bus")
  num_buses = size(bus_arr,1)
  buses = Array{Bus}(num_buses)
  bus_ref=-1
  for i in 1:num_buses
    @assert bus_arr[i,1]>0  #don't support nonpositive bus ids
    buses[i] = Bus(bus_arr[i,1:13]...)
    buses[i].Va *= pi/180
    if buses[i].bustype==3
      if bus_ref>0
        error("More than one reference bus present in the data")
      else
         bus_ref=i
      end
    end
    #println("bus ", i, " ", buses[i].Vmin, "      ", buses[i].Vmax)
  end

  #
  # load branches/lines
  #
  # branch_arr = readdlm("data/" * case_name * ".branch")
  branch_arr = readdlm(data_path * string(case_name) * ".branch")
  num_lines = size(branch_arr,1)
  lines_on = find((branch_arr[:,11].>0) .& ((branch_arr[:,1].!=lineOff.from) .| (branch_arr[:,2].!=lineOff.to)) )
  num_on   = length(lines_on)

  if lineOff.from>0 && lineOff.to>0
    println("opf_loaddata: was asked to remove line from,to=", lineOff.from, ",", lineOff.to)
    #println(lines_on, branch_arr[:,1].!=lineOff.from, branch_arr[:,2].!=lineOff.to)
  end
  if length(find(branch_arr[:,11].==0))>0
    println("opf_loaddata: ", num_lines-length(find(branch_arr[:,11].>0)), " lines are off and will be discarded (out of ", num_lines, ")")
  end

  lines = Array{Line}(num_on)

  lit=0
  for i in lines_on
    @assert branch_arr[i,11] == 1  #should be on since we discarded all other
    lit += 1
    lines[lit] = Line(branch_arr[i, 1:13]...)
    if lines[lit].angmin>-360 || lines[lit].angmax<360
      error("Bounds of voltage angles are still to be implemented.")
    end

  end
  @assert lit == num_on

  #
  # load generators
  #
  # gen_arr = readdlm("data/" * case_name * ".gen")
  gen_arr = readdlm(data_path * string(case_name) * ".gen")
  # costgen_arr = readdlm("data/" * case_name * ".gencost")
  costgen_arr = readdlm(data_path * string(case_name) * ".gencost")
  num_gens = size(gen_arr,1)

  baseMVA=100

  @assert num_gens == size(costgen_arr,1)

  gens_on=find(gen_arr[:,8]); num_on=length(gens_on)
  if num_gens-num_on>0
    println("loaddata: ", num_gens-num_on, " generators are off and will be discarded (out of ", num_gens, ")")
  end

  generators = Array{Gener}(num_on)
  i=0
  for git in gens_on
    i += 1
    generators[i] = Gener(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, Array{Int}(0)) #gen_arr[i,1:end]...)

    generators[i].bus      = gen_arr[git,1]
    generators[i].Pg       = gen_arr[git,2] / baseMVA
    generators[i].Qg       = gen_arr[git,3] / baseMVA
    generators[i].Qmax     = gen_arr[git,4] / baseMVA
    generators[i].Qmin     = gen_arr[git,5] / baseMVA
    generators[i].Vg       = gen_arr[git,6]
    generators[i].mBase    = gen_arr[git,7]
    generators[i].status   = gen_arr[git,8]
    @assert generators[i].status==1
    generators[i].Pmax     = gen_arr[git,9]  / baseMVA
    generators[i].Pmin     = gen_arr[git,10] / baseMVA
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

  admittancesAC = computeAdmittances(lines, buses, baseMVA, 0)
  admittancesDC = computeAdmittances(lines, buses, baseMVA, 2)
  #println(generators)
  #println(bus_ref)
  nbuses, nlines, ngens = length(buses), length(lines), length(generators)
  N = 1:nbuses; L = 1:nlines; G = 1:ngens
  # set bus demands
      PD = zeros(nbuses); QD = zeros(nbuses)
      for i in N
        PD[i] = buses[i].Pd / baseMVA; QD[i] = buses[i].Qd / baseMVA
      end
  # set squared bounds for bus voltage magnitude
      Wmin = zeros(nbuses); Wmax = zeros(nbuses)
      for i in N
        Wmin[i] = (buses[i].Vmin)^2; Wmax[i] = (buses[i].Vmax)^2
      end
  # set bounds for generator power generation
      Pmin = zeros(ngens); Pmax = zeros(ngens)
      Qmin = zeros(ngens); Qmax = zeros(ngens)
      for g in G
        Pmin[g] = generators[g].Pmin; Pmax[g] = generators[g].Pmax
        Qmin[g] = generators[g].Qmin; Qmax[g] = generators[g].Qmax
      end
    println("Done with initial setup.")
  # build a dictionary between buses ids and their indexes
    busIdx = mapBusIdToIdx(buses)
  # set up the fromLines and toLines for each bus
      fromLines, toLines = mapLinesToBuses(buses, lines, busIdx)
      fromBus=zeros(Int,nlines); toBus=zeros(Int,nlines)
      for l in L
        fromBus[l] = busIdx[lines[l].from]; toBus[l] = busIdx[lines[l].to]
      end
    # generators at each bus
      BusGeners = mapGenersToBuses(buses, generators, busIdx)
    # obtain entries of the admittance matrix
      Y = Dict()  # Admittances
      Y["ffR"] = admittancesAC.YffR; Y["ffI"] = admittancesAC.YffI;
      Y["ttR"] = admittancesAC.YttR; Y["ttI"] = admittancesAC.YttI;
      Y["ftR"] = admittancesAC.YftR; Y["ftI"] = admittancesAC.YftI;
      Y["tfR"] = admittancesAC.YtfR; Y["tfI"] = admittancesAC.YtfI;
      Y["shR"] = admittancesAC.YshR; Y["shI"] = admittancesAC.YshI;
  return OPFData(buses, nbuses, N, PD, QD, Wmin, Wmax, 
	lines, nlines, L, fromLines, toLines, fromBus, toBus, 
	generators, ngens, G, Pmin, Pmax, Qmin, Qmax, BusGeners,
	bus_ref, Y)
end

function  computeAdmittances(lines, buses, baseMVA, relaxType=0)
  nlines = length(lines)
  YffR=Array{Float64}(nlines)
  YffI=Array{Float64}(nlines)
  YttR=Array{Float64}(nlines)
  YttI=Array{Float64}(nlines)
  YftR=Array{Float64}(nlines)
  YftI=Array{Float64}(nlines)
  YtfR=Array{Float64}(nlines)
  YtfI=Array{Float64}(nlines)

  for i in 1:nlines
    @assert lines[i].status == 1
    Ys = 1/(lines[i].r + lines[i].x*im)
    if relaxType==2  #DC relaxation
      Ys = 1/(0.0 + lines[i].x*im)  #using DC, line resistances are set to zero
    end
    #assign nonzero tap ratio
    tap = lines[i].ratio==0?1.0:lines[i].ratio

    #add phase shifters
    tap *= exp(lines[i].angle * pi/180 * im)

    Ytt = Ys + lines[i].b/2*im
    if relaxType==2
      Ytt = Ys #for DC, line susceptances are set to zero
    end
    Yff = Ytt / (tap*conj(tap))
    Yft = -Ys / conj(tap)
    Ytf = -Ys / tap

    #split into real and imag parts
    YffR[i] = real(Yff); YffI[i] = imag(Yff)
    YttR[i] = real(Ytt); YttI[i] = imag(Ytt)
    YtfR[i] = real(Ytf); YtfI[i] = imag(Ytf)
    YftR[i] = real(Yft); YftI[i] = imag(Yft)
    #@printf("[%4d]  tap=%12.9f   %12.9f    lines[i].b=%f i\n", i, real(tap), imag(tap), lines[i].b/2);
  end

  nbuses = length(buses)
  YshR = Array{Float64}(nbuses)
  YshI = Array{Float64}(nbuses)
  for i in 1:nbuses
    YshR[i] = buses[i].Gs / baseMVA
    YshI[i] = buses[i].Bs / baseMVA
    #@printf("[%4d]   Ysh  %15.12f + %15.12f i \n", i, YshR[i], YshI[i])
  end

  @assert 0==length(find(isnan.(YffR)))+length(find(isinf.(YffR)))
  @assert 0==length(find(isnan.(YffI)))+length(find(isinf.(YffI)))
  @assert 0==length(find(isnan.(YttR)))+length(find(isinf.(YttR)))
  @assert 0==length(find(isnan.(YttI)))+length(find(isinf.(YttI)))
  @assert 0==length(find(isnan.(YftR)))+length(find(isinf.(YftR)))
  @assert 0==length(find(isnan.(YftI)))+length(find(isinf.(YftI)))
  @assert 0==length(find(isnan.(YtfR)))+length(find(isinf.(YtfR)))
  @assert 0==length(find(isnan.(YtfI)))+length(find(isinf.(YtfI)))
  @assert 0==length(find(isnan.(YshR)))+length(find(isinf.(YshR)))
  @assert 0==length(find(isnan.(YshI)))+length(find(isinf.(YshI)))

  return Admittances(YffR, YffI, YttR, YttI, YftR, YftI, YtfR, YtfI, YshR, YshI)
end

# Builds a map from lines to buses.
# For each line we store an array with zero or one element containing
# the  'From' and 'To'  bus number.
function mapLinesToBuses(buses, lines, busDict)
  nbus = length(buses)
  FromLines = [Int[] for b in 1:nbus]
  ToLines   = [Int[] for b in 1:nbus]
  for i in 1:length(lines)
      if(lines[i].status>0)
          busID = busDict[lines[i].from]
          @assert 1<= busID <= nbus
          push!(FromLines[busID], i)

          busID = busDict[lines[i].to]
          @assert 1<= busID  <= nbus
          push!(ToLines[busID], i)
      end
  end
  return FromLines,ToLines
end

# Builds a mapping between bus ids and bus indexes
#
# Returns a dictionary with bus ids as keys and bus indexes as values
function mapBusIdToIdx(buses)
  dict = Dict{Int,Int}()
  for b in 1:length(buses)
    @assert !haskey(dict,buses[b].bus_i)
    dict[buses[b].bus_i] = b
  end
  return dict
end


# Builds a map between buses and generators.
# For each bus we keep an array of corresponding generators number (as array).
#
# (Can be more than one generator per bus)
function mapGenersToBuses(buses, generators,busDict)
  gen2bus = [Int[] for b in 1:length(buses)]
  for g in 1:length(generators)
    busID = busDict[ generators[g].bus ]
    #@assert(0==length(gen2bus[busID])) #at most one generator per bus
    push!(gen2bus[busID], g)
  end
  return gen2bus
end
