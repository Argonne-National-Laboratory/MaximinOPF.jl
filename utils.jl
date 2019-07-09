using LightGraphs

function get_graph(mat::SparseMatrixCSC{Float64,Int64})
    rows = rowvals(mat)
    vals = nonzeros(mat)
    m, n = size(mat)
    g = SimpleGraph(m)
    for i = 1:n
        for j in nzrange(mat, i)
            # if abs(vals[j]) > 0
                add_edge!(g, rows[j], i)
            # end
        end
    end
    return g
end

# Get the n-by-n chordal extension
function get_chordal_extension(opfdata)
    N, L, fromBus, toBus = opfdata.N, opfdata.L, opfdata.fromBus, opfdata.toBus  
    # Laplacian graph
    I = Int64[]; J = Int64[]; V = Float64[]
    for l in L
        fid = fromBus[l]
        tid = toBus[l]
        push!(I,fid); push!(J,tid); push!(V,-1)
        push!(I,tid); push!(J,fid); push!(V,-1)
    end
    A = sparse(I,J,V)
    for i in N
        push!(I,i); push!(J,i)
        push!(V,-sum(A[i,:])+1)
    end
    A = sparse(I,J,V)
    C = sparse(cholfact(A))
    return C
end

# Get the 2n-by-2n chordal extension
function get_chordal_extension_complex(opfdata)
    N, L, fromBus, toBus = opfdata.N, opfdata.L, opfdata.fromBus, opfdata.toBus  
    # Laplacian graph
    num = length(N)
    I = Int64[]; J = Int64[]; V = Float64[]
    for l in L
        fid = fromBus[l]
        tid = toBus[l]
        push!(I,fid); push!(J,tid); push!(V,-1)
        push!(I,tid); push!(J,fid); push!(V,-1)
        push!(I,fid); push!(J,num+tid); push!(V,-1)
        push!(I,num+tid); push!(J,fid); push!(V,-1)
        push!(I,tid); push!(J,num+fid); push!(V,-1)
        push!(I,num+fid); push!(J,tid); push!(V,-1)
        push!(I,num+fid); push!(J,num+tid); push!(V,-1)
        push!(I,num+tid); push!(J,num+fid); push!(V,-1)
    end
    A = sparse(I,J,V)
    for i in 1:(2*num)
        push!(I,i); push!(J,i)
        push!(V,-sum(A[i,:])+1)
    end
    A = sparse(I,J,V,2*num,2*num)
    # λ, ϕ = eigs(A)
    # @show λ
    C = sparse(cholfact(A))
    return C
end
