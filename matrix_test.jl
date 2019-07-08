using GraphLayout

include("opfdata.jl")
include("utils.jl")

CASE_NUM = 9
opfdata = opf_loaddata(CASE_NUM)

# define data
lines, buses = opfdata.lines, opfdata.buses
nbuses, nlines = length(buses), length(lines)
N = 1:nbuses; L = 1:nlines

chordal = get_chordal_extension_complex(N, L, lines)
# @show chordal
@show full(chordal)

g = get_graph(chordal)
max_cliques = maximal_cliques(g)
@show max_cliques

A = Matrix(adjacency_matrix(g))
loc_x, loc_y = layout_spring_adj(A)
draw_layout_adj(A, loc_x, loc_y, filename="chordal.svg")

m,n=size(A)
for i=1:m
    for j in (i+1):n
        if A[i,j] > 0
            print(",($i,$j)")
        end
    end
end
