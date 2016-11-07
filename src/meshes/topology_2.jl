import Combinatorics: combinations

export MeshTopologyGeneric
export init!, connectivity!, connectivity
export build!, build, transpose!, transpose, intersection!, intersection

include("connectivity.jl")

#--------------------------#
# Type MeshTopologyGeneric #
#--------------------------#
# """
#   Stores the topology of a mesh as a set of incidence relations.
# """
type MeshTopologyGeneric{T<:AbstractFloat} <: AbstractMeshTopology
    dimension :: Int
    nodes :: Array{T, 2}
    connectivities :: Dict{Tuple{Int, Int}, MeshConnectivity}
end

# Outer Constructors
# -------------------
function MeshTopologyGeneric{T<:AbstractFloat}(dim::Integer, nodes::AbstractArray{T,2})
    connect = Dict{Tuple{Int,Int}, MeshConnectivity}()
    return MeshTopologyGeneric(dim, nodes, connect)
end

function MeshTopologyGeneric{T<:AbstractFloat,S<:Integer}(dim::Integer, nodes::AbstractArray{T,2},
                                                          cells::AbstractArray{S,2})
    mt = MeshTopologyGeneric(dim, nodes)
    init!(mt, cells)
    return mt
end

# Associated Methods
# -------------------
function Base.show(io::IO, mt::MeshTopologyGeneric)
    println(io, "MeshTopologyGeneric")
    println(io, "\tdimension : ", mt.dimension)
    println(io, "\tconnectivities : ",
            ["\n\t\t$connect" for connect in mt.connectivities]...)
end

function getConnectivity(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    haskey(mt.connectivities, (d,dd)) || @time connectivity!(mt, d, dd)
    return mt.connectivities[(d,dd)]
end

function getEntities(mt::MeshTopologyGeneric, d::Integer)
    if d == 0
        return collect(getConnectivity(mt, d, 0))
    else
        return hcat(getConnectivity(mt, d, 0)...)'
    end
end

"""

    getBoundary(mt::MeshTopologyGeneric, d::Integer)

  Determine the boundary entities of topological dimension `d`.
"""
function getBoundary(mt::MeshTopologyGeneric, d::Integer)
    D = getDim(mt) # mt.dimension
    D_Dm1 = getConnectivity(mt, D-1, D)
    bmask = [length(r) == 1 for r in D_Dm1]

    if d == D-1
        return bmask
    else
        d_Dm1 = getConnectivity(mt, d, D-1)
        bind = find(bmask)
        bmask = [any(i in bind for i in r) for r in d_Dm1]
    end
end

getBoundary(mt::MeshTopologyGeneric) = getBoundary(mt, getDim(mt)-1)

function getBoundary(mt::MeshTopologyGeneric, f::Function)
    bmask = getBoundary(mt)
    error("Not implemented yet...!")
end

"""

    init!{T<:Integer}(mt::MeshTopologyGeneric, cells::AbstractArray{T,2})

  Compute the incidence relations `D -> 0`, `0 -> D` and `D -> D` 
  from the given cell connectivity array where `D` is the dimension
  of the topology.
"""
function init!{T<:Integer}(mt::MeshTopologyGeneric, cells::AbstractArray{T,2})
    D = mt.dimension

    # compute 'D -> 0'
    ncells, nverts = size(cells)
    ind_D_0 = cells'[:]
    off_D_0 = collect(1:nverts:nverts*(ncells+1))
    
    mt.connectivities[(D,0)] = MeshConnectivity(D, 0, ind_D_0, off_D_0)

    return nothing
end

"""

    connectivity(mt::MeshTopologyGeneric, d::Integer, dd::Integer)

  Compute the incidence relation `d -> dd` by successive application
  of `build`, `transpose` and `intersection`.
"""
function connectivity!(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    d  != 0 && (haskey(mt.connectivities, (d, 0)) || build!(mt, d))
    dd != 0 && (haskey(mt.connectivities, (dd,0)) || build!(mt, dd))
    
    if !haskey(mt.connectivities, (d,dd))
        if d < dd
            connectivity!(mt, dd, d)
            transpose!(mt, dd, d)
        else
            D = mt.dimension
            ddd = (d == 0 && dd == 0) ? D : 0
            
            connectivity!(mt, d, ddd)
            connectivity!(mt, ddd, dd)
            intersection!(mt, d, dd, ddd)
        end
    end

    return mt.connectivities[(d,dd)]
end
function connectivity(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    prekeys = collect(keys(mt.connectivities))
    d_dd = connectivity!(mt, d, dd)
    postkeys = keys(mt.connectivities)
    for key in setdiff(postkeys, prekeys)
        delete!(mt.connectivities, key)
    end
    return d_dd
end

"""

    build!(mt::MeshTopologyGeneric, d::Integer)

  Compute the incidence relation `d -> 0` from `D -> 0`
  for `0 < d < D` where `D` is the dimension of the topology.
"""
function build!(mt::MeshTopologyGeneric, d::Integer)
    D = mt.dimension
    @assert 0 < d < D "0 < (d = $d) < D = $D"
    
    # compute the set of vertex sets for each cell
    print("vertex_sets: ")
    @time V = vertex_sets(mt, d)

    # take unique sorted index sets
    entities_d = unique(sort!.(vcat(V...)))
    nentities_d = length(entities_d)
    nverts_d = length(entities_d[1])
    
    # compute indices and offsets vectors
    ind_d_0 = vcat(entities_d...)
    off_d_0 = collect(1:nverts_d:nverts_d*(nentities_d+1))

    mt.connectivities[(d,0)] = MeshConnectivity(d, 0, ind_d_0, off_d_0)
end
function build(mt::MeshTopologyGeneric, d::Integer)
    prekeys = collect(keys(mt.connectivities))
    d_0 = build!(mt, d)
    postkeys = keys(mt.connectivities)
    for key in setdiff(postkeys, prekeys)
        delete!(mt.connectivities, key)
    end
    return d_0
end

"""

    transpose!(mt::MeshTopologyGeneric, d::Integer, dd::Integer)

  Compute the incidence relation `dd -> d` from `d -> dd`
  for `d > dd`.
"""
function transpose!(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    @assert d > dd
    # get connectivity d -> dd
    connect_d_dd = connectivity!(mt, d, dd)

    I, J = findn(connect_d_dd)

    # compute indices vector
    P = sortperm(J)
    ind_dd_d = similar(I)
    for i = 1:length(P)
        ind_dd_d[i] = I[P[i]]
    end

    # compute offsets vector
    off_dd_d = zeros(Int, maximum(J)+1)
    off_dd_d[1] = 1
    for j in J
        off_dd_d[j+1] += 1
    end
    cumsum!(off_dd_d, off_dd_d)

    mt.connectivities[(dd, d)] = MeshConnectivity(dd, d, ind_dd_d, off_dd_d)
end

function transpose(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    prekeys = collect(keys(mt.connectivities))
    dd_d = transpose!(mt, d, dd)
    postkeys = keys(mt.connectivities)
    for key in setdiff(postkeys, prekeys)
        delete!(mt.connectivities, key)
    end
    return dd_d
end

"""

    intersection!(mt::MeshTopologyGeneric, d::Integer, dd::Integer, ddd::Integer)

  Compute the incidence relation `d -> dd` from
  `d -> ddd` and `ddd -> dd` for `d >= dd`.
"""
function intersection(mt::MeshTopologyGeneric, d::Integer, dd::Integer, ddd::Integer)
    prekeys = collect(keys(mt.connectivities))
    d_dd = intersection!(mt, d, dd, ddd)
    postkeys = keys(mt.connectivities)
    for key in setdiff(postkeys, prekeys)
        delete!(mt.connectivities, key)
    end
    return d_dd
end

function intersection!(mt::MeshTopologyGeneric, d::Integer, dd::Integer, ddd::Integer)
    @assert d >= d
    d_ddd  = connectivity!(mt, d,   ddd)
    ddd_dd = connectivity!(mt, ddd, dd)

    # initialize indices vector for new connectivity 'd -> dd'
    ind_d_dd = [Array{Int,1}() for i = 1:size(d_ddd, 1)]
    
    if d == dd
        fill_intersection!(ind_d_dd, d_ddd, ddd_dd)
    elseif d > dd
        d_0  = connectivity!(mt, d,   0)
        dd_0 = connectivity!(mt, dd,  0)

        fill_intersection!(ind_d_dd, d_ddd, ddd_dd, d_0, dd_0)
    end

    # compute offsets
    off = vcat(1, 1 + cumsum([length(i) for i in ind_d_dd]))

    mt.connectivities[(d,dd)] = MeshConnectivity(d, dd, vcat(ind_d_dd...), off)
end

@noinline function fill_intersection!(d_dd::Array{Array{Int,1},1},
                                      d_ddd::MeshConnectivity, ddd_dd::MeshConnectivity)
    for (i, d_ddd_i) in enumerate(d_ddd)
        for k in d_ddd_i
            for j in ddd_dd[k]
                if (i != j) && !(j in d_dd[i])
                    #push!(d_dd[i], j)
                    resize!(d_dd[i], length(d_dd[i])+1)
                    d_dd[i][end] = j
                end
            end
        end
    end
    return nothing
end

@noinline function fill_intersection!(d_dd::Array{Array{Int,1},1},
                            d_ddd::MeshConnectivity, ddd_dd::MeshConnectivity,
                            d_0::MeshConnectivity, dd_0::MeshConnectivity)
    for (i, d_ddd_i) in enumerate(d_ddd)
        for k in d_ddd_i
            for j in ddd_dd[k]
                if !(j in d_dd[i]) && issubset(dd_0[j], d_0[i])
                    push!(d_dd[i], j)
                end
            end
        end
    end
    return nothing
end

"""

    vertex_sets(mt::MeshTopologyGeneric, d::Integer, sorted::Bool=true)

  Compute for each cell the set of vertex sets 
  incident to the entities of topological dimension `d`.
  If `sorted` is `true` (default) the indices will be sorted 
  in increasing order.
"""
function vertex_sets(mt::MeshTopologyGeneric, d::Integer, sorted::Bool=true)
    D = mt.dimension
    D_0 = getConnectivity(mt, D, 0)

    ncells = size(D_0, 1)
    nverts = d+1
    #combs = hcat(combinations(1:D+1, nverts)...)'
    combs = collect(combinations(1:D+1, nverts))
    ncombs = size(combs, 1)

    V = [[D_0[i,comb] for comb in combs] for i = 1:ncells]
    return V
end
export vertex_sets

