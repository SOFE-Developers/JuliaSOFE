import Combinatorics: combinations

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

    # println("Time transpose! in init!")
    # @time transpose!(mt, D, 0)       # compute '0 -> D' by transposing 'D -> 0'
    # println("Time intersection! in init!")
    # @time intersection!(mt, D, D, 0) # and 'D -> D' by intersecting 'D -> 0' and '0 -> D'

    return nothing
end

"""

    connectivity!(mt::MeshTopologyGeneric, d::Integer, dd::Integer)

  Compute the incidence relation `d -> dd` by successive application
  of `build`, `transpose` and `intersection`.
"""
function connectivity!(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    d  != 0 && (haskey(mt.connectivities, (d, 0)) || build!(mt, d))
    dd != 0 && (haskey(mt.connectivities, (dd,0)) || build!(mt, dd))
    
    if haskey(mt.connectivities, (d,dd))
        return
    end

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

"""

    build!(mt::MeshTopologyGeneric, d::Integer)

  Compute the incidence relation `D -> d` and `d -> 0`
  from `D -> 0` and `D -> D` for `0 < d < D` where `D`
  is the dimension of the topology.
"""
function build!(mt::MeshTopologyGeneric, d::Integer)
    D = mt.dimension
    @assert 0 < d < D "0 < (d = $d) < D = $D"
    
    # compute the set of vertex sets for each cell
    V = vertex_sets(mt, d)

    D_D = getConnectivity(mt, D, D)
    
    # initialize index vectors for new connectivities
    # 'D -> d' and 'd -> 0'
    nentities_D = length(D_D.offsets) - 1
    nentities_d = size(unique(hcat(vcat(V...)...), 2), 2)
    ind_D_d = [Array{Int,1}() for i = 1:nentities_D]
    ind_d_0 = [Array{Int,1}() for i = 1:nentities_d]
    
    k = 1
    for (i, D_D_i) in enumerate(D_D)
        if i == 1
            for (j, vj) in enumerate(V[i])
                push!(ind_D_d[1], j)
                ind_d_0[j] = vj
                k += 1
            end
        else
            for j in D_D_i
                if j < i
                    for vi in V[i]
                        if vi in V[j]
                            l = indexin([vi], ind_d_0)[1]
                            push!(ind_D_d[i], l)
                        else
                            if !(vi in ind_d_0)
                                push!(ind_D_d[i], k)
                                ind_d_0[k] = vi
                                k += 1
                            end
                        end
                    end
                end
            end
        end
    end
    
    # compute offsets
    off_D_d = vcat(1, 1 + Array{Int,1}(cumsum([length(i) for i in ind_D_d])))
    off_d_0 = vcat(1, 1 + Array{Int,1}(cumsum([length(i) for i in ind_d_0])))
    
    mt.connectivities[(D,d)] = MeshConnectivity(D, d, vcat(ind_D_d...), off_D_d)
    mt.connectivities[(d,0)] = MeshConnectivity(d, 0, vcat(ind_d_0...), off_d_0)
end

"""

    transpose!(mt::MeshTopologyGeneric, d::Integer, dd::Integer)

  Compute the incidence relation `dd -> d` from `d -> dd`
  for `d > dd`.
"""
# function transpose!(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
#     @assert d > dd
#     # get connectivity d -> dd
#     connect_d_dd = getConnectivity(mt, d, dd)

#     # initialize indices vector for new connectivity 'dd -> d'
#     nentities_dd = maximum(connect_d_dd.indices)
#     ind_dd_d = [Array{Int,1}() for i = 1:nentities_dd]

#     # fill indices vector
#     for (j, d_dd_j) in enumerate(connect_d_dd)
#         for i in d_dd_j
#             push!(ind_dd_d[i], j)
#         end
#     end

#     # compute offsets
#     off = vcat(1, 1 + Array{Int,1}(cumsum([length(i) for i in ind_dd_d])))

#     mt.connectivities[(dd,d)] = MeshConnectivity(dd, d, vcat(ind_dd_d...), off)
# end

function ntranspose(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    @assert d > dd
    # get connectivity d -> dd
    connect_d_dd = getConnectivity(mt, d, dd)

    I, J = findn(connect_d_dd)

    # compute indices vector
    P = sortperm(J)
    ind_dd_d = similar(I)
    for i = 1:length(P)
        ind_dd_d[i] = I[P[i]]
    end

    # compute offsets vector
    sort!(J)
    off_dd_d = zeros(Int, maximum(J)+1)
    off_dd_d[1] = 1
    for j in J
        off_dd_d[j+1] += 1
    end
    cumsum!(off_dd_d, off_dd_d)

    return MeshConnectivity(dd, d, ind_dd_d, off_dd_d)
end

function transpose!(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    @assert d > dd
    # get connectivity d -> dd
    connect_d_dd = getConnectivity(mt, d, dd)

    # initialize offsets vectors for new connectivity 'dd -> d'
    nentities_dd = maximum(connect_d_dd.indices)
    off_dd_d = ones(Int, nentities_dd+1)

    # fill offsets vector
    # (and compute number of connections 'dd -> d')
    nconnect_dd_d = fill_offsets!(off_dd_d, connect_d_dd)
    
    # initialize indices vector for new connectivity 'dd -> d'
    len_dd_d = zeros(Int, nentities_dd)
    ind_dd_d = zeros(Int, nconnect_dd_d)
    
    # fill indices vector
    fill_transpose!(ind_dd_d, len_dd_d, off_dd_d, connect_d_dd)
    
    #mt.connectivities[(dd,d)] = MeshConnectivity(dd, d, ind_dd_d, off_dd_d)
    return MeshConnectivity(dd, d, ind_dd_d, off_dd_d)
end

@noinline function fill_offsets!(off::Array{Int,1}, connect::MeshConnectivity)
    nentities = length(off)
    nconnect = 0
    for ci in connect
        for j in ci
            nconnect += 1
            for k = j+1:nentities
                off[k] += 1
            end
        end
    end
    return nconnect
end

@noinline function fill_transpose!(ind::Array{Int,1}, len::Array{Int,1}, off::Array{Int,1}, connect::MeshConnectivity)
    for (i, ci) in enumerate(connect)
        for j in ci
            ind[off[j]+len[j]] = i
            len[j] += 1
        end
    end
    return nothing
end

"""

    intersection!(mt::MeshTopologyGeneric, d::Integer, dd::Integer, ddd::Integer)

  Compute the incidence relation `d -> dd` from
  `d -> ddd` and `ddd -> dd` for `d >= dd`.
"""
# function intersection!(mt::MeshTopologyGeneric, d::Integer, dd::Integer, ddd::Integer)
#     @assert d >= d
#     d_ddd  = getConnectivity(mt, d,   ddd)
#     ddd_dd = getConnectivity(mt, ddd, dd)
#     c_ddd_dd = collect(ddd_dd)
    
#     if d > dd
#         d_0  = getConnectivity(mt, d,   0)
#         dd_0 = getConnectivity(mt, dd,  0)

#         c_d_0 = collect(d_0)
#         c_dd_0 = collect(dd_0)
#     end
    
#     # initialize indices vector for new connectivity 'd -> dd'
#     nentities_d = length(d_ddd.offsets) - 1
#     ind_d_dd = [Array{Int,1}() for i = 1:nentities_d]
    
#     for (i, d_ddd_i) in enumerate(d_ddd)
#         for k in d_ddd_i
#             for j in c_ddd_dd[k]
#                 if (d == dd && i != j) || (d > dd && issubset(c_dd_0[j], c_d_0[i]))
#                     if !(j in ind_d_dd[i])
#                         push!(ind_d_dd[i], j)
#                     end
#                 end
#             end
#         end
#     end

#     # compute offsets
#     off = vcat(1, 1 + cumsum([length(i) for i in ind_d_dd]))

#     mt.connectivities[(d,dd)] = MeshConnectivity(d, dd, vcat(ind_d_dd...), off)
# end

function intersection!(mt::MeshTopologyGeneric, d::Integer, dd::Integer, ddd::Integer)
    @assert d >= d
    d_ddd  = getConnectivity(mt, d,   ddd)
    ddd_dd = getConnectivity(mt, ddd, dd)
    c_ddd_dd = collect(ddd_dd)

    # initialize indices vector for new connectivity 'd -> dd'
    nentities_d = length(d_ddd.offsets) - 1
    ind_d_dd = [Array{Int,1}() for i = 1:nentities_d]
    
    if d == dd
        @time fill_intersection!(ind_d_dd, d_ddd, c_ddd_dd)
    elseif d > dd
        c_d_0  = collect(getConnectivity(mt, d,   0))
        c_dd_0 = collect(getConnectivity(mt, dd,  0))

        fill_intersection!(ind_d_dd, d_ddd, c_ddd_dd, c_d_0, c_dd_0)
    end

    # compute offsets
    off = vcat(1, 1 + cumsum([length(i) for i in ind_d_dd]))

    mt.connectivities[(d,dd)] = MeshConnectivity(d, dd, vcat(ind_d_dd...), off)
end

@noinline function fill_intersection!(d_dd::Array{Array{Int,1},1},
                            d_ddd::MeshConnectivity, ddd_dd::Array{Array{Int,1},1})
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
                            d_ddd::MeshConnectivity, ddd_dd::Array{Array{Int,1},1},
                            d_0::Array{Array{Int,1},1}, dd_0::Array{Array{Int,1},1})
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

    ncells = length(D_0.offsets) - 1
    combs = combinations(1:Int(D+1), Int(d+1))

    V = [Array{Array{Int,1},1}() for i in 1:ncells]
    for (i, D_0_i) in enumerate(D_0)
        for comb in combs
            if sorted
                push!(V[i], sort([D_0_i[j] for j in comb]))
            else
                push!(V[i], [D_0_i[j] for j in comb])
            end
        end
    end

    return V
end

