import Combinatorics: combinations

#-----------------------#
# Type MeshConnectivity #
#-----------------------#
"""
Stores the incidence relation 'd -> dd' for a fixed pair
of topological dimensions '(d,dd)'.
"""
type MeshConnectivity
    dims :: Tuple{Int, Int}
    indices :: Array{Int, 1}
    offsets :: Array{Int, 1}

    MeshConnectivity{T<:Integer}(dims::Tuple{Integer,Integer}, ind::Array{T,1}, off::Array{T,1}) =
        new(Tuple{Int,Int}(dims), Array{Int,1}(ind), Array{Int,1}(off))
end

# Outer Constructors
# -------------------
MeshConnectivity(dims::Tuple{Integer,Integer}) = MeshConnectivity(dims, Array{Integer,1}(), Array{Integer,1}())
MeshConnectivity(d1::Integer, d2::Integer) = MeshConnectivity((d1,d2))
MeshConnectivity{T<:Integer}(d1::Integer, d2::Integer, ind::Array{T,1}, off::Array{T,1}) =
    MeshConnectivity((d1,d2), ind, off)

# Associated Methods
# -------------------
@inline function Base.size(mc::MeshConnectivity)
    nentities_d  = max(0, length(mc.offsets) - 1)
    nentities_dd = length(mc.indices) > 0 ? maximum(mc.indices) : 0
    return (nentities_d, nentities_dd)
end
@inline Base.size(mc::MeshConnectivity, d::Integer) = size(mc)[d]

function Base.show(io::IO, mc::MeshConnectivity)
    println(io, "MeshConnectivity: ", mc.dims[1], " -> ", mc.dims[2])
end

# provide an iteration interface for `MeshConnectivity` instances.
Base.start(::MeshConnectivity) = 1

function Base.next(mc::MeshConnectivity, state::Integer)
    i, j = mc.offsets[state:state+1]
    c = view(mc.indices, i:j-1)
    # c = mc.indices[i:j-1]
    return (c, state+1)
end

Base.done(mc::MeshConnectivity, state::Integer) = state > length(mc)
Base.eltype(::Type{MeshConnectivity}) = Array{Int,1}
Base.length(mc::MeshConnectivity) = size(mc, 1)

#--------------------------#
# Type MeshTopologyGeneric #
#--------------------------#
"""
Stores the topology of a mesh as a set of incidence relations.
"""
type MeshTopologyGeneric <: AbstractMeshTopology
    dimension :: Int
    connectivities :: Dict{Tuple{Int, Int}, MeshConnectivity}

    MeshTopologyGeneric(dim::Integer, connect::Dict{Tuple{Int,Int}, MeshConnectivity}) =
        new(Int(dim), connect)
end

# Outer Constructors
# -------------------
MeshTopologyGeneric(dim::Integer) = MeshTopologyGeneric(dim, Dict{Tuple{Int,Int}, MeshConnectivity}())

# Associated Methods
# -------------------
function Base.show(io::IO, mt::MeshTopologyGeneric)
    println(io, "MeshTopologyGeneric")
    println(io, "\tdimension : ", mt.dimension)
    println(io, "\tconnectivities : ",
            ["\n\t\t$connect" for connect in mt.connectivities]...)
end

"""
Return the spatial dimension of the mesh topology,
i.e. he topological dimension of the cells.
"""
getDim(mt::MeshTopologyGeneric) = mt.dimension

"""

    getConnectivity(mt::MeshTopologyGeneric, d::Integer, dd::Integer)

Return the incidence relation `d -> dd` computing it first if necessary.
"""
function getConnectivity(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    haskey(mt.connectivities, (d,dd)) || connectivity!(mt, d, dd)
    return mt.connectivities[(d,dd)]
end

"""

    getEntities(mt::MeshTopologyGeneric, d::Integer)

Return the vertex index connectivity array for the mesh
entities of topological dimension `d`.
"""
function getEntities(mt::MeshTopologyGeneric, d::Integer)
    return collect(getConnectivity(mt, d, 0))
end

"""

    getNumber(mt::MeshTopologyGeneric, d::Integer)

Return the number of `d`-dimensional mesh entities.
"""
function getNumber(mt::MeshTopologyGeneric, d::Integer)
    # return size(getConnectivity(mt, d, 0), 1)
    return size(getEntities(mt, d), 1)
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

    transpose!(mt, D, 0)       # compute '0 -> D' by transposing 'D -> 0'
    intersection!(mt, D, D, 0) # and 'D -> D' by intersecting 'D -> 0' and '0 -> D'

    return
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
                            # l = find([v == vi for v in ind_d_0])[1]
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
function transpose!(mt::MeshTopologyGeneric, d::Integer, dd::Integer)
    @assert d > dd
    # get connectivity d -> dd
    connect_d_dd = getConnectivity(mt, d, dd)

    # initialize indices vector for new connectivity 'dd -> d'
    nentities_dd = maximum(connect_d_dd.indices)
    ind_dd_d = [Array{Int,1}() for i = 1:nentities_dd]

    # fill indices vector
    for (j, d_dd_j) in enumerate(connect_d_dd)
        for i in d_dd_j
            push!(ind_dd_d[i], j)
        end
    end

    # compute offsets
    off = vcat(1, 1 + Array{Int,1}(cumsum([length(i) for i in ind_dd_d])))

    mt.connectivities[(dd,d)] = MeshConnectivity(dd, d, vcat(ind_dd_d...), off)
end

"""

    intersection!(mt::MeshTopologyGeneric, d::Integer, dd::Integer, ddd::Integer)

Compute the incidence relation `d -> dd` from
`d -> ddd` and `ddd -> dd` for `d >= dd`.
"""
function intersection!(mt::MeshTopologyGeneric, d::Integer, dd::Integer, ddd::Integer)
    @assert d >= d
    d_ddd  = getConnectivity(mt, d,   ddd)
    ddd_dd = getConnectivity(mt, ddd, dd)
    c_ddd_dd = collect(ddd_dd)
    
    if d > dd
        d_0  = getConnectivity(mt, d,   0)
        dd_0 = getConnectivity(mt, dd,  0)

        c_d_0 = collect(d_0)
        c_dd_0 = collect(dd_0)
    end
    
    # initialize indices vector for new connectivity 'd -> dd'
    nentities_d = length(d_ddd.offsets) - 1
    ind_d_dd = [Array{Int,1}() for i = 1:nentities_d]
    
    for (i, d_ddd_i) in enumerate(d_ddd)
        for k in d_ddd_i
            for j in c_ddd_dd[k]
                if (d == dd && i != j) || (d > dd && issubset(c_dd_0[j], c_d_0[i]))
                    if !(j in ind_d_dd[i])
                        push!(ind_d_dd[i], j)
                    end
                end
            end
        end
    end

    # compute offsets
    off = vcat(1, 1 + cumsum([length(i) for i in ind_d_dd]))

    mt.connectivities[(d,dd)] = MeshConnectivity(d, dd, vcat(ind_d_dd...), off)
end

"""

    vertex_sets(mt::MeshTopologyGeneric, d::Integer)

Compute for each cell the set of vertex sets 
incident to the entities of topological dimension `d`.
"""
function vertex_sets(mt::MeshTopologyGeneric, d::Integer)
    D = mt.dimension
    D_0 = getConnectivity(mt, D, 0)

    ncells = length(D_0.offsets) - 1
    combs = combinations(1:Int(D+1), Int(d+1))

    V = [Array{Array{Int,1},1}() for i in 1:ncells]
    for (i, D_0_i) in enumerate(D_0)
        for comb in combs
            push!(V[i], [D_0_i[j] for j in comb])
        end
    end

    return V
end

