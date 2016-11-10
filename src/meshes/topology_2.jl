import Combinatorics: combinations

export MeshTopology
export incidence!, incidence, connectivity!, connectivity
#export init!,  build!, build, transpose!, transpose, intersection!, intersection

include("connectivity.jl")

typealias Float AbstractFloat

type Simp <: AbstractMeshTopology end
type Orth <: AbstractMeshTopology end

typealias TopologyTypes Union{Simp, Orth}

#-------------------#
# Type MeshTopology #
#-------------------#
"""
  Stores the topology of a mesh as a set of incidence relations.

  `MeshTopology(dim, nodes, cells)` constructs and initializes 
  a mesh topology instance.
"""
type MeshTopology{T<:TopologyTypes,Tn<:Float,Ti<:Integer} <: AbstractMeshTopology
    dimension :: Int
    nodes :: Array{Tn, 2}
    connectivities :: Dict{Tuple{Int, Int}, MeshConnectivity{Ti}}
end

# Outer Constructors
# -------------------
function MeshTopology{T<:TopologyTypes,Tn<:Float,Ti<:Integer}(::Type{T},
                                                              nodes::AbstractArray{Tn,2},
                                                              cells::AbstractArray{Ti,2})
    D = dimension(T, size(cells, 2))
    connectivities = Dict{Tuple{Int,Int}, MeshConnectivity{Ti}}()
    mt = MeshTopology{T,Tn,Ti}(D, nodes, connectivities)
    init!(mt, cells)
    return mt
end

# Associated Methods
# -------------------
function Base.show(io::IO, mt::MeshTopology)
    println(io, "MeshTopology")
    println(io, "\tdimension : ", mt.dimension)
    println(io, "\tconnectivities : ",
            ["\n\t\t$connect" for connect in mt.connectivities]...)
end

function getConnectivity(mt::MeshTopology, d::Integer, dd::Integer)
    haskey(mt.connectivities, (d,dd)) || incidence!(mt, d, dd)
    return mt.connectivities[(d,dd)][:]
end

function getEntities(mt::MeshTopology, d::Integer)
    if d == 0
        return hcat(getConnectivity(mt, d, 0)...)'
    else
        return getConnectivity(mt, d, 0)
    end
end

connectivity!(mt::MeshTopology, d::Integer, dd::Integer) = incidence!(mt, d, dd)[:]
connectivity(mt::MeshTopology, d::Integer, dd::Integer) = incidence(mt, d, dd)[:]

entities!(mt::MeshTopology, d::Integer) = connectivity!(mt, d, 0)
entities(mt::MeshTopology, d::Integer) = connectivity(mt, d, 0)

"""

    boundary(mt::MeshTopology)

  Return a boolean mask marking the boundary facets of the mesh.
"""
boundary(mt::MeshTopology) = boundary(mt::MeshTopology, dimension(mt)-1)
getBoundary(mt::MeshTopology) = boundary(mt)

"""

    boundary(mt::MeshTopology, d::Integer)

  Return a boolean mask marking the boundary entities of
  topological dimension `d`.
"""
function boundary(mt::MeshTopology, d::Integer)
    D = dimension(mt)
    facets_D = incidence(mt, D-1, D)

    mask = zeros(Bool, size(facets_D, 1))
    for i = 1:size(facets_D, 1)
        mask[i] = length(facets_D[i]) == 1
    end

    if d != D-1
        d_facets = incidence(mt, d, D-1)
        ind_facets = find(mask)

        mask = zeros(Bool, size(d_facets, 1))
        for j in eachindex(mask)
            mask[j] = any(i in ind_facets for i in d_facets[j])
        end
    end

    return mask
end

"""

    boundary(mt::MeshTopology, f::Function)

  Return a boolean mask marking the boundary facets whos
  midpoints belong to the boundary part specified by `f`.
"""
function boundary(mt::MeshTopology, f::Function)
    mask = boundary(mt)
    facets = entities(mt, dimension(mt)-1)
    centroids = mean(nodes(mt)[facets,:], 2)[:,1,:]
    mask &= f(centroids)
    return mask
end
getBoundary(mt::MeshTopology, f::Function) = boundary(mt, f)

# Initialization and incidence computation
# -----------------------------------------

"""

    init!{T<:Integer}(mt::MeshTopology, cells::AbstractArray{T,2})

  Compute the incidence relations `D -> 0` from the given cell 
  connectivity array where `D` is the dimension of the topology.
"""
function init!{T<:Integer}(mt::MeshTopology, cells::AbstractArray{T,2})
    D = dimension(mt)

    ncells, nverts = size(cells)
    
    ind_D_0 = cells'[:]
    off_D_0 = collect(1:nverts:nverts*(ncells+1))
    
    mt.connectivities[(D,0)] = MeshConnectivity(D, 0, ind_D_0, off_D_0)
end

"""

    incidence(mt::MeshTopology, d::Integer, dd::Integer)

  Compute the incidence relation `d -> dd` by successive application
  of `build`, `transpose` and `intersection`.
"""
function incidence!(mt::MeshTopology, d::Integer, dd::Integer)
    # d  != 0 && (haskey(mt.connectivities, (d, 0)) || build!(mt, d))
    # dd != 0 && (haskey(mt.connectivities, (dd,0)) || build!(mt, dd))
    haskey(mt.connectivities, (d, 0)) || build!(mt, d)
    haskey(mt.connectivities, (dd,0)) || build!(mt, dd)
    
    if !haskey(mt.connectivities, (d,dd))
        if d < dd
            incidence!(mt, dd, d)
            transpose!(mt, dd, d)
        else
            D = mt.dimension
            ddd = (d == 0 && dd == 0) ? D : 0
            
            incidence!(mt, d, ddd)
            incidence!(mt, ddd, dd)
            intersection!(mt, d, dd, ddd)
        end
    end

    return mt.connectivities[(d,dd)]
end
function incidence(mt::MeshTopology, d::Integer, dd::Integer)
    prekeys = collect(keys(mt.connectivities))
    d_dd = incidence!(mt, d, dd)
    postkeys = keys(mt.connectivities)
    for key in setdiff(postkeys, prekeys)
        delete!(mt.connectivities, key)
    end
    return d_dd
end

"""

    build!(mt::MeshTopology, d::Integer)

  Compute the incidence relation `d -> 0` from `D -> 0`
  for `0 < d < D` where `D` is the dimension of the topology.
"""
function build!(mt::MeshTopology, d::Integer)
    D = mt.dimension

    if d == 0
        nnodes = size(mt.nodes, 1)
        ind_d_0 = collect(1:nnodes)
        off_d_0 = collect(1:nnodes+1)
    elseif 0 < d < D
        # compute the set of vertex sets for each cell
        V = vertex_sets(mt, d)
        
        # take unique sorted index sets
        entities_d = unique(sort!.(vcat(V...)))
        sort!(entities_d, lt=lexless)

        nentities_d = length(entities_d)
        nverts_d = nvertices(mt, d)
        
        # compute indices and offsets vectors
        ind_d_0 = vcat(entities_d...)
        off_d_0 = collect(1:nverts_d:nverts_d*(nentities_d+1))
    elseif d == D
        warn("Cannot build cell dimension, use init! method!")
        return nothing
    end
    
    mt.connectivities[(d,0)] = MeshConnectivity(d, 0, ind_d_0, off_d_0)
end
function build(mt::MeshTopology, d::Integer)
    prekeys = collect(keys(mt.connectivities))
    d_0 = build!(mt, d)
    postkeys = keys(mt.connectivities)
    for key in setdiff(postkeys, prekeys)
        delete!(mt.connectivities, key)
    end
    return d_0
end

"""

    transpose!(mt::MeshTopology, d::Integer, dd::Integer)

  Compute the incidence relation `dd -> d` from `d -> dd`
  for `d > dd`.
"""
function transpose!(mt::MeshTopology, d::Integer, dd::Integer)
    @assert d > dd
    # get connectivity d -> dd
    d_dd = incidence!(mt, d, dd)

    I, J = findn(d_dd)

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

function transpose(mt::MeshTopology, d::Integer, dd::Integer)
    prekeys = collect(keys(mt.connectivities))
    dd_d = transpose!(mt, d, dd)
    postkeys = keys(mt.connectivities)
    for key in setdiff(postkeys, prekeys)
        delete!(mt.connectivities, key)
    end
    return dd_d
end

"""

    intersection!(mt::MeshTopology, d::Integer, dd::Integer, ddd::Integer)

  Compute the incidence relation `d -> dd` from
  `d -> ddd` and `ddd -> dd` for `d >= dd`.
"""
function intersection!(mt::MeshTopology, d::Integer, dd::Integer, ddd::Integer)
    @assert d >= dd
    d_ddd  = incidence!(mt, d,   ddd)
    ddd_dd = incidence!(mt, ddd, dd)

    nd = size(d_ddd, 1)
    ndd = size(ddd_dd, 2)
    size(d_ddd, 2) == size(ddd_dd, 1) || throw(DimensionMismatch())

    row = zeros(Int, ndd)
    nz = zeros(Bool, ndd)

    ind_d_dd = [Array{Int,1}() for i = 1:nd]

    nverts = (d == dd ? nvertices(mt, dd-1) : nvertices(mt, dd))
    
    for i = 1:nd
        fill_row!(row, nz, d_ddd[i], ddd_dd)
        # for j in d_ddd[i]
        #     for k in ddd_dd[j]
        #         row[k] += 1
        #         nz[k] || (nz[k] = true)
        #     end
        # end

        fill_ind!(ind_d_dd[i], row, nz, nverts)
        # for k = 1:length(row)
        #     if nz[k] && row[k] == nverts
        #         resize!(ind_d_dd[i], length(ind_d_dd[i])+1)
        #         ind_d_dd[i][end] = k
        #     end
        #     row[k] = 0
        #     nz[k] = false
        # end
    end

    nind = [length(i) for i in ind_d_dd]

    ind_d_dd = vcat(ind_d_dd...)
    off_d_dd = vcat(1, 1 + cumsum!(nind, nind))
    
    mt.connectivities[(d,dd)] = MeshConnectivity(d, dd, ind_d_dd, off_d_dd)
end
function intersection(mt::MeshTopology, d::Integer, dd::Integer, ddd::Integer)
    prekeys = collect(keys(mt.connectivities))
    d_dd = intersection!(mt, d, dd, ddd)
    postkeys = keys(mt.connectivities)
    for key in setdiff(postkeys, prekeys)
        delete!(mt.connectivities, key)
    end
    return d_dd
end

@noinline function fill_row!(row::Array{Int,1}, nz::Array{Bool,1}, d_ddd_i::Array{Int,1}, ddd_dd::MeshConnectivity)
    for j in d_ddd_i
        for k in ddd_dd[j]
            row[k] += 1
            nz[k] || (nz[k] = true)
        end
    end
end

@noinline function fill_ind!(ind::Array{Int,1}, row::Array{Int,1}, nz::Array{Bool,1}, nverts::Int)
    for k = 1:length(row)
        if nz[k] && row[k] == nverts # (d == dd ? dd : dd+1)
            resize!(ind, length(ind)+1)
            ind[end] = k
        end
        row[k] = 0
        nz[k] = false
    end
end

"""

    vertex_sets(mt::MeshTopology, d::Integer)

  Compute for each cell the set of vertex sets 
  incident to the entities of topological dimension `d`.
"""
function vertex_sets(mt::MeshTopology, d::Integer)
    D = mt.dimension
    D_0 = getConnectivity(mt, D, 0)

    ncells = size(D_0, 1)
    combs = vertex_combs(mt, d)
    ncombs = size(combs, 1)

    V = [[D_0[i,comb] for comb in combs] for i = 1:ncells]
    return V
end

"""

    vertex_combs(mt::MeshTopology, d::Integer)

  Return the local vertex combinations that define the
  `d`-dimensional subentities of the mesh cells. 
"""
function vertex_combs(mt::MeshTopology, d::Integer) end

"""

    nvertices(mt::MeshTopology, d::Integer)

  Return the number of vertices that define a `d`-dimensional
  entity of the mesh.
"""
function nvertices(mt::MeshTopology, d::Integer) end

"""

    dimension{T<:AbstractMeshTopology,Ti<:Integer}(::Type{T}, cells::AbstractArray{T,2})

  Determine the topological dimension of the given `cells`
  according to the given topology type `T`.
"""
function dimension{T<:TopologyTypes,Ti<:Integer}(::Type{T}, cells::AbstractArray{Ti,2})
    return dimension(T, cells)
end

#-------------------------#
# Simplex Types Specifics #
#-------------------------#
typealias MeshTopologySimp MeshTopology{Simp}

nvertices(::MeshTopology{Simp}, d::Integer) = d+1

dimension(::Type{Simp}, nverts::Integer) = nverts - 1

function vertex_combs(mt::MeshTopology{Simp}, d::Integer)
    return collect(combinations(1:dimension(mt)+1, nvertices(mt, d)))
end

#---------------------------#
# Orthotope Types Specifics #
#---------------------------#
typealias MeshTopologyOrth MeshTopology{Orth}

nvertices(::MeshTopology{Orth}, d::Integer) = 2^d

dimension(::Type{Orth}, nverts::Integer) = ispow2(nverts) ? Int(log2(nverts)) : NaN

function vertex_combs(mt::MeshTopology{Orth}, d::Integer)
    D = dimension(mt)
    if D == 1
        d == 0 && return [[1], [2]]
        d == 1 && return [[1,2]]
    elseif D == 2
        d == 0 && return [[1], [2], [3], [4]]
        d == 1 && return [[1,2], [1,3], [2,4], [3,4]]
        d == 2 && return [[1,2,3,4]]
    elseif D == 3
        d == 0 && return [[1], [2], [3], [4], [5], [6], [7], [8]]
        d == 1 && return [[1,2], [1,3], [1,5], [2,4], [2,6], [3,4],
                          [3,7], [4,8], [5,6], [5,7], [6,8], [7,8]]
        d == 2 && return [[1,2,3,4], [1,2,5,6], [2,4,6,8],
                          [3,1,7,5], [4,3,8,7], [5,6,7,8]]
        d == 3 && return [[1,2,3,4,5,6,7,8]]
    else
        return [Int[]]
    end
end

