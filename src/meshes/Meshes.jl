__precompile__()

module Meshes

using ..Elements

import ..Elements: dimension

export AbstractMesh, Mesh
export dimension, topology, nodes, entities, number, boundary
export evaluate


include("topology.jl")
using .Topology
export connectivity

#--------------------#
# Abstract Mesh Type #
#--------------------#
abstract AbstractMesh

typealias Float AbstractFloat

#-----------#
# Type Mesh #
#-----------#
type Mesh <: AbstractMesh
    dimension :: Integer
    element :: Element
    topology :: AbstractMeshTopology

    function Mesh{T<:Float, S<:Integer}(nodes::AbstractArray{T,2},
                                        cells::AbstractArray{S,2})
        dim = size(nodes, 2)
        if size(cells, 2) == dim + 1
            element = LagrangeP1(dim)
        elseif (dim == 2) & (size(cells, 2) == 4)
            element = LagrangeQ1(dim)
        else
            error()
        end
        topology = MeshTopology(nodes, cells)

        return new(dim, element, topology)
    end
end

# Associated Methods
# -------------------
"""

    dimension(m::Mesh)

  Return the spatial dimension of the mesh (nodes).
"""
@inline dimension(m::Mesh) = getfield(m, :dimension)

"""

    element(m::Mesh)

  Return the shape element of the mesh.
"""
@inline element(m::Mesh) = getfield(m, :element)

"""

    topology(m::Mesh)

  Return the topology of the mesh.
"""
@inline topology(m::Mesh) = getfield(m, :topology)

"""

    nodes(m::Mesh)

  Return the node coordinates of the mesh.
"""
Topology.nodes(m::Mesh) = nodes(topology(m))

"""

    entities(m::Mesh, d::Integer)

  Return the connectivity array for the mesh entities
  of topological dimension `d`.
"""
Topology.entities(m::Mesh, d::Integer) = entities(topology(m), d)

"""

    number(m::Mesh, d::Integer)

  Return the number of mesh entities
  of topological dimension `d`.
"""
Topology.number(m::Mesh, d::Integer) = number(topology(m), d)

Topology.boundary(m::Mesh, d::Integer) = boundary(topology(m), d)
Topology.boundary(m::Mesh, f::Function) = boundary(topology(m), f)

# Reference Maps
# ---------------
include("refmaps.jl")

# Data Evaluation
# ----------------
typealias Constant{T<:Real} Union{T, AbstractVector{T}, AbstractMatrix{T}}

function evaluate{T<:Constant,S<:Float}(c::T, points::AbstractArray{S,2}, m::Mesh)
    P = evalReferenceMaps(m, points)
    nE, nP, nW = size(P)

    if issubtype(T, Real)
        C = convert(promote_type(eltype(T),S), c)
    elseif issubtype(T, AbstractVector)
        C = convert(Vector{promote_type(eltype(T),S)}, c)
    elseif issubtype(T, AbstractMatrix)
        C = convert(Matrix{promote_type(eltype(T),S)}, c)
    end
    
    R = zeros(S, nE, nP, size(C)...)
    for ip = 1:nP
        for ie = 1:nE
            R[ie,ip,:] = C
        end
    end
    return R
end

function evaluate{T<:Float}(f::Function, points::AbstractArray{T,2}, m::Mesh)
    P = evalReferenceMaps(m, points)
    nE, nP, nW = size(P)

    R = f(reshape(P, (nE*nP,nW)))
    nF = size(R, (2,(3:ndims(R))...)...)
    R = (nF == 1) ? reshape(R, nE, nP) : reshape(R, nE, nP, nF...)

    return ndarray(R)
end

# Mesh Generation
# ----------------
include("generation.jl")

end # of module Meshes
