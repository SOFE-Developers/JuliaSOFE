__precompile__()

module Meshes

export AbstractMesh, Mesh
export dimension, topology, nodes, entities
export evalReferenceMaps, evalJacobianInverse, evalJacobianDeterminat

using ..Elements

include("topology.jl")
using .Topology
export getNodes, getConnectivity, getEntities, getNumber

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

# Reference Maps
include("refmaps.jl")

function evaluate{T<:Float}(m::Mesh, f::Function, points::AbstractArray{T,2})
    
end

function evalFunction{T<:Float}(m::Mesh, f::Function, points::AbstractArray{T,2})
    P = evalReferenceMap(m, points) # nExnPxnW
    (nE, nP, nW) = size(P)
    P = reshape(P, nE*nP, nW) # (nE*nP)xnW
    R = f(P) # (nE*nP)x[...]
    R = reshape(R, nE, nP, size(R,2))
end

# Mesh Generation
include("generation.jl")

end # of module Meshes
