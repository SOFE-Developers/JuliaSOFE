__precompile__()

module Meshes

export AbstractMesh, Mesh, TensorProductMesh, UnitSquare, UnitCube
export evalReferenceMaps, evalJacobianInverse, evalJacobianDeterminat

using ..Elements

include("topology.jl")
using .Topology

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
    element
    topology :: AbstractMeshTopology

    function Mesh{T<:Float, S<:Integer}(nodes::AbstractArray{T,2},
                                        cells::AbstractArray{S,2})
        dim = size(nodes, 2)
        if size(cells, 2) == 3
            element = LagrangeP1(dim)
        elseif size(cells, 2) == 4
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

# Reference Maps
include("refmaps.jl")

# Data Evaluation
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
