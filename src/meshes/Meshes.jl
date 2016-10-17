__precompile__()

module Meshes

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
            element = P1(dim)
        elseif size(cells, 2) == 4
            element = Q1(dim)
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
function evalReferenceMaps{T<:Float}(m::Mesh, points::AbstractArray{T,2}, deriv::Integer=0)
    nP, nD = size(points)

    nodes = getNodes(m.topology)
    basis = evalBasis(m.element, points, deriv)
    nB = size(basis, 1); nDv = size(basis, 4)

    connect = getEntities(m.topology, nD)
    nE = size(connect, 1)

    R = reshape(nodes[connect[:],:]', obj.dim*nE, nB) *
        reshape(basis, nB, nP*nDv); # (nW*nE)x(nP*[...])
    R = permutedims(reshape(R, obj.dim, nE, nP, nDv), [2 3 1 4]); # nExnPxnWx[...];

    return R
end

function evalTrafo{T<:Float}(m::Mesh, dPhi::AbstractArray{T,4})
    (nW, nD) = size(dPhi, 3, 4);
    if nW == 2
        if nD == 1
            R = sqrt(dPhi[:,:,1].^2 + dPhi[:,:,2].^2);
        elseif nD == 2
            R = dPhi[:,:,1,1].*dPhi[:,:,2,2] - dPhi[:,:,1,2].*dPhi[:,:,2,1];
        end
    elseif nW == 3
        if nD == 1
            R = sqrt(dPhi[:,:,1].^2 + dPhi[:,:,2].^2 + dPhi[:,:,3].^2);
        elseif nD == 2
            R = sqrt((dPhi[:,:,2,1].*dPhi[:,:,3,2] - dPhi[:,:,3,1].*dPhi[:,:,2,2]).^2 +
                     (dPhi[:,:,3,1].*dPhi[:,:,1,2] - dPhi[:,:,1,1].*dPhi[:,:,3,2]).^2 +
                     (dPhi[:,:,1,1].*dPhi[:,:,2,2] - dPhi[:,:,2,1].*dPhi[:,:,1,2]).^2);
        elseif nD == 3
            dPhi = dPhi[:,:,1,1].*dPhi[:,:,2,2].*dPhi[:,:,3,3] +
                dPhi[:,:,1,2].*dPhi[:,:,2,3].*dPhi[:,:,3,1] +
                dPhi[:,:,1,3].*dPhi[:,:,2,1].*dPhi[:,:,3,2] -
                dPhi[:,:,1,1].*dPhi[:,:,2,3].*dPhi[:,:,3,2] -
                dPhi[:,:,1,2].*dPhi[:,:,2,1].*dPhi[:,:,3,3] -
                dPhi[:,:,1,3].*dPhi[:,:,2,2].*dPhi[:,:,3,1];
        end
    end
    return R
end

function evalTrafoPair{T<:Float}(m::Mesh, points::AbstractArray{T,2})
    DPhi = evalReferenceMap(m, points, 1);
    trafo = evalTrafo(m, DPhi);
    return (DPhi, trafo)
end

function evalTrafoTriple{T<:Float}(m::Mesh, points::AbstractArray{T,2})
    R = evalReferenceMap(m, points, 1); # nExnPxnWx[...]
    if m.dimension == 1
        det = R; RInv = 1;
    elseif m.dimension == 2
        det = R[:,:,1,1].*R[:,:,2,2] - R[:,:,1,2].*R[:,:,2,1];
        RInv = -R;
        RInv[:,:,1,1] = R[:,:,2,2];
        RInv[:,:,2,2] = R[:,:,1,1];
    elseif m.dimension == 3
        det = R[:,:,1,1].*R[:,:,2,2].*R[:,:,3,3] +
            R[:,:,1,2].*R[:,:,2,3].*R[:,:,3,1] +
            R[:,:,1,3].*R[:,:,2,1].*R[:,:,3,2] -
            R[:,:,1,1].*R[:,:,2,3].*R[:,:,3,2] -
            R[:,:,1,2].*R[:,:,2,1].*R[:,:,3,3] -
            R[:,:,1,3].*R[:,:,2,2].*R[:,:,3,1];
        RInv[:,:,1,1] = R[:,:,2,2].*R[:,:,3,3] - R[:,:,2,3].*R[:,:,3,2];
        RInv[:,:,2,1] = -(R[:,:,2,1].*R[:,:,3,3] - R[:,:,3,1].*R[:,:,2,3]);
        RInv[:,:,3,1] = R[:,:,2,1].*R[:,:,3,2] - R[:,:,2,2].*R[:,:,3,1];
        RInv[:,:,1,2] = -(R[:,:,1,2].*R[:,:,3,3] - R[:,:,1,3].*R[:,:,3,2]);
        RInv[:,:,2,2] = R[:,:,1,1].*R[:,:,3,3] - R[:,:,1,3].*R[:,:,3,1];
        RInv[:,:,3,2] = -(R[:,:,1,1].*R[:,:,3,2] - R[:,:,1,2].*R[:,:,3,1]);
        RInv[:,:,1,3] = R[:,:,1,2].*R[:,:,2,3] - R[:,:,1,3].*R[:,:,2,2];
        RInv[:,:,2,3] = -(R[:,:,1,1].*R[:,:,2,3] - R[:,:,1,3].*R[:,:,2,1]);
        RInv[:,:,3,3] = R[:,:,1,1].*R[:,:,2,2] - R[:,:,1,2].*R[:,:,2,1];
    else
        error("Non supported dimension")
    end
    RInv = broadcast(/, RInv, det);
    return (R, RInv, det)
end

# Data Evaluation
function evalFunction{T<:Float}(m::Mesh, f::Function, points::AbstractArray{T,2})
    P = evalReferenceMap(m, points) # nExnPxnW
    (nE, nP, nW) = size(P)
    P = reshape(P, nE*nP, nW) # (nE*nP)xnW
    R = f(P) # (nE*nP)x[...]
    R = reshape(R, nE, nP, size(R,2))
end


end # of module Meshes
