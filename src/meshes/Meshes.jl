__precompile__()

module Meshes

export Mesh, TensorProductMesh, UnitSquare, UnitCube
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

"""

    evalReferenceMaps{T<:Float}(m::Mesh, points::AbstractArray{T,2}, deriv::Integer=0)

  Evaluate for each mesh element its associated reference map
  at given local points on the reference domain.
"""
function evalReferenceMaps{T<:Float}(m::Mesh, points::AbstractArray{T,2}, deriv::Integer=0)
    nP, nD = size(points)

    nodes = getNodes(m.topology)
    nW = size(nodes, 2)

    basis = evalBasis(m.element, points, deriv)
    nB = size(basis, 1)

    elems = getEntities(m.topology, nD)
    nE = size(elems, 1)
    @assert nB == size(elems, 2)

    if deriv == 0
        R = zeros(eltype(points), nE, nP, nW)
    elseif deriv == 1
        R = zeros(eltype(points), nE, nP, nW, nD)
    elseif deriv == 2
        R = zeros(eltype(points), nE, nP, nW, nD, nD)
    end

    @time fill_RefMaps!(R, nodes, elems, basis)

    return R
end

function fill_RefMaps!{T<:Float}(R::Array{T,3}, nodes::Array{T,2},
                                 elems::Array{Int,2}, basis::Array{T,2})
    for iw = 1:size(R,3) # nW
        for ip = 1:size(R,2) # nP
            for ie = 1:size(R,1) # nE
                for ib = 1:size(basis,1) # nB
                    R[ie,ip,iw] += nodes[elems[ie,ib],iw] * basis[ib,ip]
                end
            end
        end
    end
end

function fill_RefMaps!{T<:Float}(R::Array{T,4}, nodes::Array{T,2},
                                 elems::Array{Int,2}, basis::Array{T,3})
    for id = 1:size(R,4) # nD
        for iw = 1:size(R,3) # nW
            for ip = 1:size(R,2) # nP
                for ie = 1:size(R,1) # nE
                    for ib = 1:size(basis,1) # nB
                        R[ie,ip,iw,id] += nodes[elems[ie,ib],iw] * basis[ib,ip,id]
                    end
                end
            end
        end
    end
end

function fill_RefMaps!{T<:Float}(R::Array{T,5}, nodes::Array{T,2},
                                 elems::Array{Int,2}, basis::Array{T,4})
    for jd = 1:size(R,5) # nD
        for id = 1:size(R,4) # nD
            for iw = 1:size(R,3) # nW
                for ip = 1:size(R,2) # nP
                    for ie = 1:size(R,1) # nE
                        for ib = 1:size(basis,1) # nB
                            R[ie,ip,iw,id,jd] += nodes[elems[ie,ib],iw] * basis[ib,ip,id,jd]
                        end
                    end
                end
            end
        end
    end
end

"""

    evalJacobianInverse{T<:Float}(m::Mesh, points::AbstractArray{T,2})

Compute the inverse of the jacobian of each reference map 
evaluated in every given local point.
"""
function evalJacobianInverse{T<:Float}(m::Mesh, points::AbstractArray{T,2})
    jacs = evalReferenceMaps(m, points, 1)
    nE, nP, nW, nD = size(jacs)

    if nW == nD
        for ip = 1:nP
            for ie = 1:nE
                jacs[ie,ip,:,:] = inv(jacs[ie,ip,:,:])
            end
        end
    elseif nD == 1
        for ip = 1:nP
            for ie = 1:nE
                jacs[ie,ip,:,:] = 1./jacs[ie,ip,:,:]
            end
        end
    else
        error("Invalid shape of jacobians! (", nW, nD, ")")
    end

    return jacs
end

function fill_JacsInv!{T<:Float}(jacs::Array{T,4})
    issquare = (nW == nD)
    for ip = 1:size(jacs,2) # nP
        for ie = 1:size(jacs,1) # nE
            jacs[ie,ip,:,:] = issquare ? inv(jacs[ie,ip,:,:]) : 1./jacs[ie,ip,:,:]
        end
    end
end

"""

    evalJacobianInverse{T<:Float}(m::Mesh, points::AbstractArray{T,2})

Compute the determinant of the jacobian of each reference map 
evaluated in every given local point.
"""
function evalJacobianDeterminat{T<:Float}(m::Mesh, points::AbstractArray{T,2})
    jacs = evalReferenceMaps(m, points, 1)
    nE, nP, nW, nD = size(jacs)

    det_jacs = zeros(eltype(jacs), nE, nP)
    if nW == nD
        for ip = 1:nP
            for ie = 1:nE
                det_jacs[ie,ip] = det(jacs[ie,ip,:,:])
            end
        end
    elseif nW == 2 && nD == 1
        for ip = 1:nP
            for ie = 1:nE
                det_jacs[ie,ip] = norm(jacs[ie,ip,:])
            end
        end
    elseif nW == 3 && nD == 2
        for ip = 1:nP
            for ie = 1:nE
                t1 = (jacs[ie,ip,2,1] * jacs[ie,ip,3,2] - jacs[ie,ip,3,1] * jacs[ie,ip,2,2])^2
                t2 = (jacs[ie,ip,3,1] * jacs[ie,ip,1,2] - jacs[ie,ip,1,1] * jacs[ie,ip,3,2])^2
                t3 = (jacs[ie,ip,1,1] * jacs[ie,ip,2,2] - jacs[ie,ip,2,1] * jacs[ie,ip,1,2])^2
                det_jacs[ie,ip] = sqrt(t1 + t2 + t3)
            end
        end
    else
        error("Invalid shape of jacobians! (", nW, nD, ")")
    end

    return det_jacs
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

# Mesh Generation
include("generation.jl")

end # of module Meshes
