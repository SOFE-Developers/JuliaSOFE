
using CMesh: Mesh, nodes, topology, entities
using CMesh.Topology: Simplex, Segment, Triangle, Quadrilateral, Tetrahedron, Hexahedron
using CMesh.Meshes

using CMesh.Topology: AbstractEntity, Simplex, Orthotope,
                      Vertex, Segment, Triangle, Quadrilateral, Tetrahedron, Hexahedron,
                      subtype, dimension, nvertices
using CMesh.Meshes: MeshEntity, mesh, index!, coord, entity

using ..Elements

export evalReferenceMaps, evalInverseMaps, evalJacobianInverse, evalJacobianDeterminat

typealias Float AbstractFloat

"""

    evalReferenceMaps{T<:Float}(m::Mesh, points::AbstractArray{T,2}, deriv::Integer=0)

  Evaluate for each mesh element its associated reference map
  at given local points on the reference domain.
"""
function evalReferenceMaps{T<:Float}(m::Mesh, points::AbstractArray{T,2}, deriv::Integer=0)
    nP, nD = size(points)
    nE = number(m, nD)
    nW = dimension(m)

    p = zeros(T, nD)
    if deriv == 0
        r = zeros(T, nW)
        R = zeros(T, nE, nP, nW)
    elseif deriv == 1
        r = zeros(T, nW, nD)
        R = zeros(T, nE, nP, nW, nD)
    end
    
    itr = MeshEntityIterator(m, nD)
    start!(itr)
    while !done(itr)
        next!(itr)
        me = entity(itr)
        ie = index(me)
        
        for ip = 1:nP
            for id = 1:nD
                p[id] = points[ip,id]
            end

            if deriv == 0
                refmap!(r, p, me)

                for iw = 1:nW
                    R[ie,ip,iw] = r[iw]
                end
            elseif deriv == 1
                jacobian!(r, p, me)

                for id = 1:nD
                    for iw = 1:nW
                        R[ie,ip,iw,id] = r[iw,id]
                    end
                end
            end
        end
    end

    return R
end

"""

    evalInverseMaps{E<:Simplex,T<:Real,S<:Integer}(m::Mesh{E}, points::AbstractArray{T,2}, hosts::AbstractArray{S,1})

  Compute the preimage for each of the given (global) `points` by
  evaluating the inverse of the reference map of the
  corresponding mesh cell specified in `hosts`.
"""
function evalInverseMaps{E<:Simplex,T<:Real,S<:Integer}(m::Mesh{E}, points::AbstractArray{T,2}, hosts::AbstractVector{S})
    @assert size(points, 2) >= dimension(topology(m))
    
    preimg = zeros(T, size(points, 1), dimension(topology(m)))
    q = zeros(T, dimension(topology(m)))
    
    n = nodes(m)
    M = zeros(T, size(points, 2), dimension(topology(m)))
    p = zeros(T, size(points, 2))
    p0 = zeros(T, size(points, 2))
    
    #it = iter(incidence!(topology(m), dimension(topology(m)), 0))
    e = MeshEntity(m, dimension(topology(m)), 1)
    
    for ip = 1:size(points, 1)
        #e = it[hosts[ip]]
        index!(e, hosts[ip])
        for id = 1:size(points, 2)
            p[id] = points[ip,id]
            p0[id] = n[entity(e, 0, 1),id]
            for iv = 2:dimension(topology(m))+1
                M[id,iv-1] = n[entity(e, 0, iv),id] - n[entity(e, 0, 1),id]
            end
        end

        q = pinv(M) * (p - p0)

        for id = 1:dimension(topology(m))
            preimg[ip,id] = q[id]
        end
    end

    return preimg
end

"""

    evalJacobianInverse{T<:Float}(m::Mesh, points::AbstractArray{T,2})

  Compute the inverse of the jacobian of each reference map 
  evaluated in every given local point.
"""
function evalJacobianInverse{T<:Float}(m::Mesh, points::AbstractArray{T,2})
    jacs = evalReferenceMaps(m, points, 1)
    nE, nP, nW, nD = size(jacs)

    if (nW == nD) || (nD == 1)
        fill_JacsInv!(jacs)
    else
        error("Invalid shape of jacobians! (", nW, nD, ")")
    end

    return jacs
end

function fill_JacsInv!{T<:Float}(jacs::Array{T,4})
    issquare = (size(jacs,3) == size(jacs,4))
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

    D = zeros(eltype(jacs), nE, nP)
    fill_JacsDet!(D, jacs, Val{nW}, Val{nD})

    return D
end

function fill_JacsDet!{T<:Float}(D::Array{T,2}, jacs::Array{T,4}, ::Type{Val{2}}, ::Type{Val{1}})
    for ip = 1:size(D,2) # nP
        for ie = 1:size(D,1) # nE
            D[ie,ip] = sqrt(jacs[ie,ip,1,1]^2 + jacs[ie,ip,2,1]^2)
        end
    end
end

function fill_JacsDet!{T<:Float}(D::Array{T,2}, jacs::Array{T,4}, ::Type{Val{3}}, ::Type{Val{1}})
    for ip = 1:size(D,2) # nP
        for ie = 1:size(D,1) # nE
            D[ie,ip] = sqrt(jacs[ie,ip,1,1]^2 + jacs[ie,ip,2,1]^2 + jacs[ie,ip,3,1]^2)
        end
    end
end

function fill_JacsDet!{T<:Float}(D::Array{T,2}, jacs::Array{T,4}, ::Type{Val{2}}, ::Type{Val{2}})
    for ip = 1:size(D,2) # nP
        for ie = 1:size(D,1) # nE
            D[ie,ip] = jacs[ie,ip,1,1] * jacs[ie,ip,2,2] - jacs[ie,ip,2,1] * jacs[ie,ip,1,2]
        end
    end
end

function fill_JacsDet!{T<:Float}(D::Array{T,2}, jacs::Array{T,4}, ::Type{Val{3}}, ::Type{Val{2}})
    for ip = 1:size(D,2) # nP
        for ie = 1:size(D,1) # nE
            D[ie,ip] = sqrt((jacs[ie,ip,2,1] * jacs[ie,ip,3,2] - jacs[ie,ip,3,1] * jacs[ie,ip,2,2])^2
                            + (jacs[ie,ip,3,1] * jacs[ie,ip,1,2] - jacs[ie,ip,1,1] * jacs[ie,ip,3,2])^2
                            + (jacs[ie,ip,1,1] * jacs[ie,ip,2,2] - jacs[ie,ip,2,1] * jacs[ie,ip,1,2])^2)
        end
    end
end

function fill_JacsDet!{T<:Float}(D::Array{T,2}, jacs::Array{T,4}, ::Type{Val{3}}, ::Type{Val{3}})
    for ip = 1:size(D,2) # nP
        for ie = 1:size(D,1) # nE
            D[ie,ip] = jacs[ie,ip,1,1] * jacs[ie,ip,2,2] * jacs[ie,ip,3,3] +
                jacs[ie,ip,1,2] * jacs[ie,ip,2,3] * jacs[ie,ip,3,1] +
                jacs[ie,ip,1,3] * jacs[ie,ip,2,1] * jacs[ie,ip,3,2] -
                jacs[ie,ip,1,1] * jacs[ie,ip,2,3] * jacs[ie,ip,3,2] -
                jacs[ie,ip,1,2] * jacs[ie,ip,2,1] * jacs[ie,ip,3,3] -
                jacs[ie,ip,1,3] * jacs[ie,ip,2,2] * jacs[ie,ip,3,1]
        end
    end
end

######################################################################################################
######################################################################################################

"""

    refmap{T<:Real,E<:AbstractEntity}(p::AbstractVector{T}, me::MeshEntity{E})

Evaluate the reference map corresponding to the mesh
entity `me` in the local point `p`, i.e.
map the point `p` from the reference domain
to the global mesh entity `me`.
"""
function refmap{T<:Real,E<:Simplex}(p::AbstractVector{T}, me::MeshEntity{E})
    return refmap!(zeros(T, dimension(mesh(me))), p, me)
end
function refmap!{T<:Real,E<:Simplex}(q::Vector{T}, p::AbstractVector{T}, me::MeshEntity{E})
    dimw = dimension(mesh(me))
    dimp = length(p); @assert dimp == dimension(me)

    s = sum(p)
    for id = 1:dimw
        q[id] = (1 - s) * coord(me, 1, id)
        for iv = 2:nvertices(E, dimp)
            q[id] += p[iv-1] * coord(me, iv, id)
        end
    end
    
    return q
end

"""

    jacobian{T<:Real,E<:AbstractEntity}(p::AbstractVector{T}, me::MeshEntity{E})

Evaluate the jacobian of the reference map corresponding 
to the mesh entity `me` in the local point `p`.
"""
function jacobian{T<:Real,E<:Simplex}(p::AbstractVector{T}, me::MeshEntity{E})
    return jacobian!(zeros(T, dimension(mesh(me)), dimension(me)), p, me)
end
function jacobian!{T<:Real,E<:Simplex}(q::Array{T,2}, p::AbstractVector{T}, me::MeshEntity{E})
    dimw = dimension(mesh(me))
    dimp = length(p); @assert dimp == dimension(me)
    
    for id = 1:dimw
        for jd = 1:dimp
            q[id,jd] = -coord(me, 1, id) # ∂iφ1 == -1
        end
        for iv = 2:nvertices(E, dimp)
            q[id,iv-1] += coord(me, iv, id) # ∂iφj == δij
        end
    end

    return q
end
