module Spaces

using CMesh.Meshes

using ..Elements

import ..Elements: nDoF
import ..CMeshExtensions: evaluate

export AbstractFESpace, FESpace, MixedFESpace
export mesh, element, domain, domain!, shift, shift!
export dofMap, nDoF, dofIndices, dofMask, extractDoFs, fixedDoF, freeDoF
export interpolate, evaluate

typealias Float AbstractFloat

abstract AbstractFESpace

#--------------#
# Type FESpace #
#--------------#
type FESpace{Tm<:AbstractMesh, Te<:AbstractElement} <: AbstractFESpace
    mesh :: Tm
    element :: Te
    domain :: Function
    shift :: Function
end

function FESpace{Tm<:AbstractMesh,Te<:AbstractElement}(mesh::Tm, element::Te)
    domain = x -> falses(size(x,1))
    shift = x -> zeros(size(x,1))
    return FESpace{Tm,Te}(mesh, element, domain, shift)
end

# Associated Methods
# -------------------
"""

    mesh(fes::FESpace)

  Return the mesh of the finite element space.
"""
mesh(fes::FESpace) = getfield(fes, :mesh)

"""

    element(fes::FESpace)

  Return the reference element of the finite element space.
"""
element(fes::FESpace) = getfield(fes, :element)

"""

    domain(fes::FESpace)

  Return the domain function of the finite element space
  that specifies the constrained boundary.
"""
domain(fes::FESpace) = getfield(fes, :domain)
domain!(fes::FESpace, f::Function) = setfield!(fes, :domain, f)

"""

    shift(fes::FESpace)

  Return the shift function that specifies a value on 
  the constrained domain of the finite element space.
"""
shift(fes::FESpace) = getfield(fes, :shift)
shift!(fes::FESpace, f::Function) = setfield!(fes, :shift, f)

"""

    fixedDoF(fes::FESpace)

  Return a boolean mask marking the constrained
  degrees of freedom of the finite element space.
"""
function fixedDoF(fes::FESpace)
    dim = dimension(mesh(fes)) - 1
    bmask = boundary(mesh(fes), domain(fes))
    return extractDoF(fes, dim, bmask)
end

"""

    freeDoF(fes::FESpace)

  Return a boolean mask marking the unconstrained
  degrees of freedom of the finite element space.
"""
freeDoF(fes::FESpace) = !fixedDoF(fes)

"""

    dofMap(fes::FESpace, d::Integer, mask::AbstractVector)

  Return the degrees of freedom mapping that connects the global mesh
  entities of topological dimension `d` to the local reference element.

  Establishes the connection between the local and global degrees of freedom
  via a connectivity array `C` where `C[i,j] = k` connects the `i`-th
  local basis function  on the reference element to the `k`-th global basis 
  function on the `j`-th element.

  # Arguments
  * `d::Integer`: The topological dimension of the entities
                    for which to compute the dof map
  * `mask::Vector{T<:Integer}`: A mask marking specific entities
                    for which to compute the dof map
"""
function dofMap{T<:Integer}(fes::FESpace, d::Integer, mask::AbstractArray{T,1})
    dofs = generateDoFs(element(fes), mesh(fes), d)
    ndof = nDoF(element(fes), d)

    #M = zeros(T, sum(ndof), number(mesh(fes), d))
    M = zeros(Int, sum(ndof), number(mesh(fes), d))

    rb = 0
    for dd = 0:d-1
        d_dd = connectivity(topology(mesh(fes)), d, dd)
        ra = rb + 1; rb += ndof[dd+1]
        fill_dofMap!(view(M, ra:rb, :), dofs[dd+1], d_dd)
    end

    for j = 1:size(dofs[d+1], 2)
        for i = 1:size(dofs[d+1], 1)
            M[rb+i,j] = dofs[d+1][i,j]
        end
    end

    return M[:,mask]
end
dofMap(fes::FESpace, d::Integer) = dofMap(fes, d, 1:number(mesh(fes), d))
dofMap(fes::FESpace) = dofMap(fes, dimension(element(fes)))

function fill_dofMap!{T<:Integer}(M::AbstractArray{T,2}, dofs::AbstractArray{T,2}, inc::AbstractArray{T,2})
    ptr = 1
    for i = 1:size(inc, 1)
        for j = 1:size(inc, 2)
            for k = 1:size(dofs, 1)
                M[ptr,i] = dofs[k, inc[i,j]]
                ptr += 1
            end
        end
        ptr = 1
    end
    return nothing
end

function generateDoFs{T<:AbstractElement}(elem::T, mesh::Mesh, d::Integer)
    doftuple = dofTuple(elem)
    nentities = [number(mesh, dd) for dd = 0:d]
    ndofs = map(*, doftuple, nentities)
    dofrange = [1, 1 + cumsum(ndofs)...]
    dofs = [reshape(dofrange[i]:dofrange[i+1]-1, doftuple[i], nentities[i])
            for i = 1:d+1]

    return dofs
end

function assembleDoFs{T<:Integer}(dofs::AbstractArray{T,1}, d::Integer)
end

"""

    nDoF(fes::FESpace)

  Return the total number of degrees of freedom for the
  finite element space `fes`.
"""
function nDoF(fes::FESpace)
    return maxabs(dofMap(fes, dimension(mesh(fes))))
end

"""

    dofIndices(fes::FESpace, args...)

  Return the indices of the degrees of freedom for the finite
  element space `fes` associated with the mesh entities of 
  topological dimension `d`.

  # Arguments
  * `d::Integer`: The topological dimension of the entities
                    for which to compute the dof indices
  * `mask::Vector{T<:Integer}`: A mask marking specific entities
                    for which to compute the dof indices
"""
function dofIndices(fes::FESpace, args...)
    M = dofMap(fes, args...)
    return sort!(unique(M))
end

"""

    dofMask(fes::FESpace, args...)

  Return a boolean mask specifying the degrees of freedom
  for the finite element space `fes` associated with the 
  mesh entities of topological dimension `d`.

  # Arguments
  * `d::Integer`: The topological dimension of the entities
                    for which to compute the dof mask
  * `mask::Vector{T<:Integer}`: A mask marking specific entities
                    for which to compute the dof mask
"""
function dofMask(fes::FESpace, args...)
    mask = zeros(Bool, nDoF(fes))
    for i in dofIndices(fes, args...)
        mask[i] = true
    end
    return mask
end

extractDoF(fes::FESpace, args...) = dofMask(fes, args...)

function interpolate(m::Mesh, el::Element, f::Function)
    @assert isnodal(el)

    d = dimension(el)
    p = order(el)

    v = nodes(m)
    e = p > 1 ? evalReferenceMap(m, linspace(0,1,p+1)[2:end-1]) : zeros(0,d)
    i = (p > 2 && d > 1) ? evalReferenceMap(m, lagrangeNodesP(d,p)[p*(d+1):end,:]) : zeros(0,d)

    n = vcat(v, e, i)

    return f(n)
end
interpolate(fes::FESpace, f::Function) = interpolate(fes.mesh, fes.element, f)

"""

    evaluate{T<:AbstractFloat}(fes::FESpace, dofs::AbstractVector{T}, deriv::Integer=0)

  Evaluate the linear combination of the finite element
  space's basis functions or their derivatives in the 
  given local `points` w.r.t. the given `dof` values.
"""
function evaluate{T<:Float,S<:Float}(fes::FESpace, dofs::AbstractVector{S},
                                     points::AbstractArray{T,2}, deriv::Integer=0)
    nP, dimP = size(points)
    dofmap = dofMap(fes, dimP)
    nB, nE = size(dofmap)
    basis = evalBasis(element(fes), points, deriv)
    nB, nP, nC, nD = size(basis, 1:4...)

    dofs = convert(Vector{promote_type(T, S)}, dofs)
    
    if deriv == 0
        @assert nD == 1
        U = zeros(T,nE,nP,nC)
    elseif deriv == 1
        U = zeros(T,nE,nP,nC,nD)
    elseif deriv == 2
        U = zeros(T,nE,nP,nC,nD,nD)
    end

    fill_dofvalues!(U, basis, dofs, dofmap)

    return U
end

function fill_dofvalues!{T<:AbstractFloat}(U::AbstractArray{T,3}, basis::Array{T,3},
                                   dofs::AbstractVector{T}, dofmap::Array{Int,2})
    nE, nP, nC = size(U)
    nB, nP, nC = size(basis)
    nB, nE = size(dofmap)
    
    for ic = 1:nC
        for ip = 1:nP
            for ie = 1:nE
                for ib = 1:nB
                    U[ie,ip,ic] += dofs[dofmap[ib,ie]] * basis[ib,ip,ic]
                end
            end
        end
    end

    return nothing
end

function fill_dofvalues!{T<:AbstractFloat}(U::Array{T,4}, basis::Array{T,4},
                                   dofs::AbstractVector{T}, dofmap::Array{Int,2})
    nE, nP, nC, nD = size(U)
    nB, nP, nC, nD = size(basis)
    nB, nE = size(dofmap)
    
    for id = 1:nD
        fill_dofvalues!(view(U, :, :, :, id), basis[:,:,:,id], dofs, dofmap)
    end

    return nothing
end

function fill_dofvalues!{T<:AbstractFloat}(U::Array{T,5}, basis::Array{T,5},
                                   dofs::AbstractVector{T}, dofmap::Array{Int,2})
    nE, nP, nC, nDi, nDj = size(U)
    nB, nP, nC, nDi, nDj = size(basis)
    nB, nE = size(dofmap)

    for jd = 1:nDj
        for id = 1:nDi
            fill_dofvalues!(view(U, :, :, :, id, jd), basis[:,:,:,id,jd], dofs, dofmap)
        end
    end

    return nothing
end

# Mixed finite element spaces
# ----------------------------
include("mixed.jl")

end # of module Spaces
