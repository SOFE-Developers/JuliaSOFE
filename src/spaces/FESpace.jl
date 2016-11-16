__precompile__()

module Spaces

import ..Elements: nDoF

using ..Elements
using ..Meshes

export AbstractFESpace, FESpace, MixedFESpace
export mesh, element, domain, domain!, shift, shift!
export dofMap, nDoF, dofIndices, dofMask, extractDoFs, fixedDoF, freeDoF
export interpolate

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
    return extractDoF(fes, d=dim, mask=bmask)
end

"""

    freeDoF(fes::FESpace)

  Return a boolean mask marking the unconstrained
  degrees of freedom of the finite element space.
"""
freeDoF(fes::FESpace) = !fixedDoF(fes)

"""

    dofMap(fes::FESpace, d::Integer)

  Return the degrees of freedom mapping that connects the global mesh
  entities of topological dimension `d` to the local reference element.

  Establishes the connection between the local and global degrees of freedom
  via a connectivity array `C` where `C[i,j] = k` connects the `i`-th
  local basis function  on the reference element to the `k`-th global basis 
  function on the `j`-th element.

  # Keyword Arguments
  * `d::Integer`: The topological dimension of the entities
                    for which to compute the dof map
  * `mask::Vector{T<:Integer}`: A mask marking specific entities
                    for which to compute the dof map
"""
function dofMap{T<:Integer}(fes::FESpace;
                d::Integer = dimension(mesh(fes)),
                mask::AbstractArray{T,1} = 1:number(mesh(fes), d))
    dofTuple = Elements.dofTuple(element(fes))
    dofPerDim = nDoF(element(fes), d)
    nEntities = [getNumber(fes.mesh.topology, dd) for dd = 0:d]
    dofsNeeded = [nEntities[dd+1] * dofTuple[dd+1] for dd = 0:d]
    ndofs = [0, cumsum(dofsNeeded)...]
    dofs = [reshape(ndofs[i]+1:ndofs[i+1], dofTuple[i], nEntities[i]) for i = 1:d+1]

    M = [zeros(Int, dofPerDim[i], nEntities[d+1]) for i = 1:d+1]

    # first, iterate over subdims
    for dd = 0:d-1
        d_dd = getConnectivity(topology(mesh(fes)), d, dd)
        for i = 1:size(d_dd, 1)
            for j = 1:size(d_dd, 2)
                r = (j-1)*dofTuple[dd+1]+1 : j*dofTuple[dd+1]
                M[dd+1][r,i] = dofs[dd+1][:,d_dd[i,j]]
            end
        end
    end

    # set dofs of query dim
    M[d+1] = dofs[d+1]

    return vcat(M...)[:,mask]
end

"""

    nDoF(fes::FESpace)

  Return the total number of degrees of freedom for the
  finite element space `fes`.
"""
function nDoF(fes::FESpace)
    return maxabs(dofMap(fes, d=dimension(mesh(fes))))
end

"""

    dofIndices(fes::FESpace; kwargs...)

  Return the indices of the degrees of freedom for the finite
  element space `fes` associated with the mesh entities of 
  topological dimension `d`.

  # Keyword Arguments
  * `d::Integer`: The topological dimension of the entities
                    for which to compute the dof indices
  * `mask::Vector{T<:Integer}`: A mask marking specific entities
                    for which to compute the dof indices
"""
function dofIndices(fes::FESpace; kwargs...)
    M = dofMap(fes; kwargs...)
    return sort!(unique(M))
end

"""

    dofMask(fes::FESpace; kwargs...)

  Return a boolean mask specifying the degrees of freedom
  for the finite element space `fes` associated with the 
  mesh entities of topological dimension `d`.

  # Keyword Arguments
  * `d::Integer`: The topological dimension of the entities
                    for which to compute the dof mask
  * `mask::Vector{T<:Integer}`: A mask marking specific entities
                    for which to compute the dof mask
"""
function dofMask(fes::FESpace; kwargs...)
    mask = zeros(Bool, nDoF(fes))
    for i in dofIndices(fes; kwargs...)
        mask[i] = true
    end
    return mask
end

extractDoF(fes::FESpace; kwargs...) = dofMask(fes; kwargs...)

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

#-------------------#
# Type MixedFESpace #
#-------------------#
type MixedFESpace <: AbstractFESpace
    subspaces :: Array{FESpace, 1}
end
MixedFESpace(spaces::FESpace...) = MixedFESpace([space for space in spaces])

"""

    subspaces(mfes::MixedFESpace)

  Return the each subspace of the mixed finite element space.
"""
subspaces(mfes::MixedFESpace) = getfield(mfes, :subspaces)
subspace(mfes::MixedFESpace, i::Integer) = subspaces(mfes)[i]

fixedDoF(mfes::MixedFESpace) = mapreduce(fixedDoF, vcat, subspaces(mfes))
freeDoF(mfes::MixedFESpace) = mapreduce(freeDoF, vcat, subspaces(mfes))
nDoF(mfes::MixedFESpace) = mapreduce(nDoF, +, subspaces(mfes))

fixedDoF(mfes::MixedFESpace, i::Integer) = fixedDoF(subspace(mfes, i))
freeDoF(mfes::MixedFESpace, i::Integer) = freeDoF(freeDoF(subspace(mfes, i)))
nDoF(mfes::MixedFESpace, i::Integer) = nDoF(subspace(mfes, i))

Base.getindex(mfes::MixedFESpace, i::Integer) = subspace(mfes, i)
Base.length(mfes::MixedFESpace) = length(subspaces(mfes))
Base.start(::MixedFESpace) = 1
Base.next(mfes::MixedFESpace, state::Integer) = (subspace(mfes, state), state+1)
Base.done(mfes::MixedFESpace, state::Integer) = state > length(mfes)
Base.eltype(::Type{MixedFESpace}) = FESpace
Base.endof(mfes::MixedFESpace) = length(mfes)

end # of module Spaces
