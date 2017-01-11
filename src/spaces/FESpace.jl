module Spaces

using CMesh.Meshes

using ..Elements
using ..CMeshExtensions: evalJacobianDeterminat

import ..Elements: nDoF

export AbstractFESpace, FESpace, MixedFESpace
export mesh, element, domain, domain!, shift, shift!
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

# DoF Management
# ---------------
include("dofman.jl")

# Interpolation
# --------------
include("interpolation.jl")

# Projection
# -----------
include("projection.jl")

# Evaluation
# -----------
include("evaluation.jl")
# Mixed finite element spaces
# ----------------------------
include("mixed.jl")

end # of module Spaces
