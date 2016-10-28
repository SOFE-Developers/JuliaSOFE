__precompile__()

module Operators

using ..Elements
using ..Meshes
using ..Spaces
using ..Quadrature

export AbstractOperator, coeff, space, matrix, vector, assemble!
export Operator, IdId, GradGrad
export Functional, Id

typealias Float AbstractFloat

include("coefficients.jl")

#--------------------#
# Abstract Operators #
#--------------------#
abstract AbstractOperator

# General Operator Methods
"""

    space(op::AbstractOperator)

  Return the finite element space of the operator `op`.
"""
@inline space{T<:AbstractOperator}(op::T) = op.fes
@inline getFESpace{T<:AbstractOperator}(op::T) = space(op)

"""

    coeff(op::AbstractOperator)

  Return the coefficient of the operator `op`.
"""
@inline coeff{T<:AbstractOperator}(op::T) = op.coeff
@inline getCoefficient{T<:AbstractOperator}(op::T) = coeff(op)

# Includes
include("bilinear.jl")
include("linear.jl")

end # of module Operators

