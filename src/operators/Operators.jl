__precompile__()

module Operators

using ..Elements
using ..Meshes
using ..Spaces
using ..Quadrature

export AbstractCoefficient, ConstantCoefficient, ScalarCoefficient,
    VectorCoefficient, MatrixCoefficient, FunctionCoefficient
export value, evaluate
export AbstractOperator, coeff, space, matrix, vector, assemble!
export Operator, IdId, GradGrad
export Functional, Id

typealias Float AbstractFloat

include("coefficients.jl")

#--------------------#
# Abstract Operators #
#--------------------#
abstract AbstractOperator

include("basic.jl")

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

function assemble(op::Operator, d::Integer)
    qpoints, qweights = quadData(op.quadrule, d)
    
    fes = space(op)
    dMap = dofMap(fes, d)
    ndof = Spaces.nDoF(fes)
    
    nB, nE = size(dMap)
    nP, nD = size(qpoints)

    dofI = zeros(Int, nE, nB, nB)
    dofJ = zeros(Int, nE, nB, nB)

    for jb = 1:nB
        for ib = 1:nB
            for ie = 1:nE
                dofI[ie,ib,jb] = dMap[ib,ie]
                dofJ[ie,ib,jb] = dMap[jb,ie]
            end
        end
    end
    
    C, U, V = evaluate(op, d)
    Det = abs(evalJacobianDeterminat(fes.mesh, qpoints))
    
    entries = zeros(eltype(qpoints), nE, nB, nB)
    fill_entries!(op, entries, C, U, V, qweights, Det)

    A = sparse(dofI[:], dofJ[:], entries[:], ndof, ndof)
    return A
end

function assemble(fnc::Functional, d::Integer)
    qpoints, qweights = quadData(fnc.quadrule, d)

    fes = space(fnc)
    dMap = dofMap(fes, d)
    ndof = Spaces.nDoF(fes)

    nB, nE = size(dMap)
    nP, nD = size(qpoints)

    dofI = zeros(Int, nE, nB)
    
    for ib = 1:nB
        for ie = 1:nE
            dofI[ie,ib] = dMap[ib,ie]
        end
    end
    
    C, V = evaluate(fnc, d)
    Det = abs(evalJacobianDeterminat(fes.mesh, qpoints))
    
    entries = zeros(eltype(qpoints), nE, nB)
    fill_entries!(fnc, entries, C, V, qweights, Det)

    L = sparsevec(dofI[:], entries[:], ndof)
    return L
end

end # of module Operators

