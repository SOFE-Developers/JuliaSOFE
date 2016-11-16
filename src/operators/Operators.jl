__precompile__()

module Operators

using ..Elements
using ..Meshes
using ..Spaces
using ..Quadrature
using ..Helpers

export AbstractOperator
export space, coeff, assemble, assemble!

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
space{T<:AbstractOperator}(op::T) = getfield(op, :fes)

"""

    coeff(op::AbstractOperator)

  Return the coefficient of the operator `op`.
"""
coeff{T<:AbstractOperator}(op::T) = getfield(op, :coeff)

# Includes
include("bilinear.jl")
include("linear.jl")

#---------------------#
# Assembling Routines #
#---------------------#
assemble(op::AbstractOperator) = assemble(op, dimension(mesh(space(op))))
assemble!(op::AbstractOperator) = assemble!(op, dimension(mesh(space(op))))
assemble!(op::Operator, d::Integer) = matrix!(op, assemble(op, d))
assemble!(fnc::Functional, d::Integer) = vector!(fnc, assemble(fnc, d))

function assemble(op::Operator, d::Integer)
    qpoints, qweights = quadData(op.quadrule, d)
    
    fes = space(op)
    dMap = dofMap(fes, d=d)
    ndof = nDoF(fes)
    
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
    dMap = dofMap(fes, d=d)
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
    return full(L)
end

include("entries.jl")

end # of module Operators

