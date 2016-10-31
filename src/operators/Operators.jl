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

abstract ScalarOperator <: AbstractOperator
abstract VectorOperator <: AbstractOperator
abstract MatrixOperator <: AbstractOperator

typealias op ScalarOperator
typealias Op VectorOperator
typealias OP MatrixOperator

typealias data ScalarCoefficient
typealias Data VectorCoefficient
typealias DATA MatrixCoefficient

#-------------------------#
# Identity Operator Types #
#-------------------------#
type id <: ScalarOperator end
type Id <: VectorOperator end
type ID <: MatrixOperator end

#-------------------------#
# Gradient Operator Types #
#-------------------------#
type grad <: ScalarOperator end
type Grad <: VectorOperator end
type GRAD <: MatrixOperator end

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
#include("linear.jl")

function assemble(op::Operator, d::Integer)
    fes = space(op)
    qpoints, qweights = quadData(op.quadrule, d)
    
    #C = evaluate(coeff(op), fes.mesh, qpoints)
    C = convert(eltype(qpoints), value(coeff(op)))
    U, V = evaluate(op, d)
    Det = evalJacobianDeterminat(fes.mesh, qpoints)
    Det = abs(Det)

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
    
    entries = zeros(eltype(qpoints), nE, nB, nB)
    fill_entries!(op, entries, C, U, V, qweights, Det)

    M = sparse(dofI[:], dofJ[:], entries[:], ndof, ndof)
    return M
end

# Constant Coefficients
function fill_entries!{T<:Real,Tc<:data,Tu<:op,Tv<:op}(::Operator{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3},
                                                       C::T,
                                                       U::AbstractArray{T,2},
                                                       V::AbstractArray{T,2},
                                                       w::AbstractArray{T,1},
                                                       D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nP = size(w, 1)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    E[ie,ib,jb] += C * U[jb,ip] * V[ib,ip] * w[ip] * D[ie,ip]
                end
            end
        end
    end
    return nothing
end

function fill_entries!{T<:Real,Tc<:data,Tu<:Op,Tv<:Op}(::Operator{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3},
                                                       C::T,
                                                       U::AbstractArray{T,4},
                                                       V::AbstractArray{T,4},
                                                       w::AbstractArray{T,1},
                                                       D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for id = 1:nW
                        E[ie,ib,jb] += C * U[ie,jb,ip,id] * V[ie,ib,ip,id] * w[ip] * D[ie,ip]
                    end
                end
            end
        end
    end
    return nothing
end

function fill_entries!{T<:Real,Tc<:Data,Tu<:Op,Tv<:op}(::Operator{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3},
                                                       C::AbstractArray{T,1},
                                                       U::AbstractArray{T,4},
                                                       V::AbstractArray{T,2},
                                                       w::AbstractArray{T,1},
                                                       D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for id = 1:nW
                        E[ie,ib,jb] += C[id] * U[ie,jb,ip,id] * V[ib,ip] * w[ip] * D[ie,ip]
                    end
                end
            end
        end
    end
    return nothing
end

function fill_entries!{T<:Real,Tc<:DATA,Tu<:Op,Tv<:Op}(::Operator{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3},
                                                       C::AbstractArray{T,2},
                                                       U::AbstractArray{T,4},
                                                       V::AbstractArray{T,4},
                                                       w::AbstractArray{T,1},
                                                       D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for jd = 1:nW
                        for id = 1:nW
                            E[ie,ib,jb] += C[id,jd] * U[ie,jb,ip,jd] * V[ie,ib,ip,id] * w[ip] * D[ie,ip]
                        end
                    end
                end
            end
        end
    end
    return nothing
end

# Non-Constant Coefficients
function fill_data_op_op!{T<:Real}(E::AbstractArray{T,3},
                                   C::AbstractArray{T,2}, U::AbstractArray{T,2}, V::AbstractArray{T,2},
                                   w::AbstractArray{T,1}, D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    E[ie,ib,jb] += C[ie,ip] * U[jb,ip] * V[ib,ip] * w[ip] * D[ie,ip]
                end
            end
        end
    end
    return nothing
end

function fill_data_Op_Op!{T<:Real}(E::AbstractArray{T,3},
                                   C::AbstractArray{T,2}, U::AbstractArray{T,4}, V::AbstractArray{T,4},
                                   w::AbstractArray{T,1}, D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for id = 1:nW
                        E[ie,ib,jb] += C[ie,ip] * U[ie,jb,ip,id] * V[ie,ib,ip,id] * w[ip] * D[ie,ip]
                    end
                end
            end
        end
    end
    return nothing
end

function fill_Data_Op_op!{T<:Real}(E::AbstractArray{T,3},
                                   C::AbstractArray{T,3}, U::AbstractArray{T,4}, V::AbstractArray{T,2},
                                   w::AbstractArray{T,1}, D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for id = 1:nW
                        E[ie,ib,jb] += C[ie,ip,id] * U[ie,jb,ip,id] * V[ib,ip] * w[ip] * D[ie,ip]
                    end
                end
            end
        end
    end
    return nothing
end

function fill_DATA_Op_Op!{T<:Real}(E::AbstractArray{T,3},
                                   C::AbstractArray{T,4}, U::AbstractArray{T,4}, V::AbstractArray{T,4},
                                   w::AbstractArray{T,1}, D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for jd = 1:nW
                        for id = 1:nW
                            E[ie,ib,jb] += C[ie,ip,id,jd] * U[ie,jb,ip,jd] * V[ie,ib,ip,id] * w[ip] * D[ie,ip]
                        end
                    end
                end
            end
        end
    end
    return nothing
end
end # of module Operators

