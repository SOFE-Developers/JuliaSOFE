__precompile__()

module Operators

using ..Elements
using ..Meshes
using ..Spaces
using ..Quadrature

export AbstractOperator, Operator, IdId, GradGrad, Id
export coeff, space, matrix, assemble!

typealias Float AbstractFloat

# Abstract Operators
abstract AbstractOperator
abstract BilinearOperator <: AbstractOperator
abstract LinearOperator <: AbstractOperator

#-------------------#
# Coefficient Types #
#-------------------#
abstract AbstractCoefficient

typealias CoeffType Union{Real, AbstractVector, AbstractArray, Function}

type Coefficient{T<:CoeffType} <: AbstractCoefficient
    c :: T
end

typealias ScalarCoefficient{T<:Real} Coefficient{T}
typealias VectorCoefficient{T<:AbstractVector} Coefficient{T}
typealias MatrixCoefficient{T<:AbstractMatrix} Coefficient{T}

#----------------#
# Operator Types #
#----------------#
type Operator{O<:BilinearOperator,C<:AbstractCoefficient} <: BilinearOperator
    fes :: FESpace
    coeff :: C
    matrix :: Matrix
    quadrule :: QuadRule
end

# Outer Constructors
# -------------------
function Operator{O<:BilinearOperator,C<:AbstractCoefficient}(::Type{O}, fes::FESpace, coeff::C)
    if issubtype(type_(fes.element), PElement)
        qrule = QuadRuleSimp1()
    else
        error("Currently only simplical elements supported...")
    end

    return Operator{O,C}(fes, coeff, Matrix(), qrule)
end

function Operator{O<:BilinearOperator, T<:Real}(::Type{O}, fes::FESpace, coeff::T)
    return Operator(O, fes, ScalarCoefficient(coeff))
end

function Operator{O<:BilinearOperator, T<:AbstractVector}(::Type{O}, fes::FESpace, coeff::T)
    return Operator(O, fes, VectorCoefficient(coeff))
end
    
function Operator{O<:BilinearOperator, T<:AbstractMatrix}(::Type{O}, fes::FESpace, coeff::T)
    return Operator(O, fes, MatrixCoefficient(coeff))
end

#------------------#
# Functional Types #
#------------------#
type Functional{F<:LinearOperator,C<:AbstractCoefficient} <: LinearOperator
    fes :: FESpace
    coeff :: C
    vector :: Vector
    quadrule :: QuadRule
end

# Outer Constructors
# -------------------
function Functional{F<:LinearOperator,C<:AbstractCoefficient}(::Type{F}, fes::FESpace, coeff::C)
    if issubtype(type_(fes.element), PElement)
        qrule = QuadRuleSimp1()
    else
        error("Currently only simplical elements supported...")
    end

    return Functional{F,C}(fes, coeff, Vector(), qrule)
end

function Functional{F<:LinearOperator, T<:Real}(::Type{F}, fes::FESpace, coeff::T)
    return Functional(F, fes, ScalarCoefficient(coeff))
end

function Functional{F<:LinearOperator, T<:AbstractVector}(::Type{F}, fes::FESpace, coeff::T)
    return Functional(F, fes, VectorCoefficient(coeff))
end
    
function Functional{F<:LinearOperator, T<:AbstractMatrix}(::Type{F}, fes::FESpace, coeff::T)
    return Functional(F, fes, MatrixCoefficient(coeff))
end


# Associated Methods
# -------------------
@inline coeff{T<:AbstractOperator}(op::T) = op.coeff
@inline space{T<:AbstractOperator}(op::T) = op.fes
@inline matrix{T<:BilinearOperator}(a::T) = a.matrix
@inline vector{T<:LinearOperator}(l::T) = l.vector

@inline getCoefficient{T<:AbstractOperator}(op::T) = coeff(op)
@inline getFESpace{T<:AbstractOperator}(op::T) = space(op)
@inline getMatrix{T<:BilinearOperator}(a::T) = matrix(a)
function setMatrix!{T<:BilinearOperator,S<:Float}(a::T, M::AbstractMatrix{S})
    a.matrix = M;
end
@inline getVector{T<:LinearOperator}(l::T) = vector(op)
function setVector!{T<:LinearOperator,S<:Float}(l::T, V::AbstractVector{S})
    l.vector = V;
end

#-----------#
# Type IdId #
#-----------#
type IdId <: BilinearOperator
end
IdId(fes::FESpace, coeff) = Operator(IdId, fes, coeff)

function assemble!(op::Operator{IdId}, d::Integer)
    fes = space(op)
    qpoints, qweights = quadData(op.quadrule, d)
    #C = evalFunction(fes.mesh, coeff(fes), qpoints)
    dMap = dofMap(fes, d)
    ndof = Spaces.nDoF(fes)
    jdet = evalJacobianDeterminat(fes.mesh, qpoints)
    basis = evalBasis(fes.element, qpoints, 0)

    nB, nE = size(dMap)
    nP, nD = size(qpoints)

    C = ones(nE,nP)
    
    entries = zeros(eltype(qpoints), nE, nB, nB)
    dofI = zeros(Int, nE, nB, nB)
    dofJ = zeros(Int, nE, nB, nB)

    for jb = 1:nB
        for ib = 1:nB
            for ie = 1:nE
                dofI[ie,ib,jb] = dMap[ib,ie]
                dofJ[ie,ib,jb] = dMap[jb,ie]
                for ip = 1:nP
                    entries[ie,ib,jb] += C[ie,ip] * basis[jb,ip] * basis[ib,ip] * qweights[ip] * jdet[ie,ip]
                end
            end
        end
    end

    M = sparse(dofI[:], dofJ[:], entries[:], ndof, ndof)
    setMatrix!(op, M)
end

#---------------#
# Type GradGrad #
#---------------#
type GradGrad <: BilinearOperator
end
GradGrad(coeff, fes) = Operator(GradGrad, coeff, fes)

function assemble!(op::Operator{GradGrad})
    fes = space(op)
    qpoints, qweights = quadData(fes, dim(fes.mesh))
    C = evalFunction(fes.mesh, coeff(fes), qpoints)
    dMap = dofMap(fes, dim(fes.mesh))
    ndof = nDoF(fes)
    DPhi = evalReferenceMaps(fes.mesh, qpoints, 1)
    jdet = evalJacobianDeterminat(fes.mesh, qpoints)

    dbasis = evalBasis(fes.element, points, 1)

    nB, nE = size(dMap)
    nP, nD = size(qpoints)

    entries = zeros(eltype(qpoints), nE, nB, nB)
    dofI = zeros(Int, nE, nB, nB)
    dofJ = zeros(Int, nE, nB, nB)

    error("Not Implemented yet...")
    
    for jb = 1:nB
        for ib = 1:nB
            for ie = 1:nE
                dofI[ie,ib,jb] = dMap[ib,ie]
                dofJ[ie,ib,jb] = dMap[jb,ie]
                for ip = 1:nP
                    dphi = DPhi[ie,ip,:,:]
                    grad = dphi'\hcat(dbasis[ib,ip,:], dbasis[jb,ip,:])
                    entries[ie,ib,jb] += C[ie,ip] * basis[jb,ip] * basis[ib,ip] * qweights[ip] * jdet[ie,ip]
                end
            end
        end
    end

    M = sparse(dofI[:], dofJ[:], entries[:], ndof, ndof)
    setMatrix!(op, M)
end

#-----------#
# Type IdId #
#-----------#
type Id <: LinearOperator
end
Id(coeff, fes) = Functional(Id, coeff, fes)

function assemble!(op::Functional{Id}, d::Integer)
    fes = space(op)
    qpoints, qweights = quadData(fes, d)
    C = evalFunction(fes.mesh, coeff(fes), qpoints)
    dMap = dofMap(fes, d)
    ndof = nDoF(fes)
    jdet = evalJacobianDeterminat(fes.mesh, points)
    basis = evalBasis(fes.element, points, 0)

    nB, nE = size(dMap)
    nP, nD = size(qpoints)

    entries = zeros(eltype(qpoints), nE, nB)
    dofI = zeros(Int, nE, nB)

    for ib = 1:nB
        for ie = 1:nE
            dofI[ie,ib] = dMap[ib,ie]
            for ip = 1:nP
                entries[ie,ib] += C[ie,ip] * basis[ib,ip] * qweights[ip] * jdet[ie,ip]
            end
        end
    end

    F = sparsevec(dofI[:], entries[:], ndof)
    setMatrix!(op, F)
end

end # of module Operators

