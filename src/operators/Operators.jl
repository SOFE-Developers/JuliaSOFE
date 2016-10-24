__precompile__()

module Operators

using ..Elements
using ..Meshes
using ..Spaces
using ..Quadrature

export AbstractOperator, Operator, IdId, GradGrad, Id
export coeff, space, matrix, assemble!

typealias Float AbstractFloat

abstract AbstractOperator
abstract BilinearOperator <: AbstractOperator
abstract LinearOperator <: AbstractOperator

#---------------#
# Type Operator #
#---------------#
type Operator{O<:AbstractOperator} <: AbstractOperator
    fes :: FESpace
    coeff :: Function
    matrix :: AbstractArray
    quadrule :: QuadRule
end
function Operator{O<:AbstractOperator}(::Type{O}, fes::FESpace, coeff::Function)
    if issubtype(type_(fes.element), PElement)
        qrule = QuadRuleSimp1()
    else
        error("Currently only simplical elements supported...")
    end

    return Operator{O}(fes, coeff, [], qrule)
end
function Operator{O<:AbstractOperator, T<:Real}(::Type{O}, fes::FESpace, coeff::T)
    return Operator(O, fes,
                    x->convert(coeff, typeof(x)) * ones(typeof(x), size(x,1)))
end

# Associated Methods
# -------------------
@inline coeff(op::Operator) = op.coeff
@inline space(op::Operator) = op.fes
@inline matrix(op::Operator) = op.matrix
@inline getCoefficient(op::Operator) = coeff(op)
@inline getFESpace(op::Operator) = space(op)
@inline getMatrix(op::Operator) = matrix(op)
function setMatrix!{T<:Float}(op::Operator, M::AbstractArray{T,2})
    op.matrix = M;
end

#-----------#
# Type IdId #
#-----------#
type IdId <: BilinearOperator
    IdId(fes::FESpace, coeff) = Operator(IdId, fes, coeff)
end

function assemble!(op::Operator{IdId}, d::Integer)
    fes = space(op)
    qpoints, qweights = quadData(op.quadrule, d)
    #C = evalFunction(fes.mesh, coeff(fes), qpoints)
    dMap = dofMap(fes, d)
    ndof = Spaces.nDoF(fes)
    jdet = evalJacobianDeterminat(fes.mesh, points)
    basis = evalBasis(fes.element, points, 0)

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
Id(coeff, fes) = Operator(Id, coeff, fes)

function assemble!(op::Operator{Id}, d::Integer)
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
