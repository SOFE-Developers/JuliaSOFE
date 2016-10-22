__precompile__()

module Operators

using ..Elements
using ..Meshes
using ..FESpace

typealias Float AbstractFloat

abstract AbstractOperator
abstract BilinearOperator <: AbstractOperator
abstract LinearOperator <: AbstractOperator

#---------------#
# Type Operator #
#---------------#
type Operator{O<:AbstractOperator, T<:Float} <: AbstractOperator
    coeff :: Function
    fes :: FESpace
    matrix :: AbstractArray{T,2}
end
function Operator{O<:AbstractOperator}(::Type{O}, coeff::Function, fes::FESpace)
    return Operator{O}(coeff, fes)
end
function Operator{O<:AbstractOperator, T<:Real}(::Type{O}, coeff::T, fes::FESpace)
    return Operator{O}(x->convert(coeff, typeof(x))) * ones(typeof(x), size(x,1), fes)
end

# Associated Methods
# -------------------
@inline coeff(op::Operator) = op.coeff
@inline space(op::Operator) = op.fes
@inline matrix(op::Operator) = op.matrix
@inline getCoefficient(op::Operator) = coeff(op)
@inline getFESpace(op::Operator) = space(op)
@inline getMatrix(op::Operator) = matrix(op)

#-----------#
# Type IdId #
#-----------#
type IdId <: BilinearOperator
end
IdId(coeff, fes) = Operator(IdId, coeff, fes)

function assemble(op::Element{IdId}, d::Integer)
    fes = space(op)
    qpoints, qweights = quadData(fes, d)
    C = evalFunction(fes.mesh, coeff(fes), qpoints)
    dMap = dofMap(fes, d)
    ndof = nDoF(fes)
    jdet = evalJacobianDeterminat(fes.mesh, points)
    basis = evalBasis(fes.element, points, 0)

    nB, nE = size(dMap)
    nP, nD = size(qpoints)

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

    return sparse(dofI[:], dofJ[:], entries[:], ndof, ndof)
end

#---------------#
# Type GradGrad #
#---------------#
type GradGrad <: BilinearOperator
end
GradGrad(coeff, fes) = Operator(GradGrad, coeff, fes)

function assemble(op::Element{GradGrad})
    fes = space(op)
    qpoints, qweights = quadData(fes, dim(fes.mesh))
    C = evalFunction(fes.mesh, coeff(fes), qpoints)
    dMap = dofMap(fes, dim(fes.mesh))
    ndof = nDoF(fes)
    jdet = evalJacobianDeterminat(fes.mesh, points)

    error("Not implemented yet!") # -> glob derivatives
    
    dbasis = evalBasis(fes.element, points, 1)

    nB, nE = size(dMap)
    nP, nD = size(qpoints)

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

    return sparse(dofI[:], dofJ[:], entries[:], ndof, ndof)
end

#-----------#
# Type IdId #
#-----------#
type Id <: LinearOperator
end
Id(coeff, fes) = Operator(Id, coeff, fes)

function assemble(op::Element{Id}, d::Integer)
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

    return sparsevec(dofI[:], entries[:], ndof)
end

end # of module Operators

