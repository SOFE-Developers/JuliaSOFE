abstract BilinearOperator <: AbstractOperator

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

# Associated methods
# -------------------
"""

    matrix(a::BilinearOperator)

  Return the discretized version of the bilinear operator `a`.
"""
@inline matrix{T<:BilinearOperator}(a::T) = a.matrix
@inline getMatrix{T<:BilinearOperator}(a::T) = matrix(a)
function setMatrix!{T<:BilinearOperator,S<:Float}(a::T, M::AbstractMatrix{S})
    a.matrix = M;
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

