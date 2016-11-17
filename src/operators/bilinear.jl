# EXPORTS
export BilinearOperator, Operator
export idid
export GradGrad
export matrix, matrix!, evaluate
export getMatrix, setMatrix!

#-------------------------#
# Bilinear Operator Types #
#-------------------------#
abstract BilinearOperator <: AbstractOperator

typealias AbsOp AbstractOperator
typealias AbsCoeff AbstractCoefficient

#----------------#
# Operator Types #
#----------------#
type Operator{C<:AbsCoeff,U<:AbsOp,V<:AbsOp} <: BilinearOperator
    fes :: FESpace
    coeff :: C
    matrix :: SparseMatrixCSC
    quadrule :: QuadRule
end

# Outer Constructors
# -------------------
function Operator{U<:AbsOp,V<:AbsOp,C<:AbsCoeff}(::Type{U}, ::Type{V}, fes::FESpace, coeff::C)
    if issubtype(Elements.type_(fes.element), PElement)
        qrule = QuadRuleSimp2()
    else
        error("Currently only simplical elements supported...")
    end

    return Operator{C,U,V}(fes, coeff, Matrix(), qrule)
end

function Operator{U<:AbsOp,V<:AbsOp,T<:Real}(::Type{U}, ::Type{V}, fes::FESpace, coeff::T)
    return Operator(U, V, fes, ScalarCoefficient(coeff))
end

function Operator{U<:AbsOp,V<:AbsOp,T<:AbstractVector}(::Type{U}, ::Type{V}, fes::FESpace, coeff::T)
    return Operator(U, V, fes, VectorCoefficient(coeff))
end
    
function Operator{U<:AbsOp,V<:AbsOp,T<:AbstractMatrix}(::Type{U}, ::Type{V}, fes::FESpace, coeff::T)
    return Operator(U, V, fes, MatrixCoefficient(coeff))
end

function Operator{U<:AbsOp,V<:AbsOp,T<:Function}(::Type{U}, ::Type{V}, fes::FESpace, coeff::T)
    return Operator(U, V, fes, FunctionCoefficient(coeff))
end

# Associated methods
# -------------------
"""

    matrix(a::BilinearOperator)

  Return the discretized version of the bilinear operator `a`.
"""
@inline matrix{T<:BilinearOperator}(a::T) = getfield(a, :matrix)
@inline getMatrix{T<:BilinearOperator}(a::T) = matrix(a)
@inline matrix!{T<:BilinearOperator}(a::T, A::SparseMatrixCSC) = setfield!(a, :matrix, A)
@inline setMatrix!{T<:BilinearOperator}(a::T, A::SparseMatrixCSC) = matrix!(a, A)

#-------------------------#
# L2 Scalar Product Types #
#-------------------------#
type idid <: BilinearOperator
end
idid(fes::FESpace, coeff) = Operator(id, id, fes, coeff)

function evaluate{C<:ConstantCoefficient}(op::Operator{C,id,id}, d::Integer)
    points = qpoints(op.quadrule, d)
    c = value(eltype(points), coeff(op))
    basis = evalBasis(element(space(op)), points, 0)
    return c, basis[:,:,1], basis[:,:,1]
end

function evaluate{C<:FunctionCoefficient}(op::Operator{C,id,id}, d::Integer)
    points = qpoints(op.quadrule, d)
    c = evaluate(coeff(op), points, mesh(space(op)))
    basis = evalBasis(element(space(op)), points, 0)
    return c, basis[:,:,1], basis[:,:,1]
end

#-----------------#
# Diffusion Types #
#-----------------#
type GradGrad <: BilinearOperator
end
GradGrad(fes, coeff) = Operator(Grad, Grad, fes, coeff)

function evaluate{C<:ConstantCoefficient}(op::Operator{C,Grad,Grad}, d::Integer)
    points = qpoints(op.quadrule, d)

    c = value(eltype(points), coeff(op))
    #c = evaluate(coeff(op), points, mesh(space(op)))
    
    dbasis = evalBasis(element(space(op)), points, 1)
    invdphi = evalJacobianInverse(mesh(space(op)), points)

    nE, nP, nW, nD = size(invdphi)
    nB, nP, nC, nD = size(dbasis)
    @assert nC == 1
    
    grad = zeros(eltype(points), nE, nB, nP, nD)

    for id = 1:nD
        for iw = 1:nW
            for ip = 1:nP
                for ib = 1:nB
                    for ie = 1:nE
                        grad[ie,ib,ip,id] += invdphi[ie,ip,iw,id] * dbasis[ib,ip,1,id]
                    end
                end
            end
        end
    end

    return c, grad, grad
end

function evaluate{C<:FunctionCoefficient}(op::Operator{C,Grad,Grad}, d::Integer)
    points = qpoints(op.quadrule, d)

    c = evaluate(coeff(op), points, mesh(space(op)))
    
    dbasis = evalBasis(element(space(op)), points, 1)
    invdphi = evalJacobianInverse(mesh(space(op)), points)

    nE, nP, nW, nD = size(invdphi)
    nB, nP, nC, nD = size(dbasis)
    @assert nC == 1

    grad = zeros(eltype(points), nE, nB, nP, nD)

    for id = 1:nD
        for iw = 1:nW
            for ip = 1:nP
                for ib = 1:nB
                    for ie = 1:nE
                        grad[ie,ib,ip,id] += invdphi[ie,ip,iw,id] * dbasis[ib,ip,1,id]
                    end
                end
            end
        end
    end

    return c, grad, grad
end

