abstract BilinearOperator <: AbstractOperator

typealias AbsOp AbstractOperator
typealias AbsCoeff AbstractCoefficient

#----------------#
# Operator Types #
#----------------#
type Operator{C<:AbsCoeff,U<:AbsOp,V<:AbsOp} <: BilinearOperator
    fes :: FESpace
    coeff :: C
    matrix :: Matrix
    quadrule :: QuadRule
end

# Outer Constructors
# -------------------
function Operator{C<:AbsCoeff,U<:AbsOp,V<:AbsOp}(::Type{U}, ::Type{V}, fes::FESpace, coeff::C)
    if issubtype(type_(fes.element), PElement)
        qrule = QuadRuleSimp1()
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
# Type idid #
#-----------#
type idid <: BilinearOperator
end
idid(fes::FESpace, coeff) = Operator(id, id, fes, coeff)

function evaluate{C<:AbsCoeff}(op::Operator{C,id,id}, d::Integer)
    points = qpoints(op.quadrule, d)
    fes = space(op)
    basis = evalBasis(fes.element, points, 0)
    return basis, basis
end

#---------------#
# Type GradGrad #
#---------------#
type GradGrad <: BilinearOperator
end
GradGrad(fes, coeff) = Operator(Grad, Grad, fes, coeff)

function evaluate{C<:AbsCoeff}(op::Operator{C,Grad,Grad}, d::Integer)
    points = qpoints(op.quadrule, d)
    fes = space(op)
    dbasis = evalBasis(fes.element, points, 1)
    invdphi = evalJacobianInverse(fes.mesh, points)

    nE, nP, nW, nW = size(invdphi)
    nB, nP, nW = size(dbasis)

    grad = zeros(eltype(points), nE, nB, nP, nW)

    for jd = 1:nW
        for id = 1:nW
            for ib = 1:nB
                for ip = 1:nP
                    for ie = 1:nE
                        grad[ie,ip,ib,id] += invdphi[ie,ip,jd,id] * dbasis[ib,ip,id]
                    end
                end
            end
        end
    end
    return grad, grad
end

