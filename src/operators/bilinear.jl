# EXPORTS
export BilinearOperator, Operator
export idid
export GradGrad
export matrix, matrix!, evaluate

#-------------------------#
# Bilinear Operator Types #
#-------------------------#
typealias AbsOp AbstractOperator
typealias AbsCoeff AbstractCoefficient

type BilinearForm{C<:AbsCoeff,U<:AbsOp,V<:AbsOp} <: AbstractVariationalOperator
    trialspace :: FESpace
    testspace :: FESpace
    coeff :: C
    matrix :: Nullable{SparseMatrixCSC}
    quadrule :: QuadRule
end

# Outer Constructors
# -------------------
function BilinearForm{U<:AbsOp,V<:AbsOp,C<:AbsCoeff}(::Type{U}, ::Type{V}, fes::FESpace, coeff::C)
    if issimp(element(fes))
        qrule = QuadRuleSimp2()
    else
        error("Currently only simplical elements supported...")
    end

    return BilinearForm{C,U,V}(fes, coeff, Matrix(), qrule)
end

function BilinearForm{U<:AbsOp,V<:AbsOp,T<:Real}(::Type{U}, ::Type{V}, fes::FESpace, coeff::T)
    return BilinearForm(U, V, fes, ScalarCoefficient(coeff))
end

function BilinearForm{U<:AbsOp,V<:AbsOp,T<:AbstractVector}(::Type{U}, ::Type{V}, fes::FESpace, coeff::T)
    return BilinearForm(U, V, fes, VectorCoefficient(coeff))
end
    
function BilinearForm{U<:AbsOp,V<:AbsOp,T<:AbstractMatrix}(::Type{U}, ::Type{V}, fes::FESpace, coeff::T)
    return BilinearForm(U, V, fes, MatrixCoefficient(coeff))
end

function BilinearForm{U<:AbsOp,V<:AbsOp,T<:Function}(::Type{U}, ::Type{V}, fes::FESpace, coeff::T)
    return BilinearForm(U, V, fes, FunctionCoefficient(coeff))
end

# Associated methods
# -------------------
trialspace(op::BilinearForm) = getfield(op, :trialspace)
testspace(op::BilinearForm) = getfield(op, :testspace)

"""

    matrix(a::BilinearForm)

  Return the discretized version of the bilinear operator `a`.
"""
@inline matrix(a::BilinearForm) = getfield(a, :matrix)
@inline matrix!(a::BilinearForm, A::SparseMatrixCSC) = setfield!(a, :matrix, A)

function evaluate(a::BilinearForm, d::Integer)
    points = qpoints(op.quadrule, d)

    C = evaluate(coeff(op), points, mesh(fespace(trialop(a))))
end

#-------------------------#
# L2 Scalar Product Types #
#-------------------------#
type idid <: BilinearForm
end
idid(fes::FESpace, coeff) = Operator(id, id, fes, coeff)

function evaluate{C<:AbstractCoefficient}(op::Operator{C,id,id}, d::Integer)
    points = qpoints(op.quadrule, d)
    c = evaluate(coeff(op), points, mesh(space(op)))
    basis = evalBasis(element(space(op)), points, 0)
    @assert ndims(basis) == 3 && size(basis, 3) == 1
    basis = view(basis, :, :, 1)
    return c, basis, basis
end

#-----------------#
# Diffusion Types #
#-----------------#
type GradGrad <: BilinearForm
end
GradGrad(fes, coeff) = Operator(Grad, Grad, fes, coeff)

function evaluate{C<:AbstractCoefficient}(op::Operator{C,Grad,Grad}, d::Integer)
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

