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
function Operator{C<:AbsCoeff,U<:AbsOp,V<:AbsOp}(::Type{U}, ::Type{V}, fes::FESpace, coeff::C)
    if issubtype(type_(fes.element), PElement)
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

#-----------#
# Type idid #
#-----------#
type idid <: BilinearOperator
end
idid(fes::FESpace, coeff) = Operator(id, id, fes, coeff)

function evaluate{C<:ConstantCoefficient}(op::Operator{C,id,id}, d::Integer)
    points = qpoints(op.quadrule, d)
    c = value(eltype(points), coeff(op))
    basis = evalBasis(element(space(op)), points, 0)
    return c, basis, basis
end

#---------------#
# Type GradGrad #
#---------------#
type GradGrad <: BilinearOperator
end
GradGrad(fes, coeff) = Operator(Grad, Grad, fes, coeff)

function evaluate{C<:ConstantCoefficient}(op::Operator{C,Grad,Grad}, d::Integer)
    points = qpoints(op.quadrule, d)

    c = value(eltype(points), coeff(op))
    
    dbasis = evalBasis(element(space(op)), points, 1)
    invdphi = evalJacobianInverse(mesh(space(op)), points)

    nE, nP, nW, nD = size(invdphi)
    nB, nP, nW = size(dbasis)

    grad = zeros(eltype(points), nE, nB, nP, nW)

    for jd = 1:nW
        for id = 1:nW
            for ip = 1:nP
                for ib = 1:nB
                    for ie = 1:nE
                        grad[ie,ib,ip,id] += invdphi[ie,ip,jd,id] * dbasis[ib,ip,id]
                    end
                end
            end
        end
    end

    return c, grad, grad
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
