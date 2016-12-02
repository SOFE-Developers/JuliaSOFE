# EXPORTS
export BilinearForm
export trialspace, testspace, coeff, matrix, matrix!, evaluate

#---------------------#
# Bilinear Form Types #
#---------------------#
typealias AbsOp AbstractOperator
typealias AbsCoeff AbstractCoefficient

type BilinearForm{C<:AbsCoeff,U<:AbsOp,V<:AbsOp} <: AbstractVariationalOperator
    trialspace :: FESpace
    testspace :: FESpace
    coeff :: C
    quadrule :: QuadRule
    matrix :: Nullable{SparseMatrixCSC}
end

# Outer Constructors
# -------------------
function BilinearForm{U<:AbsOp,V<:AbsOp,C<:AbsCoeff}(::Type{U}, ::Type{V}, fes::FESpace, coeff::C)
    if issimp(element(fes))
        qrule = QuadRuleSimp2()
    else
        error("Currently only simplical elements supported...")
    end
    mat = Nullable{SparseMatrixCSC}()
    return BilinearForm{C,U,V}(fes, fes, coeff, qrule, mat)
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
trialspace(a::BilinearForm) = getfield(a, :trialspace)
testspace(a::BilinearForm) = getfield(a, :testspace)
coeff(a::BilinearForm) = getfield(a, :coeff)

"""

    matrix(a::BilinearForm)

  Return the discretized version of the bilinear operator `a`.
"""
matrix(a::BilinearForm) = get(getfield(a, :matrix))
matrix!(a::BilinearForm, A::SparseMatrixCSC) = setfield!(a, :matrix, Nullable{SparseMatrixCSC}(A))

function evaluate{C<:AbsCoeff,U<:AbsOp,V<:AbsOp}(a::BilinearForm{C,U,V}, d::Integer)
    points = qpoints(a.quadrule, d)
    
    c = evaluate(coeff(a), points, mesh(trialspace(a)))
    u = evaluate(U, points, trialspace(a))
    v = evaluate(V, points, testspace(a))

    return c, u, v
end


