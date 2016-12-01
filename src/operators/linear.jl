# EXPORTS
export LinearForm
export testspace, coeff, vector, vector!, evaluate

#-------------------#
# Linear Form Types #
#-------------------#
type LinearForm{C<:AbstractCoefficient,V<:AbstractOperator} <: AbstractVariationalOperator
    testspace :: FESpace
    coeff :: C
    quadrule :: QuadRule
    vector :: Nullable{Vector}
end

# Outer Constructors
# -------------------
function LinearForm{C<:AbstractCoefficient,V<:AbstractOperator}(::Type{V}, fes::FESpace, coeff::C)
    if issimp(element(fes))
        qrule = QuadRuleSimp2()
    else
        error("Currently only simplical elements supported...")
    end
    vec = Nullable{Vector}()
    return LinearForm{C,V}(fes, coeff, qrule, vec)
end

function LinearForm{V<:AbstractOperator, T<:Real}(::Type{V}, fes::FESpace, coeff::T)
    return LinearForm(V, fes, ScalarCoefficient(coeff))
end

function LinearForm{V<:AbstractOperator, T<:AbstractVector}(::Type{V}, fes::FESpace, coeff::T)
    return LinearForm(V, fes, VectorCoefficient(coeff))
end
    
function LinearForm{V<:AbstractOperator, T<:AbstractMatrix}(::Type{V}, fes::FESpace, coeff::T)
    return LinearForm(V, fes, MatrixCoefficient(coeff))
end

function LinearForm{V<:AbstractOperator, T<:Function}(::Type{V}, fes::FESpace, coeff::T)
    return LinearForm(V, fes, FunctionCoefficient(coeff))
end

# Associated Methods
# -------------------
testspace(l::LinearForm) = getfield(l, :testspace)
coeff(l::LinearForm) = getfield(l, :coeff)

"""

    vector(l::LinearForm)

  Return the discretized version of the linear operator `l`.
"""
vector{T<:LinearForm}(l::T) = getfield(l, :vector)
vector!{T<:LinearForm}(l::T, L::Vector) = setfield!(l, :vector, Nullable{Vector}(L))

function evaluate{C<:AbstractCoefficient,V<:AbstractOperator}(l::LinearForm{C,V}, d::Integer)
    points = qpoints(l.quadrule, d)

    c = evaluate(coeff(l), points, mesh(testspace(l)))
    v = evaluate(V, points, testspace(l))

    return c, v
end
