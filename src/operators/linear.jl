abstract LinearOperator <: AbstractOperator

#------------------#
# Functional Types #
#------------------#
type Functional{C<:AbstractCoefficient,V<:AbstractOperator} <: LinearOperator
    fes :: FESpace
    coeff :: C
    vector :: Vector
    quadrule :: QuadRule
end

# Outer Constructors
# -------------------
function Functional{C<:AbstractCoefficient,V<:AbstractOperator}(::Type{V}, fes::FESpace, coeff::C)
    if issubtype(type_(fes.element), PElement)
        qrule = QuadRuleSimp2()
    else
        error("Currently only simplical elements supported...")
    end

    return Functional{C,V}(fes, coeff, zeros(nDoF(fes)), qrule)
end

function Functional{V<:AbstractOperator, T<:Real}(::Type{V}, fes::FESpace, coeff::T)
    return Functional(V, fes, ScalarCoefficient(coeff))
end

function Functional{V<:AbstractOperator, T<:AbstractVector}(::Type{V}, fes::FESpace, coeff::T)
    return Functional(V, fes, VectorCoefficient(coeff))
end
    
function Functional{V<:AbstractOperator, T<:AbstractMatrix}(::Type{V}, fes::FESpace, coeff::T)
    return Functional(V, fes, MatrixCoefficient(coeff))
end

function Functional{V<:AbstractOperator, T<:Function}(::Type{V}, fes::FESpace, coeff::T)
    return Functional(V, fes, FunctionCoefficient(coeff))
end

# Associated Methods
# -------------------
"""

    vector(l::LinearOperator)

  Return the discretized version of the linear operator `l`.
"""
@inline vector{T<:LinearOperator}(l::T) = getfield(l, :vector)
@inline getVector{T<:LinearOperator}(l::T) = vector(op)
@inline vector!{T<:LinearOperator}(l::T, L::Vector) = setfield!(l, :vector, L)
setVector!{T<:LinearOperator}(l::T, L::Vector) = vector!(l, L)

#------------#
# Type fid   #
#------------#
type fid <: LinearOperator
end
fid(fes::FESpace, coeff) = Functional(id, fes, coeff)

function evaluate{C<:ConstantCoefficient}(fnc::Functional{C,id}, d::Integer)
    points = qpoints(fnc.quadrule, d)
    c = value(eltype(points), coeff(fnc))
    basis = evalBasis(element(space(fnc)), points, 0)
    return c, basis
end

function evaluate{C<:FunctionCoefficient}(fnc::Functional{C,id}, d::Integer)
    points = qpoints(fnc.quadrule, d)
    c = evaluate(coeff(fnc), points, mesh(space(fnc)))
    basis = evalBasis(element(space(fnc)), points, 0)
    return c, basis
end

