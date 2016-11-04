# EXPORTS
export AbstractCoefficient
export ConstantCoefficient, ScalarCoefficient, VectorCoefficient, MatrixCoefficient,
export FunctionCoefficient
export value, func, evaluate

#-------------------#
# Coefficient Types #
#-------------------#
"""
Abstract base type for all operator coefficients. 
"""
abstract AbstractCoefficient

# Associated Methods
# -------------------
"""
    evaluate{T<:Float}(coeff::AbstractCoefficient, points::AbstractArray{T})

  Evaluate the coefficient `coeff` in given (global) `points`.
"""
function evaluate{T<:Float}(coeff::AbstractCoefficient, points::AbstractArray{T,2})
end

"""
    evaluate{T<:Float}(coeff::AbstractCoefficient, points::AbstractArray{T}, m::Mesh)

  Evaluate the coefficient `coeff` in global points of the mesh `m` 
  resulting from the evaluation of its referefence maps in given local `points`.
"""
function evaluate{T<:Float}(coeff::AbstractCoefficient, points::AbstractArray{T,2}, m::Mesh)
end

#--------------------------#
# Type ConstantCoefficient #
#--------------------------#
typealias ConstType Union{Real, AbstractVector, AbstractArray}

"""
Type for constant coefficients such as scalars, vectors or matrices.
"""
type ConstantCoefficient{T<:ConstType} <: AbstractCoefficient
    value :: T
end

typealias ScalarCoefficient{T<:Real} ConstantCoefficient{T}
typealias VectorCoefficient{T<:Real, S<:AbstractVector{T}} ConstantCoefficient{S}
typealias MatrixCoefficient{T<:Real, S<:AbstractMatrix{T}} ConstantCoefficient{S}

# Associated methods
# -------------------
"""

    value(c::ConstantCoefficient)

  Return the constant value of the coefficient `c`.
"""
@inline value(c::ConstantCoefficient) = getfield(c, :value)

"""

    value{T<:Real}(::Type{T), c::ConstantCoefficient)

  Return the constant value of the coefficient `c` where its 
  element type is promoted to type `T` if necessary.
"""
value{T<:Real,S<:Real}(::Type{T}, c::ConstantCoefficient{S}) =
    convert(promote_type(T, S), value(c))
value{T<:Real,S<:AbstractVector}(::Type{T}, c::ConstantCoefficient{S}) =
    convert(Vector{promote_type(T,eltype(S))}, value(c))
value{T<:Real,S<:AbstractMatrix}(::Type{T}, c::ConstantCoefficient{S}) =
    convert(Matrix{promote_type(T,eltype(S))}, value(c))

function (c::ConstantCoefficient)(x)
    return value(eltype(x), c)
end

function evaluate{T<:Float,C<:ConstantCoefficient}(c::C, points::AbstractArray{T,2}, m::Mesh)
    P = evalReferenceMaps(m, points)
    nE, nP, nW = size(P)

    val = value(T, c)
    R = zeros(T, nE, nP, size(val)...)
    for ip = 1:nP
        for ie = 1:nE
            R[ie,ip,:] = val
        end
    end
    return R
end

function evaluate{T<:Float,C<:ConstantCoefficient}(c::C, points::AbstractArray{T,2})
    nP, nW = size(points)
    val = value(T, c)
    R = zeros(T, nP, size(val)...)
    for ip = 1:nP
        R[ip,:] = val
    end
    return R
end

#--------------------------#
# Type FunctionCoefficient #
#--------------------------#
"""
Type for non-constant callable coefficients.
"""
type FunctionCoefficient{T<:Function} <: AbstractCoefficient
    func :: T
end

# Associated Methods
# -------------------
"""

    func(f::FunctionCoefficient)

  Return the callable function instance of the coefficient `f`.
"""
@inline func(f::FunctionCoefficient) = getfield(f, :func)
    
function (f::FunctionCoefficient){T<:Real}(xs::T...)
    return func(fc)(xs...)
end
(f::FunctionCoefficient){T<:Real}(x::AbstractVector{T}) = f(x...)
function (f::FunctionCoefficient){T<:Real}(X::AbstractArray{T,2})
    Y = [view(X, i, :) for i = 1:size(X,1)]
    return f.(Y)
end

function evaluate{F<:FunctionCoefficient,T<:Float}(f::F, points::AbstractArray{T,2})
    return f(points)
end

function evaluate{F<:FunctionCoefficient,T<:Float}(f::F, points::AbstractArray{T,2}, m::Mesh)
    P = evalReferenceMaps(m, points)
    nE, nP, nW = size(P)

    R = f(reshape(P, (nE*nP,nW)))
    nF = size(R, (2,(3:ndims(R))...)...)
    R = (nF == 1) ? reshape(R, nE, nP) : reshape(R, nE, nP, nF...)

    #return R
    return ndarray(R)
end

