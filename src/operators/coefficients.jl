#-------------------#
# Coefficient Types #
#-------------------#
abstract AbstractCoefficient

#--------------------------#
# Type ConstantCoefficient #
#--------------------------#
typealias ConstType Union{Real, AbstractVector, AbstractArray}
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

function (c::ConstantCoefficient)(x)
    return value(c)
end

function evaluate{T<:Float,C<:ConstantCoefficient}(m::Mesh, c::C, points::AbstractArray{T,2})
    P = evalReferenceMaps(m, points)
    nE, nP, nW = size(P)

    val = value(c)
    R = zeros(T, nE, nP, size(val)...)
    for ip = 1:nP
        for ie = 1:nE
            R[ie,ip,:] = val
        end
    end
    return R
end

#--------------------------#
# Type FunctionCoefficient #
#--------------------------#
type FunctionCoefficient{T<:Function} <: AbstractCoefficient
    func :: T
end

# Associated Methods
# -------------------
function (f::FunctionCoefficient)(x)
    return getfield(f, :func)(x)
end

function evaluate{T<:Float,F<:FunctionCoefficient}(m::Mesh, f::F, points::AbstractArray{T,2})
    P = evalReferenceMaps(m, points)
    nE, nP, nW = size(P)

    R = f(reshape(P, (nE*nP,nW)))
    return reshape(R, nE, nP, size(R,(2:ndims(R))...)...)
end

