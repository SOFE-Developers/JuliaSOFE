__precompile__()

module PDE

using ..Operators

typealias Float AbstractFloat

abstract AbstractPDE

#----------#
# Type PDE #
#----------#
type PDE{T<:AbstractPDE} <: AbstractPDE
    lhs :: Array{Operator,1}
    rhs :: Array{Functional,1}
    solution :: Array{Float,1}
end

# Outer Constructors
# -------------------
function PDE{T<:AbstractPDE}(::Type{T}, lhs::AbstractArray{Operator,1}, rhs::AbstractArray{Functional,1})
    return PDE{T}(lhs, rhs, [])
end

# Associated Methods
# -------------------


end # of module PDE
