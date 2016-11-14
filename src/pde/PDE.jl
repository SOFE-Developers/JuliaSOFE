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
function solve(pde::AbstractPDE)
    assemble!(pde)
    compute(pde)
end

function assemble!(pde::AbstractPDE)
    map(assemble!, pde.lhs)
    # for operator in pde.lhs
    #     assemble!(operator)
    # end

    map(assemble!, pde.rhs)
    # for functional in pde.rhs
    #     assemble!(functional)
    # end
end

function compute(pde::AbstractPDE)
    free = space(pde).freeDof
    w = interpolate(space(pde), shift(space(pde)))

    A, b = system(pde)

    u = zeros(free)
    u[!free] = w[!free]
    u[free] = w[free] + A[free,free]\(b-A*w)[free]

    return u
end

function system(pde::AbstractPDE)
    A = matrix(pde.lhs[1])
    for k = 2:length(pde.lhs)
        A += matrix(pde.lhs[k])
    end

    b = vector(pde.rhs[1])
    for k = 2:length(pde.rhs)
        b += vector(pde.rhs[k])
    end

    return A, b
end

end # of module PDE
