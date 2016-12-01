__precompile__()

module Problems

using ..Spaces
using ..Operators

import ..Operators: trialspace, testspace, assemble, assemble!

export AbstractPDE, PDE
export solve, assemble!, compute, system

typealias Float AbstractFloat

abstract AbstractPDE

#----------#
# Type PDE #
#----------#
type PDE{T<:AbstractPDE} <: AbstractPDE
    lhs :: Array{BilinearForm,1}
    rhs :: Array{BilinearForm,1}
    solution :: Array{Float,1}
end

type GenericPDE <: AbstractPDE
end

# Outer Constructors
# -------------------
function PDE{T<:AbstractPDE,A<:BilinearForm,L<:LinearForm}(::Type{T},
                                                           lhs::AbstractArray{A,1},
                                                           rhs::AbstractArray{L,1})
    return PDE{T}(lhs, rhs, [])
end
PDE{T<:AbstractPDE}(::Type{T}, lhs::BilinearForm, rhs::LinearForm) = PDE([lhs], [rhs])
PDE{A<:BilinearForm,L<:LinearForm}(lhs::AbstractArray{A,1}, rhs::AbstractArray{L,1}) = PDE(GenericPDE, lhs, rhs)
PDE(lhs::BilinearForm, rhs::LinearForm) = PDE([lhs], [rhs])

# Associated Methods
# -------------------
lhs(pde::AbstractPDE) = getfield(pde, :lhs)
rhs(pde::AbstractPDE) = getfield(pde, :rhs)

trialspace(pde::AbstractPDE) = MixedFESpace(map(trialspace, lhs(pde))...)
testspace(pde::AbstractPDE) = MixedFESpace(map(testspace, lhs(pde))...)

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
    free = freeDoF(space(pde))
    #w = interpolate(space(pde), shift(space(pde)))
    w = vcat([interpolate(fes, shift(fes)) for fes in testspace(pde)]...)
    
    A, b = system(pde)

    u = zeros(Float32, nDoF(trialspace(pde)))
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

end # of module Problems
