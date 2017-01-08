__precompile__()

module Problems

using ..Spaces
using ..Operators

import ..Operators: trialspace, testspace, assemble, assemble!

export AbstractPDE, PDE
export solve, assemble!, compute, system

typealias Float AbstractFloat

abstract AbstractPDE

typealias NullableBF Nullable{BilinearForm}
typealias NullableLF Nullable{LinearForm}
typealias None nothing
typealias BF BilinearForm
typealias LF LinearForm
#----------#
# Type PDE #
#----------#
type PDE{T<:AbstractPDE} <: AbstractPDE
    lhs :: Matrix{Vector{BF}}
    rhs :: Vector{Vector{LF}}
    trialspace :: Vector{AbstractFESpace}
    testspace :: Vector{AbstractFESpace}
    solution :: Array{Float,1}
end

# Generic PDE Type
type GenericPDE <: AbstractPDE end
Base.show(io::IO, ::Type{GenericPDE}) = print(io, "GenericPDE")

# Outer Constructors
# -------------------
function PDE{T<:AbstractPDE,B<:BilinearForm,L<:LinearForm}(::Type{T},
                                                           lhs::Matrix{Vector{B}},
                                                           rhs::Vector{Vector{L}})
    neq = length(rhs); @assert size(lhs) == (neq, neq)
    trial = Vector{AbstractFESpace}(neq)
    test = Vector{AbstractFESpace}(neq)
    for i = 1:neq
        for j = 1:neq
            if !isempty(lhs[i,j])
                if !isdefined(trial, j)
                    trial[j] = trialspace(lhs[i,j][1])
                else
                    for k = 1:length(lhs[i,j])
                        @assert is(trial[j], trialspace(lhs[i,j][k]))
                    end
                end
                if !isdefined(test, i)
                    test[i] = testspace(lhs[i,j][1])
                else
                    for k = 1:length(lhs[i,j])
                        @assert is(test[i], testspace(lhs[i,j][k]))
                    end
                    for k = 1:length(rhs[i])
                        @assert is(test[i], testspace(rhs[i][k]))
                    end
                end
            end
        end
    end
    
    return PDE{T}(lhs, rhs, trial, test, [])
end
PDE{B<:BilinearForm, L<:LinearForm}(lhs::Matrix{Vector{B}}, rhs::Vector{Vector{L}}) = PDE(GenericPDE, lhs, rhs)
PDE(lhs::BilinearForm, rhs::LinearForm) = PDE(reshape([[lhs]], 1, 1), [[rhs]])
PDE{T<:AbstractPDE}(::Type{T}, lhs::BilinearForm, rhs::LinearForm) = PDE(T, reshape([[lhs]], 1, 1), [[rhs]])

# Associated Methods
# -------------------
Base.show{T<:AbstractPDE}(io::IO, pde::PDE{T}) = print(io, "PDE{", T, "}")

lhs(pde::AbstractPDE) = getfield(pde, :lhs)
rhs(pde::AbstractPDE) = getfield(pde, :rhs)

trialspace(pde::AbstractPDE, j::Integer) = getfield(pde, :trialspace)[j]
testspace(pde::AbstractPDE, i::Integer) = getfield(pde, :testspace)[i]

function solve(pde::AbstractPDE)
    assemble!(pde)
    compute(pde)
end

function assemble!(pde::AbstractPDE)
    neq = length(rhs(pde))
    for i = 1:neq
        for j = 1:neq
            map(assemble!, lhs(pde)[i,j])
        end
        map(assemble!, rhs(pde)[i])
    end
end

function system(pde::AbstractPDE)
    cumndofs_trial = [0, cumsum([nDoF(fes) for fes in pde.trialspace])...]
    cumndofs_test = [0, cumsum([nDoF(fes) for fes in pde.testspace])...]

    I = Int[]; J = Int[]; V = Float64[]
    b = Float64[]

    neq = length(rhs(pde))
    for i = 1:neq
        for j = 1:neq
            for k = 1:length(lhs(pde)[i,j])
                II, JJ, VV = findnz(matrix(lhs(pde)[i,j][k]))
                I = vcat(I, cumndofs_test[i] + II)
                J = vcat(J, cumndofs_trial[j] + JJ)
                V = vcat(V, VV)
            end
        end
        b = vcat(b, mapreduce(vector, +, rhs(pde)[i]))
    end
    A = sparse(I, J, V)

    return A, b
end

function compute(pde::AbstractPDE)
    free = mapreduce(freeDoF, vcat, pde.trialspace)
    w = vcat([interpolate(fes, shift(fes)) for fes in pde.trialspace]...)
    
    A, b = system(pde)

    u = zeros(Float32, mapreduce(nDoF, +, pde.trialspace))
    u[!free] = w[!free]
    u[free] = w[free] + A[free,free]\(b-A*w)[free]

    return u
end

end # of module Problems
