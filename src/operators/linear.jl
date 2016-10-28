abstract LinearOperator <: AbstractOperator

#------------------#
# Functional Types #
#------------------#
type Functional{F<:LinearOperator,C<:AbstractCoefficient} <: LinearOperator
    fes :: FESpace
    coeff :: C
    vector :: Vector
    quadrule :: QuadRule
end

# Outer Constructors
# -------------------
function Functional{F<:LinearOperator,C<:AbstractCoefficient}(::Type{F}, fes::FESpace, coeff::C)
    if issubtype(type_(fes.element), PElement)
        qrule = QuadRuleSimp1()
    else
        error("Currently only simplical elements supported...")
    end

    return Functional{F,C}(fes, coeff, Vector(), qrule)
end

function Functional{F<:LinearOperator, T<:Real}(::Type{F}, fes::FESpace, coeff::T)
    return Functional(F, fes, ScalarCoefficient(coeff))
end

function Functional{F<:LinearOperator, T<:AbstractVector}(::Type{F}, fes::FESpace, coeff::T)
    return Functional(F, fes, VectorCoefficient(coeff))
end
    
function Functional{F<:LinearOperator, T<:AbstractMatrix}(::Type{F}, fes::FESpace, coeff::T)
    return Functional(F, fes, MatrixCoefficient(coeff))
end


# Associated Methods
# -------------------
"""

    matrix(l::LinearOperator)

  Return the discretized version of the linear operator `l`.
"""
@inline vector{T<:LinearOperator}(l::T) = l.vector
@inline getVector{T<:LinearOperator}(l::T) = vector(op)
function setVector!{T<:LinearOperator,S<:Float}(l::T, V::AbstractVector{S})
    l.vector = V;
end

#-----------#
# Type Id   #
#-----------#
type Id <: LinearOperator
end
Id(coeff, fes) = Functional(Id, coeff, fes)

function assemble!(op::Functional{Id}, d::Integer)
    fes = space(op)
    qpoints, qweights = quadData(fes, d)
    C = evalFunction(fes.mesh, coeff(fes), qpoints)
    dMap = dofMap(fes, d)
    ndof = nDoF(fes)
    jdet = evalJacobianDeterminat(fes.mesh, points)
    basis = evalBasis(fes.element, points, 0)

    nB, nE = size(dMap)
    nP, nD = size(qpoints)

    entries = zeros(eltype(qpoints), nE, nB)
    dofI = zeros(Int, nE, nB)

    for ib = 1:nB
        for ie = 1:nE
            dofI[ie,ib] = dMap[ib,ie]
            for ip = 1:nP
                entries[ie,ib] += C[ie,ip] * basis[ib,ip] * qweights[ip] * jdet[ie,ip]
            end
        end
    end

    F = sparsevec(dofI[:], entries[:], ndof)
    setMatrix!(op, F)
end
