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

    return Functional{C,V}(fes, coeff, zeros(Spaces.nDoF(fes)), qrule)
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

    vector(l::LinearOperator)

  Return the discretized version of the linear operator `l`.
"""
@inline vector{T<:LinearOperator}(l::T) = getfield(l, :vector)
@inline getVector{T<:LinearOperator}(l::T) = vector(op)
function setVector!{T<:LinearOperator,S<:Float}(l::T, V::AbstractVector{S})
    l.vector = V;
end

#------------#
# Type fid   #
#------------#
type fid <: LinearOperator
end
fid(fes::FESpace, coeff) = Functional(id, fes, coeff)

function evaluate{C<:AbstractCoefficient}(fnc::Functional{C,id}, d::Integer)
    points = qpoints(fnc.quadrule, d)
    c = value(eltype(points), coeff(fnc))
    basis = evalBasis(element(space(fnc)), points, 0)
    return c, basis
end

#---------------------#
# Assembling Routines #
#---------------------#

# Constant Coefficients
function fill_entries!{T<:Real,Tc<:data,Tv<:op}(::Functional{Tc,Tv},
                                                E::AbstractArray{T,2},
                                                C::T,
                                                V::AbstractArray{T,2},
                                                w::AbstractArray{T,1},
                                                D::AbstractArray{T,2})
    nE, nB = size(E)
    nP = size(w, 1)
    for ib = 1:nB
        for ie = 1:nE
            for ip = 1:nP
                E[ie,ib] += C * V[ib,ip] * w[ip] * D[ie,ip]
            end
        end
    end

    return nothing
end

