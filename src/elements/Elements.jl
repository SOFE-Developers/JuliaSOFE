__precompile__()

module Elements

import CMesh: dimension

export AbstractElement, PElement, QElement, Element
export dimension, order, nBasis, nVertices, dofTuple, nDoF, evalBasis, isnodal
export issimp, isorth

abstract AbstractElement{T}

abstract PElement <: AbstractElement
abstract QElement <: AbstractElement

typealias ElementTypes Union{PElement, QElement}

typealias Float AbstractFloat

#--------------#
# Type Element #
#--------------#
type Element{T<:ElementTypes} <: AbstractElement{T}
    dimension :: Int
end

# Outer Constructors
# -------------------
Element{T<:ElementTypes}(::Type{T}, dim::Integer) = Element{T}(dim)

# Associated methods
# -------------------
"""

    dimension(el::Element)

  The topological dimension of the element.
"""
dimension(el::Element) = getfield(el, :dimension)

issimp{T<:PElement}(::AbstractElement{T}) = true
issimp{T<:QElement}(::AbstractElement{T}) = false
isorth{T<:PElement}(::AbstractElement{T}) = false
isorth{T<:QElement}(::AbstractElement{T}) = true

"""

    order(el::AbstractElement)

  The polynomial order of the element's basis (shape) functions.
"""
function order(el::AbstractElement)
end

"""

    nBasis(el::Element, [d::Integer])

  The number of basis (shape) functions associated with the 
  `d`-dimensional entities of the element.
"""
@inline nBasis(el::Element, d::Integer) = nBasis(el)[d]

"""

    nVertices(el::Element [, d::Integer])

  The number of vertices that define the `d`-dimensional
  entities of the element.
"""
nVertices(el::AbstractElement, d::Integer) = nVertices(el)[d]
nVertices{T<:PElement}(el::AbstractElement{T}) = tuple(2:(dimension(el)+1)...)
nVertices{T<:QElement}(el::AbstractElement{T}) = tuple(2.^(1:dimension(el))...)

nDoF(el::AbstractElement) = nDoF(el, dimension(el))
function nDoF{T<:PElement}(el::AbstractElement{T}, d::Integer)
    #return map(*, dofTuple(el), binomial(dimension(el)+1, k) for k = 1:dimension(el)+1)
    return map(*, dofTuple(el), binomial(d+1, k) for k = 1:d+1)
end
function nDoF{T<:QElement}(el::AbstractElement{T})
    if dimension(el) == 1
        return map(*, dofTuple(el), (2,1))
    elseif dimension(el) == 2
        return map(*, dofTuple(el), (4,4,1))
    elseif dimension(el) == 3
        return map(*, dofTuple(el), (8,12,6,1))
    else
        error("Invalid element dimension ", dimension(el))
    end
end

"""

    evalBasis{T<:AbstractFloat}(el::Element, points::AbstractArray{T,1}, deriv::Integer=0)

  Evaluate the element's shape (basis) functions at given `points`
  on the reference domain.
"""
function evalBasis{T<:AbstractFloat}(el::Element,
                                     points::AbstractArray{T,2},
                                     deriv::Integer=0)
    nP, nW = size(points)
    nB = nBasis(el, nW)
    nC = 1 #ndims(el)

    if deriv == 0
        B = zeros(T, nB, nP, nC)
        return evalD0Basis!(el, points, B) # nB x nP x nC
    elseif deriv == 1
        B = zeros(T, nB, nP, nC, nW)
        return evalD1Basis!(el, points, B) # nB x nP x nC x nW
    elseif deriv == 2
        B = zeros(T, nB, nP, nC, nW, nW)
        return evalD2Basis!(el, points, B) # nB x nP x nC x nW x nW
    else
        error("Invalid derivation order! ($deriv)")
    end
end

include("vector.jl")
include("mixed.jl")

include("lagrange.jl")

include("hierarchic.jl")

end # of module Elements

