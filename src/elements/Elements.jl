__precompile__()

module Elements

import ..Helpers: dimension

export AbstractElement, PElement, QElement, Element
export LagrangeP1, LagrangeQ1
export dimension, order, nBasis, nVertices, dofTuple, nDoF, evalBasis, isnodal
export issimplical, isorthotopic
export ndims

abstract AbstractElement
abstract PElement <: AbstractElement
abstract QElement <: AbstractElement

typealias ElementTypes Union{PElement, QElement}

typealias Float AbstractFloat

#--------------#
# Type Element #
#--------------#
# type Element{T<:ElementTypes} <: AbstractElement
#     dimension :: Int
# end
type Element{T<:ElementTypes,N} <: AbstractElement
    dimension :: Int
end
# type FiniteElement{T<:ElementTypes,N} <: AbstractElement
#     dimension :: Int
# end

# typealias Element{T} FiniteElement{T,1}
# typealias VectorElement{T,N} FiniteElement{T,N}

# Outer Constructors
# -------------------
Element{T<:ElementTypes}(::Type{T}, dim::Integer) = Element{T,1}(dim)
VectorElement{T<:ElementTypes}(::Type{T}, dim::Integer, cmp::Integer) = Element{T,Int(cmp)}(dim)

# Associated methods
# -------------------
"""

    dimension(el::Element)

  The topological dimension of the element.
"""
@inline dimension(el::Element) = el.dimension

@inline Base.ndims{T,N}(::Element{T,N}) = isa(N, Integer) ? N : throw(TypeError)

@inline issimplical{T<:PElement}(::Element{T}) = true
@inline issimplical{T<:QElement}(::Element{T}) = false
@inline isorthotopic{T<:PElement}(::Element{T}) = false
@inline isorthotopic{T<:QElement}(::Element{T}) = true


"""

    order(el::Element)

  The polynomial order of the element's basis (shape) functions.
"""
function order(el::Element)
end

"""

    nBasis(el::Element, [d::Integer])

  The number of basis (shape) functions associated with the 
  `d`-dimensional entities of the element.
"""
@inline nBasis(el::Element, d::Integer) = nBasis(el)[d]

"""

    nVertices(el::Element, [d::Integer])

  The number of vertices that define the `d`-dimensional
  entities of the element.
"""
@inline nVertices(el::Element, d::Integer) = nVertices(el)[d]
@inline nVertices{T<:PElement}(el::Element{T}) = tuple(2:(dimension(el)+1)...)
@inline nVertices{T<:QElement}(el::Element{T}) = tuple(2.^(1:dimension(el))...)

@inline nDoF(el::Element) = nDoF(el, dimension(el))
function nDoF{T<:PElement}(el::Element{T}, d::Integer)
    #return map(*, dofTuple(el), binomial(dimension(el)+1, k) for k = 1:dimension(el)+1)
    return map(*, dofTuple(el), binomial(d+1, k) for k = 1:d+1)
end
function nDoF{T<:QElement}(el::Element{T})
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
    nC = ndims(el)

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

include("lagrange.jl")

include("mixed.jl")

end # of module Elements

