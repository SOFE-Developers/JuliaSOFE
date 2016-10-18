__precompile__()

module Elements

export AbstractElement, LagrangeP1, LagrangeQ1

abstract AbstractElement

#--------------#
# Type Element #
#--------------#
type Element <: AbstractElement
    dim :: Integer
    order :: Integer
    nBasis :: Tuple
    nVertices :: Tuple

    Element(dim::Integer, order::Integer, nB, nV) = new(dim, order, tuple(nB...), tuple(nV...))
end

# Associated methods
# -------------------
Base.size(el::Element) = el.nBasis
Base.size(el::Element, d::Integer) = size(el)[d]

"""

    evalBasis{T<:AbstractFloat}(el::AbstractElement, points::AbstractArray{T,1}, deriv::Integer=0)

Evaluate the element's shape (basis) functions at given `points`
on the reference domain.
"""
function evalBasis{T<:AbstractFloat}(el::AbstractElement,
                                     points::AbstractArray{T,1},
                                     deriv::Integer=0)
    if deriv == 0
        return evalD0Basis(el, points) # nB x nP
    elseif deriv == 1
        return evalD1Basis(el, points) # nB x nP x 1 x nD
    elseif deriv == 2
        return evalD2Basis(el, points) # nB x nP x 1 x nD x nD
    else
        error("Invalid derivation order! ($deriv)")
    end
end

#-----------------#
# Type LagrangeP1 #
#-----------------#
type LagrangeP1 <: AbstractElement
    super :: Element

    LagrangeP1(dim::Integer) = new(Element(dim, 1, collect(1:dim+1), collect(1:dim+1)))
end

# Associated methods
# -------------------
function evalD0Basis{T<:AbstractFloat}(el::LagrangeP1, points::AbstractArray{T,2})
    nP, nD = size(points)
    nB = size(el.super, nD)
    B = zeros(nB, nP)

    B[1,:]     = 1 - sum(points, 2)
    B[2:end,:] = points'

    return B
end

function evalD1Basis{T<:AbstractFloat}(el::LagrangeP1, points::AbstractArray{T,2})
    nP, nD = size(points)
    nB = size(el.super, nD)
    B = zeros(nB, nP, 1, nD)

    B[1,:,:,:] = -1
    for k = 1:nD
        B[k+1,:,:,k] = 1
    end

    return B
end

function evalD2Basis{T<:AbstractFloat}(el::LagrangeP1, points::AbstractArray{T,2})
    nP, nD = size(points)
    nB = size(el.super, nD)
    B = zeros(nB, nP, 1, nD, nD)

    return B
end

#-----------------#
# Type LagrangeQ1 #
#-----------------#
type LagrangeQ1 <: AbstractElement
    super :: Element

    LagrangeQ1(dim::Integer) = new(Element(dim, collect(1:dim+1), collect(1:dim+1), 1))
end


end # of module Elements

