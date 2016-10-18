__precompile__()

module Elements

export AbstractElement, Element, LagrangeP1

abstract AbstractElement

#--------------#
# Type Element #
#--------------#
type Element{E<:AbstractElement}
    type_ :: E
    dimension :: Integer
end

function Element{E<:AbstractElement}(::Type{E}, dimension)
    return Element{E}(E(), dimension)
end

# Associated Methods
# -------------------
@inline dimension(el::Element) = el.dimension
@inline nbasis(el::Element, d::Integer) = nbasis(el)[d]

function evalBasis(el::Element, points, deriv::Integer=0)
    if deriv == 0
        return evalD0Basis(el, points)
    elseif deriv == 1
        return evalD1basis(el, points)
    elseif deriv == 2
        return evalD2Basis(el, points)
    else
        error("Invalid derivation order $deriv")
    end
end

function evalD0Basis{T<:Real}(el::Element, points::AbstractArray{T,2})
    nP, nD = size(points)
    nB = nbasis(el, nD)

    B = zeros(nB, nP)
    for i = 1:nP
        B[:,i] = evalD0Basis(el, points[i,:]...)
    end
    
    return B
end

function evalD1Basis{T<:Real}(el::Element, points::AbstractArray{T,2})
    nP, nD = size(points)
    nB = nbasis(el, nD)

    B = zeros(nB, nP, nD)
    for i = 1:nP
        B[:,i,:] = evalD1Basis(el, points[i,:]...)
    end
    
    return B
end

function evalD2Basis{T<:Real}(el::Element, points::AbstractArray{T,2})
    nP, nD = size(points)
    nB = nbasis(el, nD)

    B = zeros(nB, nP, nD, nD)
    for i = 1:nP
        B[:,i,:,:] = evalD2Basis(el, points[i,:]...)
    end
    
    return B
end

#-----------------#
# Type LagrangeP1 #
#-----------------#
type LagrangeP1 <: AbstractElement
end

order(::Element{LagrangeP1}) = 1
nbasis(el::Element{LagrangeP1}) = tuple(2:(dimension(el)+1)...)
nvertices(el::Element{LagrangeP1}) = tuple(2:(dimension(el)+1)...)

evalD0Basis(el::Element{LagrangeP1}, x::Real) = [1 - x, x]
evalD0Basis(el::Element{LagrangeP1}, x::Real, y::Real) = [1 - (x + y), x, y]    
evalD0Basis(el::Element{LagrangeP1}, x::Real, y::Real, z::Real) = [1 - (x + y + z), x, y, z]    

evalD1Basis(el::Element{LagrangeP1}, x::Real) = [-1, 1]
evalD1Basis(el::Element{LagrangeP1}, x::Real, y::Real) = [-1 -1; 1  0; 0 1]    
evalD1Basis(el::Element{LagrangeP1}, x::Real, y::Real, z::Real) = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1]

evalD2Basis(el::Element{LagrangeP1}, x::Real) = zeros(nbasis(el,1), dimension(el), dimension(el))
evalD2Basis(el::Element{LagrangeP1}, x::Real, y::Real) = zeros(nbasis(el,2), dimension(el), dimension(el))
evalD2Basis(el::Element{LagrangeP1}, x::Real, y::Real, z::Real) = zeros(nbasis(el,3), dimension(el), dimension(el))

end # of module Elements
