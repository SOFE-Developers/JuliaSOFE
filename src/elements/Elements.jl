__precompile__()

module Elements

export AbstractElement, Element, LagrangeP1, LagrangeQ1
export dimension, order, nBasis, nVertices, dofTuple, nDoF, evalBasis

abstract AbstractElement
abstract PElement <: AbstractElement
abstract QElement <: AbstractElement

typealias Float AbstractFloat

#--------------#
# Type Element #
#--------------#
type Element{E<:AbstractElement} <: AbstractElement
    dimension :: Int
end
Element{E<:AbstractElement}(::Type{E}, dim::Integer) = Element{E}(dim)

# Associated methods
# -------------------
"""

    dimension(el::Element)

  The topological dimension of the element.
"""
@inline dimension(el::Element) = el.dimension

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
@inline nVertices(el::AbstractElement, d::Integer) = nVertices(el)[d]
@inline nVertices{E<:PElement}(el::Element{E}) = tuple(2:(dim(el)+1)...)
@inline nVertices{E<:QElement}(el::Element{E}) = tuple(2.^(1:dim(el))...)

function nDoF{E<:PElement}(el::Element{E})
    return map(*, dofTuple(el), binomial(dim(el)+1, k) for k = 1:dim(el)+1)
end
function nDoF{E<:QElement}(el::Element{E})
    if dim(el) == 1
        return map(*, dofTuple(el), (2,1))
    elseif dim(el) == 2
        return map(*, dofTuple(el), (4,4,1))
    elseif dim(el) == 3
        return map(*, dofTuple(el), (8,12,6,1))
    else
        error("Invalid element dimension ", dim(el))
    end
end
@inline nDoF(el::AbstractElement, d::Integer) = nDoF(el)[d]

"""

    evalBasis{T<:AbstractFloat}(el::Element, points::AbstractArray{T,1}, deriv::Integer=0)

  Evaluate the element's shape (basis) functions at given `points`
  on the reference domain.
"""
function evalBasis{T<:AbstractFloat}(el::Element,
                                     points::AbstractArray{T,2},
                                     deriv::Integer=0)
    nP, nD = size(points)
    nB = nBasis(el, nD)

    if deriv == 0
        B = zeros(T, nB, nP)
        return evalD0Basis!(el, points, B) # nB x nP
    elseif deriv == 1
        B = zeros(T, nB, nP, nD)
        return evalD1Basis!(el, points, B) # nB x nP x nD
    elseif deriv == 2
        B = zeros(T, nB, nP, nD, nD)
        return evalD2Basis!(el, points, B) # nB x nP x nD x nD
    else
        error("Invalid derivation order! ($deriv)")
    end
end

#-----------------#
# Type LagrangeP1 #
#-----------------#
type LagrangeP1 <: PElement
    LagrangeP1(dim::Integer) = Element(LagrangeP1, dim)
end

order(::Element{LagrangeP1}) = 1
nBasis(el::Element{LagrangeP1}) = tuple(2:(dim(el)+1)...)
dofTuple(el::Element{LagrangeP1}) = (1, 0, 0, 0)[1:dim(el)+1]

# Associated Methods
# -------------------
function evalD0Basis!{T<:Float}(el::Element{LagrangeP1}, points::AbstractArray{T,2}, out::Array{T,2})
    out[1,:]     = 1 - sum(points, 2)
    out[2:end,:] = points'
    return out
end

function evalD1Basis!{T<:Float}(el::Element{LagrangeP1}, points::AbstractArray{T,2}, out::Array{T,3})
    out[1,:,:] = -1
    for k = 1:size(out, 1)-1
        out[k+1,:,k] = 1
    end
    return out
end

function evalD2Basis!{T<:Float}(el::Element{LagrangeP1}, points::AbstractArray{T,2}, out::Array{T,4})
    out[:] = 0
    return out
end

#-----------------#
# Type LagrangeQ1 #
#-----------------#
type LagrangeQ1 <: QElement
    LagrangeQ1(dim::Integer) = Element(LagrangeQ1, dim)
end

order(::Element{LagrangeQ1}) = 1
nBasis(el::Element{LagrangeQ1}) = tuple(2.^(1:dim(el))...)
dofTuple(el::Element{LagrangeQ1}) = (1, 0, 0, 0)[1:dim(el)+1]

# Associated Methods
# -------------------
function evalD0Basis!{T<:Float}(el::Element{LagrangeQ1}, points::AbstractArray{T,2}, out::Array{T,3})
    nD = size(points, 2)
    if nD == 1
        out[1,:] = 1-points[:,1];
        out[2,:] = points[:,1];
    elseif nD == 2
        out[1,:] = (1-points[:,1]).*(1-points[:,2]);
        out[2,:] = points[:,1].*(1-points[:,2]);
        out[3,:] = (1-points[:,1]).*points[:,2];
        out[4,:] = points[:,1].*points[:,2];
    elseif nD == 3
        out[1,:] = (1-points[:,1]).*(1-points[:,2]).*(1-points[:,3]);
        out[2,:] = points[:,1].*(1-points[:,2]).*(1-points[:,3]);
        out[3,:] = (1-points[:,1]).*points[:,2].*(1-points[:,3]);
        out[4,:] = points[:,1].*points[:,2].*(1-points[:,3]);
        out[5,:] = (1-points[:,1]).*(1-points[:,2]).*points[:,3];
        out[6,:] = points[:,1].*(1-points[:,2]).*points[:,3];
        out[7,:] = (1-points[:,1]).*points[:,2].*points[:,3];
        out[8,:] = points[:,1].*points[:,2].*points[:,3];
    else
        error("Invalid point dimension ", nD)
    end
    return out
end

function evalD1Basis!{T<:Float}(el::Element{LagrangeQ1}, points::AbstractArray{T,2}, out::Array{T,4})
    nD = size(points, 2)
    if nD == 1
        out[1,:,:,1] = -1;
        out[2,:,:,1] = 1;
    elseif nD == 2
        out[1,:,:,1] = -(1-points[:,2]);
        out[1,:,:,2] = -(1-points[:,1]);
        out[2,:,:,1] = 1-points[:,2];
        out[2,:,:,2] = -points[:,1];
        out[3,:,:,1] = -points[:,2];
        out[3,:,:,2] = 1-points[:,1];
        out[4,:,:,1] = points[:,2];
        out[4,:,:,2] = points[:,1];
    elseif nD == 3
        out[1,:,:,1] = -(1-points[:,2]).*(1-points[:,3]);
        out[1,:,:,2] = -(1-points[:,1]).*(1-points[:,3]);
        out[1,:,:,3] = -(1-points[:,1]).*(1-points[:,2]);

        out[2,:,:,1] = (1-points[:,2]).*(1-points[:,3]);
        out[2,:,:,2] = -points[:,1].*(1-points[:,3]);
        out[2,:,:,3] = -points[:,1].*(1-points[:,2]);

        out[3,:,:,1] = -points[:,2].*(1-points[:,3]);
        out[3,:,:,2] = (1-points[:,1]).*(1-points[:,3]);
        out[3,:,:,3] = -(1-points[:,1]).*points[:,2];

        out[4,:,:,1] = points[:,2].*(1-points[:,3]);
        out[4,:,:,2] = points[:,1].*(1-points[:,3]);
        out[4,:,:,3] = -points[:,1].*points[:,2];

        out[5,:,:,1] = -(1-points[:,2]).*points[:,3];
        out[5,:,:,2] = -(1-points[:,1]).*points[:,3];
        out[5,:,:,3] = (1-points[:,1]).*(1-points[:,2]);
        
        out[6,:,:,1] = (1-points[:,2]).*points[:,3];
        out[6,:,:,2] = -points[:,1].*points[:,3];
        out[6,:,:,3] = points[:,1].*(1-points[:,2]);
        
        out[7,:,:,1] = -points[:,2].*points[:,3];
        out[7,:,:,2] = (1-points[:,1]).*points[:,3];
        out[7,:,:,3] = (1-points[:,1]).*points[:,2];
        
        out[8,:,:,1] = points[:,2].*points[:,3];
        out[8,:,:,2] = points[:,1].*points[:,3];
        out[8,:,:,3] = points[:,1].*points[:,2];
    else
        error("Invalid point dimension ", nD)
    end
    return out
end

function evalD2Basis!{T<:Float}(el::Element{LagrangeQ1}, points::AbstractArray{T,2}, out::Array{T,5})
  out[:] = 0
  return out
end

end # of module Elements

