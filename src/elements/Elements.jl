__precompile__()

module Elements

export AbstractElement, LagrangeP1, LagrangeQ1

abstract AbstractElement

typealias Float AbstractFloat

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
                                     points::AbstractArray{T,2},
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

# Associated Methods
# -------------------
function evalD0Basis{T<:Float}(el::LagrangeQ1, points::AbstractArray{T,2})
    nP, nD = size(points)
    size(el.super, nD)
    B = zeros(nB, nP, 1)
    if nD == 1
        B[1,:] = 1-points[:,1];
        B[2,:] = points[:,1];
    elseif nD == 2
        B[1,:] = (1-points[:,1]).*(1-points[:,2]);
        B[2,:] = points[:,1].*(1-points[:,2]);
        B[3,:] = (1-points[:,1]).*points[:,2];
        B[4,:] = points[:,1].*points[:,2];
    elseif nD == 3
        B[1,:] = (1-points[:,1]).*(1-points[:,2]).*(1-points[:,3]);
        B[2,:] = points[:,1].*(1-points[:,2]).*(1-points[:,3]);
        B[3,:] = (1-points[:,1]).*points[:,2].*(1-points[:,3]);
        B[4,:] = points[:,1].*points[:,2].*(1-points[:,3]);
        B[5,:] = (1-points[:,1]).*(1-points[:,2]).*points[:,3];
        B[6,:] = points[:,1].*(1-points[:,2]).*points[:,3];
        B[7,:] = (1-points[:,1]).*points[:,2].*points[:,3];
        B[8,:] = points[:,1].*points[:,2].*points[:,3];
    else
        error("Not supported size of points")
    end
    return B
end

function evalD1Basis{T<:Float}(el::LagrangeQ1, points::AbstractArray{T,2})
    (nP, nD) = size(points)
    B = zeros(el.super.nB[nD], nP, 1, nD)
    if nD == 1
        B[1,:,:,1] = -1;
        B[2,:,:,1] = 1;
    elseif nD == 2
        B[1,:,:,1] = -(1-points[:,2]);
        B[1,:,:,2] = -(1-points[:,1]);
        B[2,:,:,1] = 1-points[:,2];
        B[2,:,:,2] = -points[:,1];
        B[3,:,:,1] = -points[:,2];
        B[3,:,:,2] = 1-points[:,1];
        B[4,:,:,1] = points[:,2];
        B[4,:,:,2] = points[:,1];
    elseif nD == 3
        B[1,:,:,1] = -(1-points[:,2]).*(1-points[:,3]);
        B[1,:,:,2] = -(1-points[:,1]).*(1-points[:,3]);
        B[1,:,:,3] = -(1-points[:,1]).*(1-points[:,2]);

        B[2,:,:,1] = (1-points[:,2]).*(1-points[:,3]);
        B[2,:,:,2] = -points[:,1].*(1-points[:,3]);
        B[2,:,:,3] = -points[:,1].*(1-points[:,2]);

        B[3,:,:,1] = -points[:,2].*(1-points[:,3]);
        B[3,:,:,2] = (1-points[:,1]).*(1-points[:,3]);
        B[3,:,:,3] = -(1-points[:,1]).*points[:,2];

        B[4,:,:,1] = points[:,2].*(1-points[:,3]);
        B[4,:,:,2] = points[:,1].*(1-points[:,3]);
        B[4,:,:,3] = -points[:,1].*points[:,2];

        B[5,:,:,1] = -(1-points[:,2]).*points[:,3];
        B[5,:,:,2] = -(1-points[:,1]).*points[:,3];
        B[5,:,:,3] = (1-points[:,1]).*(1-points[:,2]);
        
        B[6,:,:,1] = (1-points[:,2]).*points[:,3];
        B[6,:,:,2] = -points[:,1].*points[:,3];
        B[6,:,:,3] = points[:,1].*(1-points[:,2]);
        
        B[7,:,:,1] = -points[:,2].*points[:,3];
        B[7,:,:,2] = (1-points[:,1]).*points[:,3];
        B[7,:,:,3] = (1-points[:,1]).*points[:,2];
        
        B[8,:,:,1] = points[:,2].*points[:,3];
        B[8,:,:,2] = points[:,1].*points[:,3];
        B[8,:,:,3] = points[:,1].*points[:,2];
    else
        error("Not supported size of points")
    end
    return B
end

function evalD2Basis{T<:Float}(el::LagrangeQ1, points::AbstractArray{T,2})
  (nP, nD) = size(points)
  B = zeros(el.super.nB[nD], nP, 1, nD, nD)
  return B
end


end # of module Elements

