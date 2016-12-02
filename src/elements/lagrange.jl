#-----------------#
# Type LagrangeP1 #
#-----------------#
type LagrangeP1 <: PElement
end
LagrangeP1(dim::Integer) = Element(LagrangeP1, dim)

isnodal(::Element{LagrangeP1}) = true
order(::Element{LagrangeP1}) = 1
nBasis(el::Element{LagrangeP1}) = tuple(2:(dimension(el)+1)...)
dofTuple(el::Element{LagrangeP1}) = (1, 0, 0, 0)[1:dimension(el)+1]

# Associated Methods
# -------------------
function evalD0Basis!{T<:Float}(el::Element{LagrangeP1}, points::AbstractArray{T,2}, out::Array{T,3})
    for ip = 1:size(points, 1)
        out[1,ip,:] = one(T)
        for id = 1:size(points, 2)
            out[1,ip,:] -= points[ip,id]
            out[id+1,ip,:] = points[ip,id]
        end
    end
    return out
end

function evalD1Basis!{T<:Float}(el::Element{LagrangeP1}, points::AbstractArray{T,2}, out::Array{T,4})
    out[1,:,:,:] = -one(T)
    for k = 1:size(out, 1)-1
        out[k+1,:,:,k] = one(T)
    end
    return out
end

function evalD2Basis!{T<:Float}(el::Element{LagrangeP1}, points::AbstractArray{T,2}, out::Array{T,5})
    out[:] = zero(T)
    return out
end

#-----------------#
# Type LagrangeQ1 #
#-----------------#
type LagrangeQ1 <: QElement
end
LagrangeQ1(dim::Integer) = Element(LagrangeQ1, dim)

isnodal(::Element{LagrangeQ1}) = true
order(::Element{LagrangeQ1}) = 1
nBasis(el::Element{LagrangeQ1}) = tuple(2.^(1:dimension(el))...)
dofTuple(el::Element{LagrangeQ1}) = (1, 0, 0, 0)[1:dimension(el)+1]

# Associated Methods
# -------------------
function evalD0Basis!{T<:Float}(el::Element{LagrangeQ1}, points::AbstractArray{T,2}, out::Array{T,3})
    nP, nD = size(points)
    I = one(T)
    if nD == 1
        for ip = 1:nP
            out[1,ip,:] = I - points[ip,1]
            out[2,ip,:] = points[ip,1]
        end
    elseif nD == 2
        for ip = 1:nP
            out[1,ip,:] = (I - points[ip,1]) * (I - points[ip,2])
            out[2,ip,:] = points[ip,1] * (I - points[ip,2])
            out[3,ip,:] = (I - points[ip,1]) * points[ip,2]
            out[4,ip,:] = points[ip,1] * points[ip,2]
        end
    elseif nD == 3
        for ip = 1:nP
            out[1,ip,:] = (I - points[ip,1]) * (I - points[ip,2]) * (1-points[ip,3])
            out[2,ip,:] = points[ip,1] * (I - points[ip,2]) * (1-points[ip,3])
            out[3,ip,:] = (I - points[ip,1]) * points[ip,2] * (I - points[ip,3])
            out[4,ip,:] = points[ip,1] * points[ip,2] * (I - points[ip,3])
            out[5,ip,:] = (I - points[ip,1]) * (I - points[ip,2]) * points[ip,3]
            out[6,ip,:] = points[ip,1] * (I - points[ip,2]) * points[ip,3]
            out[7,ip,:] = (I - points[ip,1]) * points[ip,2] * points[ip,3]
            out[8,ip,:] = points[ip,1] * points[ip,2] * points[ip,3]
        end
    else
        error("Invalid point dimension ", nD)
    end
    return out
end

function evalD1Basis!{T<:Float}(el::Element{LagrangeQ1}, points::AbstractArray{T,2}, out::Array{T,4})
    nP, nD = size(points)
    I = one(T)
    if nD == 1
        for ip = 1:nP
            out[1,ip,:,1] = -I
            out[2,ip,:,1] = I
        end
    elseif nD == 2
        for ip = 1:nP
            out[1,ip,:,1] = -(I - points[ip,2])
            out[1,ip,:,2] = -(I - points[ip,1])
            out[2,ip,:,1] = I - points[ip,2]
            out[2,ip,:,2] = -points[ip,1]
            out[3,ip,:,1] = -points[ip,2]
            out[3,ip,:,2] = I - points[ip,1]
            out[4,ip,:,1] = points[ip,2]
            out[4,ip,:,2] = points[ip,1]
        end
    elseif nD == 3
        for ip = 1:nP
            out[1,ip,:,1] = -(I - points[ip,2]) * (I - points[ip,3])
            out[1,ip,:,2] = -(I - points[ip,1]) * (I - points[ip,3])
            out[1,ip,:,3] = -(I - points[ip,1]) * (I - points[ip,2])
            
            out[2,ip,:,1] = (I - points[ip,2]) * (I - points[ip,3])
            out[2,ip,:,2] = -points[ip,1] * (I - points[ip,3])
            out[2,ip,:,3] = -points[ip,1] * (I - points[ip,2])
            
            out[3,ip,:,1] = -points[ip,2] * (I - points[ip,3])
            out[3,ip,:,2] = (1-points[ip,1]) * (I - points[ip,3])
            out[3,ip,:,3] = -(I - points[ip,1]) * points[ip,2]
            
            out[4,ip,:,1] = points[ip,2] * (I - points[ip,3])
            out[4,ip,:,2] = points[ip,1] * (I - points[ip,3])
            out[4,ip,:,3] = -points[ip,1] * points[ip,2]
            
            out[5,ip,:,1] = -(I - points[ip,2]) * points[ip,3]
            out[5,ip,:,2] = -(I - points[ip,1]) * points[ip,3]
            out[5,ip,:,3] = (I - points[ip,1]) * (I - points[ip,2])
            
            out[6,ip,:,1] = (I - points[ip,2]) * points[ip,3]
            out[6,ip,:,2] = -points[ip,1] * points[ip,3]
            out[6,ip,:,3] = points[ip,1] * (I - points[ip,2])
            
            out[7,ip,:,1] = -points[ip,2] * points[ip,3]
            out[7,ip,:,2] = (I - points[ip,1]) * points[ip,3]
            out[7,ip,:,3] = (I - points[ip,1]) * points[ip,2]
            
            out[8,ip,:,1] = points[ip,2] * points[ip,3]
            out[8,ip,:,2] = points[ip,1] * points[ip,3]
            out[8,ip,:,3] = points[ip,1] * points[ip,2]
        end
    else
        error("Invalid point dimension ", nD)
    end
    return out
end

function evalD2Basis!{T<:Float}(el::Element{LagrangeQ1}, points::AbstractArray{T,2}, out::Array{T,5})
  out[:] = zero(T)
  return out
end

