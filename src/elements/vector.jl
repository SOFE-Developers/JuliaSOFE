export VectorElement
export subelem, components

#--------------------#
# Type VectorElement #
#--------------------#
type VectorElement{T<:ElementTypes} <: AbstractElement{T}
    subelem :: Element{T}
    components :: Int
end

# Outer Constructors
# -------------------
function VectorElement{T<:ElementTypes}(el::T, cmp::Integer)
    return VectorElement{T}(el, Int(cmp))
end

function VectorElement{T<:ElementTypes}(::Type{T}, dim::Integer, cmp::Integer)
    return VectorElement(Element{T}(dim), cmp)
end

# Associated Methods
# -------------------
subelem(el::VectorElement) = getfield(el, :subelem)
components(el::VectorElement) = getfield(el, :components)

dimension(el::VectorElement) = dimension(subelem(el))

isnodal(el::VectorElement) = isnodal(subelem(el))

order(el::VectorElement) = order(subelem(el))

nBasis(el::VectorElement) = tuple((components(el) * nb for nb in nBasis(subelem(el)))...)

dofTuple(el::VectorElement) = tuple((components(el) * nd for nd in dofTuple(subelem(el)))...)

function evalBasis{T<:AbstractFloat}(el::VectorElement, points::AbstractArray{T,2}, deriv::Integer=0)
    B = evalBasis(subelem(el), points, deriv)
    nC = components(el)
    
    if deriv == 0
        return repeat(B, outer=(1,1,nC))
    elseif deriv == 1
        return repeat(B, outer=(1,1,nC,1))
    elseif deriv == 0
        return repeat(B, outer=(1,1,nC,1,1))
    end
end

