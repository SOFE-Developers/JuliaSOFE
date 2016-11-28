export MixedElement
export subelements

#-------------------#
# Type MixedElement #
#-------------------#
type MixedElement{N} <: AbstractElement
    dimension :: Int
    subelements :: NTuple{N,Element}

    function MixedElement{N}(es::NTuple{N,Element})
        dim = dimension(es[1])
        @assert all(dimension(e) == dim for e in es)
        return new(dim, es)
    end
end
MixedElement(es::Element...) = MixedElement{length(es)}(es)

# Associated Methods
# -------------------
dimension(el::MixedElement) = getfield(el, :dimension)
subelements(el::MixedElement) = getfield(el, :subelements)

function evalBasis{T<:AbstractFloat}(el::MixedElement, points::AbstractArray{T,2}, deriv::Integer=0)
    return tuple((evalBasis(e, points, deriv) for e in el)...)
end

# Iteration Interface
# --------------------
Base.getindex(el::MixedElement, i::Integer) = getindex(subelements(el), i)
Base.length(el::MixedElement) = length(subelements(el))
Base.eltype(::MixedElement) = Element

Base.start(::MixedElement) = 1
Base.next(el::MixedElement, status::Integer) = (el[status], status+1)
Base.done(el::MixedElement, status::Integer) = status > length(el)

