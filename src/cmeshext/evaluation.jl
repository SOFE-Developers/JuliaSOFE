export evaluate

typealias Constant{T<:Real} Union{T, AbstractVector{T}, AbstractMatrix{T}}

function evaluate{T<:Constant,S<:Float}(c::T, points::AbstractArray{S,2}, m::Mesh)
    P = evalReferenceMaps(m, points)
    nE, nP, nW = size(P)

    if issubtype(T, Real)
        C = convert(promote_type(eltype(T),S), c)
    elseif issubtype(T, AbstractVector)
        C = convert(Vector{promote_type(eltype(T),S)}, c)
    elseif issubtype(T, AbstractMatrix)
        C = convert(Matrix{promote_type(eltype(T),S)}, c)
    end
    
    R = zeros(S, nE, nP, size(C)...)
    for ip = 1:nP
        for ie = 1:nE
            R[ie,ip,:] = C
        end
    end
    return R
end

function evaluate{T<:Float}(f::Function, points::AbstractArray{T,2}, m::Mesh)
    P = evalReferenceMaps(m, points)
    nE, nP, nW = size(P)

    R = f(reshape(P, (nE*nP,nW)))
    nF = size(R, (2,(3:ndims(R))...)...)
    R = (nF == 1) ? reshape(R, nE, nP) : reshape(R, nE, nP, nF...)

    return ndarray(R)
end
