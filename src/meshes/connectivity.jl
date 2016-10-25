#-----------------------#
# Type MeshConnectivity #
#-----------------------#
"""
  Stores the incidence relation 'd -> dd' for a fixed pair
  of topological dimensions '(d,dd)'.
"""
type MeshConnectivity{T<:Integer} <: AbstractSparseMatrix{T,T}
    dims :: Tuple{Int, Int}
    indices :: Vector{T}
    offsets :: Vector{T}

    function MeshConnectivity(d1::Integer, d2::Integer, ind::Vector{T}, off::Vector{T})
        return new((Int(d1), Int(d2)), ind, off)
    end
end

# Outer Constructors
# -------------------
function MeshConnectivity(d1::Integer, d2::Integer, ind::Vector, off::Vector)
    T = promote_type(eltype(ind), eltype(off))
    return MeshConnectivity{T}(d1, d2, ind, off)
end
MeshConnectivity(dims::Tuple{Integer, Integer}, ind::Vector, off::Vector) =
    MeshConnectivity(dims..., ind, off)
MeshConnectivity(d1::Integer, d2::Integer) = MeshConnectivity(d1, d2, Int[], Int[])
MeshConnectivity(dims::Tuple{Integer,Integer}) = MeshConnectivity(dims...)

# Associated Methods
# -------------------
function Base.show(io::IO, ::MIME"text/plain", mc::MeshConnectivity)
    show(io, mc)
end

function Base.show(io::IO, mc::MeshConnectivity)
    print(io, "MeshConnectivity: ", mc.dims[1], " -> ", mc.dims[2])
end

@inline function Base.size(mc::MeshConnectivity)
    n_d  = isempty(mc.offsets) ? 0 : length(mc.offsets)-1
    n_dd = isempty(mc.indices) ? 0 : maximum(mc.indices)
    return (n_d, n_dd)
end
@inline Base.size(mc::MeshConnectivity, d::Integer) = size(mc)[d]

# implement some of the Base.SparseMatrixCSC methods
Base.nnz(mc::MeshConnectivity) = Int(mc.offsets[end]-1)

Base.nzrange(mc::MeshConnectivity, row::Integer) = mc.offsets[row]:(mc.offsets[row+1]-1)

function Base.find(mc::MeshConnectivity)
    sz = size(mc)
    I, J = findn(mc)
    return sub2ind(sz, I, J)
end

function Base.findn{T}(mc::MeshConnectivity{T})
    numnz = nnz(mc)
    I = Array{T}(numnz)
    J = Array{T}(numnz)

    count = 1
    @inbounds for row = 1:size(mc,1), r = nzrange(mc, row)
        I[count] = row
        J[count] = mc.indices[r]
        count += 1
    end
    count -= 1

    return (I, J)
end

# getindex
Base.getindex(mc::MeshConnectivity, I::Tuple{Integer,Integer}) = getindex(mc, I[1], I[2])
Base.getindex(mc::MeshConnectivity, i::Integer, ::Colon) = getindex(mc, i, nzrange(mc,i))
Base.getindex(mc::MeshConnectivity, ::Colon, j::Integer) = getindex(mc, 1:size(mc,1), j)

function Base.getindex(mc::MeshConnectivity, i::Integer, j::Integer)
    r1 = Int(mc.offsets[i])
    r2 = Int(mc.offsets[i+1])

    if !(1 <= i <= size(mc,1) && 1 <= j <= r2-r1); throw(BoundsError()); end

    return mc.indices[r1+j-1]
end

function Base.getindex{T<:Integer}(mc::MeshConnectivity{T}, I::AbstractVector, J::AbstractVector)
    m = size(mc, 1)
    nI = length(I)
    nJ = length(J)
    C = Array{T}(nI, nJ)

    for (ic, i) in enumerate(I)
        (1 <= i <= m) || throw(BoundsError(mc.offsets[1:end-1], i))
        r1 = Int(mc.offsets[i]-1)
        r2 = Int(mc.offsets[i+1]-1)
        for (jc, j) in enumerate(J)
            (1 <= j <= r2-r1) || throw(BoundsError(mc.indices[r1+1:r2], j))
            C[ic,jc] = mc.indices[r1+j]
        end
    end
    
    return C
end


# provide an iteration interface for `MeshConnectivity` instances.
Base.start(::MeshConnectivity) = 1

function Base.next(mc::MeshConnectivity, state::Integer)
    #return (mc.indices[mc.offsets[state]:mc.offsets[state+1]-1], state+1)
    return (view(mc.indices, mc.offsets[state]:mc.offsets[state+1]-1), state+1)
end

Base.done(mc::MeshConnectivity, state::Integer) = state > length(mc)
Base.eltype(::Type{MeshConnectivity}) = Array{Int,1}
Base.length(mc::MeshConnectivity) = size(mc, 1)

#Base.getindex(mc::MeshConnectivity, i::Integer) = mc.indices[mc.offsets[i]:mc.offsets[i+1]-1]
Base.getindex(mc::MeshConnectivity, i::Integer) = view(mc.indices, mc.offsets[i]:mc.offsets[i+1]-1)
Base.endof(mc::MeshConnectivity) = length(mc)

