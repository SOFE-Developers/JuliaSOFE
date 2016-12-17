module Helpers

export dimension
export tensorprod, lagrangeNodesP, lagrangeNodesQ, ndarray

function dimension
end

function tensorprod{T}(gridx::AbstractArray{T,1})
    return hcat([i for i in gridx])
end

function tensorprod{T}(gridx::AbstractArray{T,1}, gridy::AbstractArray{T,1})
    return hcat([j for i in gridx for j in gridy],
                [i for i in gridx for j in gridy])
end

function tensorprod{T}(gridx::AbstractArray{T,1}, gridy::AbstractArray{T,1}, gridz::AbstractArray{T,1})
    return hcat([k for i in gridx for j in gridy for k in gridz],
                [j for i in gridx for j in gridy for k in gridz],
                [i for i in gridx for j in gridy for k in gridz])
end

"""
Return the interpolation nodes for the lagrange Pp element
of dimension `d` and polynomial order `p`.
"""
function lagrangeNodesP(d::Integer, p::Integer; sorted::Bool=false)
    @assert d in (1, 2, 3)
    @assert p >= 0

    if p == 0
        nodes = 1/3 * ones(1, d)
    else
        if d == 1
            ls = linspace(0, 1, p+1)
            nodes = hcat(ls[1], ls[end], ls[2:end-1])'
        elseif d == 2
            # 3 vertex nodes
            v = [0. 0.; 1. 0.; 0. 1.]
            # 3(p-1) edge nodes
            ls = linspace(0, 1, (p-1)+2)[2:end-1]
            e = vcat(hcat(ls, zeros(ls)),
                     hcat(zeros(ls), ls),
                     hcat(reverse(ls), ls))
            # (p-1)(p-2)/2 interior nodes
            i = tensorprod(ls, ls)
            h = 1/p
            i = i[sum(i, 2)[:] .< 1-h/2,:]

            nodes = vcat(v, e, i)
        elseif d == 3
            # 4 vertex nodes
            v = [0. 0. 0.; 1. 0. 0.; 0. 1. 0.; 0. 0. 1]
            # 6(p-1) edge nodes
            ls = linspace(0, 1, (p-1)+2)[2:end-1]
            zs = zeros(ls)
            rs = reverse(ls)
            e = vcat(hcat(ls, zs, zs),
                     hcat(zs, ls, zs),
                     hcat(zs, zs, ls),
                     hcat(rs, ls, zs),
                     hcat(rs, zs, ls),
                     hcat(zs, rs, ls))
            # (p-1)(p-2)(p-3)/3 interior nodes
            i = tensorprod(ls, ls, ls)
            h = 1/p
            i = i[sum(i, 2)[:] .< 1-h/2,:]

            nodes = vcat(v, e, i)
        end
    end

    if sorted
        nodes = sortrows(nodes, by=reverse, lt=lexless)
    end

    return nodes
end

function lagrangeNodesQ(d::Integer, p::Integer; sorted::Bool=false)
    @assert d in (1, 2, 3)
    @assert p >= 0

    if p == 0
        nodes = 0.5 * ones(1,d)
    else
        ls = linspace(0, 1, p+1)
        ls = vcat(ls[1], ls[end], ls[2:end-1])
        
        if d == 1
            nodes = ls''
        elseif d == 2
            # 4 vertex nodes
            vs = ls[1:2]
            v = tensorprod(vs, vs)

            # 4(p-1) edge nodes
            es = ls[3:end]
            zs = zeros(es)
            os = ones(es)
            e = vcat(hcat(es, zs),
                     hcat(zs, es),
                     hcat(os, es),
                     hcat(es, os))

            # (p-1)(p-1) interior nodes
            i = tensorprod(es, es)

            nodes = vcat(v, e, i)
        elseif d == 3
            # 8 vertex nodes
            vs = ls[1:2]
            v = tensorprod(vs, vs, vs)

            # 12(p-1) edge nodes
            es = ls[3:end]
            zs = zeros(es)
            os = ones(es)
            e = vcat(hcat(es, zs, zs), hcat(zs, es, zs), hcat(zs, zs, es),
                     hcat(os, es, zs), hcat(os, zs, es), hcat(es, os, zs),
                     hcat(zs, os, es), hcat(os, os, es), hcat(es, zs, os),
                     hcat(zs, es, os), hcat(os, es, os), hcat(es, os, os))

            # (p-1)(p-1)(p-1) interior nodes
            i = tensorprod(es, es, es)

            nodes = vcat(v, e, i)
        end
    end

    if sorted
        nodes = sortrows(nodes, by=reverse, lt=lexless)
    end
    
    return nodes
end    

#function ndarray{T<:Real}(A::AbstractArray{AbstractArray{T}})
function ndarray{T<:AbstractArray}(A::AbstractArray{T})
    sza = size(A)
    sza1 = size(A[1])
    szb = tuple(sza..., sza1...)
    #B = Array{T}(szb...)
    B = Array{eltype(A[1])}(szb...)

    for i in eachindex(A)
        ii = ind2sub(sza, i)
        for j in eachindex(A[i])
            k = tuple(ii..., ind2sub(sza1, j)...)
            B[k...] = A[i][j]
        end
    end

    return B
end
ndarray{T<:Real}(A::AbstractArray{T}) = A

end # of module Helpers
