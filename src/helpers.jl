function tensorprod{T}(gridx::AbstractArray{T,1}, gridy::AbstractArray{T,1})
    return hcat([i for i in gridx for j in gridy],
                [j for i in gridx for j in gridy])
end

function tensorprod{T}(gridx::AbstractArray{T,1}, gridy::AbstractArray{T,1}, gridz::AbstractArray{T,1})
    return hcat([i for i in gridx for j in gridy for k in gridz],
                [j for i in gridx for j in gridy for k in gridz],
                [k for i in gridx for j in gridy for k in gridz])
end

"""
Return the interpolation nodes for the lagrange Pp element
of dimension `d` and polynomial order `p`.
"""
function lagrangeNodesP(d::Integer, p::Integer)
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
end
