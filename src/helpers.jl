module Helpers

using CMesh.Utils

export lagrangeNodesP, lagrangeNodesQ

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
            if p == 1
                nodes = hcat(ls[1], ls[end])'
            else
                nodes = hcat(ls[1], ls[end], ls[2:end-1])'
            end
        elseif d == 2
            # 3 vertex nodes
            v = [0. 0.; 1. 0.; 0. 1.]
            if p > 1
                # 3(p-1) edge nodes
                ls = linspace(0, 1, (p-1)+2)[2:end-1]
                e = vcat(hcat(ls, zeros(ls)),
                         hcat(zeros(ls), ls),
                         hcat(reverse(ls), ls))
                # (p-1)(p-2)/2 interior nodes
                i = tensorprod(ls, ls)
                h = 1/p
                i = i[sum(i, 2)[:] .< 1-h/2,:]
            else
                e = i = zeros(0, d)
            end

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

function ndomain(f::Function)
    n = 0
    while n < 10
        try
            x = rand(n)
            y = f(x)
            return n
        catch
            n += 1
        end
    end
end

function ncodomain(f::Function)
    n = ndomain(f)
    return size(f(rand(n)))
end

dims(f::Function) = (ndomain(f), ncodomain(f))

end # of module Helpers
