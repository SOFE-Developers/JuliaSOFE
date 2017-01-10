# Lobatto Shape Functions
# ------------------------
include("shapefun.jl")

#----------#
# Type PpH #
#----------#
type PpH{p} <: PElement
end
PpH(dim::Integer, p::Integer) = Element(PpH{p}, dim)

ishierarchic(::Element{PpH}) = true
order{p}(::Element{PpH{p}}) = Int(p)
function nBasis{p}(el::Element{PpH{p}})
    nB = (p+1, div((p+1)*(p+2), 2), div((p+1)*(p+2)*(p+3), 6) + div(4*2*(p-1)*(p-2),2))
    return tuple(nB[1:dimension(el)]...)
end
function dofTuple{p}(el::Element{PpH{p}})
    dT = ones(Int, dimension(el)+1)
    for i = 0:dimension(el)
        for j = 1:i
            dT[i+1] *= (p - j)
        end
        dT[i+1] /= factorial(i)
    end
    return dT
end

function evalD0Basis!{T<:Real,p}(el::Element{PpH{p}}, points::AbstractArray{T,2}, out::Array{T,3})
    nP, nD = size(points)
    #p = order(el)

    if nD == 1
        for ip = 1:nP
            for ib = 1:p+1
                out[ib,ip,1] = poly_lobatto(ib-1)(2*points[ip,1] - 1)
            end
        end
    elseif nD == 2
        for ip = 1:nP
            # barycentric coordinates
            l1 = one(T) - (points[ip,1] + points[ip,2]); l2 = points[ip,1]; l3 = points[ip,2]
            # vertex functions
            out[1,ip,1] = l1; out[2,ip,1] = l2; out[3,ip,1] = l3
            # edge functions (2 <= k <= pe)
            off = 4
            for k = 0:p-2
                out[off,ip,1] = l1 * l2 * poly_lobatto_kernel(k)(l2 - l1); off += 1
                out[off,ip,1] = l1 * l3 * poly_lobatto_kernel(k)(l3 - l1); off += 1
                out[off,ip,1] = l2 * l3 * poly_lobatto_kernel(k)(l3 - l2); off += 1
            end
            # interior functions (1 <= n1,n2; n1+n2 <= pi-1)
            for i = 3:p
                for j = 0:i-3
                    out[off,ip,1] = l1 * l2 * l3 *
                        poly_lobatto_kernel(j)(l2 - l1) * poly_lobatto_kernel(i-j-3)(l3 - l1)
                    off += 1
                end
            end
        end
    elseif nD == 3
        error("Not implemented, yet!")
    end

    return out
end
