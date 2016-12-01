using FastGaussQuadrature

using ..Helpers

export GaussQuadOrth, GaussQuadSimp
export gaussquadorth, gaussquadsimp

#---------------------#
# Gauss Legendre Type #
#---------------------#
type Gauss <: AbstractQuadRule
end

function GaussQuadOrth(order::Integer, dim::Integer)
    nw = (gaussquadorth(order, d) for d = 1:dim)
    nodes, weights = collect(zip(nw...))

    return QuadRule{Gauss}(order, nodes, weights)
end

function GaussQuadSimp(order::Integer, dim::Integer)
    nw = (gaussquadsimp(order, d) for d = 1:dim)
    nodes, weights = collect(zip(nw...))

    return QuadRule{Gauss}(order, nodes, weights)
end

function gaussquadorth(order::Integer, dim::Integer, ntype::Symbol=:legendre)
    n = Int(ceil((order+1)/2))

    if ntype == :chebyshev
        nodes, weights = gausschebyshev(n)
    elseif ntype == :hermite
        nodes, weights = gausshermite(n)
    elseif ntype == :jacobi
        nodes, weights = gaussjacobi(n)
    elseif ntype == :laguerre
        nodes, weights = gausslaguerre(n)
    elseif ntype == :legendre
        nodes, weights = gausslegendre(n)
    elseif ntype == :lobatto
        nodes, weights = gausslobatto(n)
    elseif ntype == :radau
        nodes, weights = gaussradau(n)
    else
        error("Invalid quadrature node type: ", ntype)
    end

    nodes *= 0.5; nodes += 0.5
    weights *= 0.5

    if dim == 1
        nodes = tensorprod(nodes)
    elseif dim == 2
        nodes = tensorprod(nodes, nodes)
        weights = prod(tensorprod(weights, weights), 2)[:]
    elseif dim == 3
        nodes = tensorprod(nodes, nodes, nodes)
        weights = prod(tensorprod(weights, weights, weights), 2)[:]
    end

    return nodes, weights
end

function gaussquadsimp(order::Integer, dim::Integer, ntype::Symbol=:legendre)
    nodes, weights = gaussquadorth(order, dim, ntype)

    if dim == 1
        #...
    elseif dim == 2
        I = one(eltype(nodes))
        weights .*= (I - nodes[:,2])
        nodes[:,1] .*= (I - nodes[:,2])
    elseif dim == 3
        error("Not implemented, yet!")
    end

    return nodes, weights
end
