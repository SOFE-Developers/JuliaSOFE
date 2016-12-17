using Jacobi

using ..Helpers

#-------------------------#
# Newton Cotes Quadrature #
#-------------------------#
type NewtonCotes <: AbstractQuadRule
end

function NewtonCotesQuadOrth(order::Integer, dim::Integer)
    
end

function newtoncotesquadorth(order::Integer, dim::Integer)
    nodes = lagrangeNodesQ(dim, order, sorted=true)
end
