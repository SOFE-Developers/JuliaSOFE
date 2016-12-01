__precompile__()

module Quadrature

export AbstractQuadRule, QuadRule
export qorder, qnodes, qweights, quadData

typealias Float AbstractFloat

abstract AbstractQuadRule
abstract QuadRuleSimp <: AbstractQuadRule
abstract QuadRuleOrth <: AbstractQuadRule

#---------------#
# Type QuadRule #
#---------------#
type QuadRule{T<:AbstractQuadRule} <: AbstractQuadRule
    order :: Int
    nodes :: Tuple
    weights :: Tuple
end

function QuadRule{T<:AbstractQuadRule}(::Type{T}, order::Integer,
                                       points, weights)
    return QuadRule{T}(order, points, weights)
end

# Associated Methods
# -------------------
"""

    qorder(qr::AbstractQuadRule)

  Return the maximum polynomial order for which the quadrature rule is exact.
"""
qorder(qr::AbstractQuadRule) = getfield(qr, :order)

"""

    qnodes(qr::AbstractQuadRule, d::Integer)

  Return the quadrature nodes for the entities of
  topological dimension `d`.
"""
qnodes(qr::AbstractQuadRule, d::Integer) = getfield(qr, :nodes)[d]

"""

    qweights(qr::AbstractQuadRule, d::Integer)

  Return the quadrature weights for the entities of
  topological dimension `d`.
"""
qweights(qr::AbstractQuadRule, d::Integer) = getfield(qr, :weights)[d]

"""

    quadData(qr::AbstractQuadRule, d::Integer)

  Return the quadrature points and weights for the entities
  of topological dimension `d`.
"""
quadData(qr::AbstractQuadRule, d::Integer) = (qnodes(qr, d), qweights(qr, d))

# Gauss Quadrature
# -----------------
include("gauss.jl")

# #------------#
# # Type Simp2 #
# #------------#
# type Simp2 <: AbstractQuadRule
# end
# typealias QuadRuleSimp2 QuadRule{Simp2}

# function QuadRuleSimp2()
#     order = 2
#     points = ([0.11270167,  0.5       ,  0.88729833]'',
#               [0.16666667  0.16666667;
#                0.66666667  0.16666667;
#                0.16666667  0.66666667],
#               [0.1381966  0.1381966  0.1381966; 
#                0.5854102  0.1381966  0.1381966;
#                0.1381966  0.5854102  0.1381966;
#                0.1381966  0.1381966  0.5854102])
#     weights = ([0.27777778, 0.44444444, 0.27777778],
#                [0.16666667, 0.16666667, 0.16666667],
#                [0.04166667, 0.04166667, 0.04166667, 0.04166667])
    
#     return QuadRule(Simp2, order, points, weights)
# end


end # of module Quadrature
