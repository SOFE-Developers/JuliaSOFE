__precompile__()

module Quadrature

typealias Float AbstractFloat

export AbstractQuadRule, QuadRule, QuadRuleSimp1
export order, qpoints, qweights, quadData

abstract AbstractQuadRule

#---------------#
# Type QuadRule #
#---------------#
type QuadRule{X<:AbstractQuadRule} <: AbstractQuadRule
    order :: Int
    points :: Tuple
    weights :: Tuple
end

function QuadRule{X<:AbstractQuadRule}(::Type{X}, order::Integer,
                                       points, weights)
    return QuadRule{X}(order, points, weights)
end

# Associated Methods
# -------------------
"""

    order(qr::AbstractQuadRule)

  Return the maximum polynomial order for which the quadrature rule is exact.
"""
@inline order(qr::AbstractQuadRule) = qr.order

"""

    qpoints(qr::AbstractQuadRule, d::Integer)

  Return the quadrature points for the entities of
  topological dimension `d`.
"""
@inline qpoints(qr::AbstractQuadRule, d::Integer) = qr.points[d]

"""

    qweights(qr::AbstractQuadRule, d::Integer)

  Return the quadrature weights for the entities of
  topological dimension `d`.
"""
@inline qweights(qr::AbstractQuadRule, d::Integer) = qr.weights[d]

"""

    quadData(qr::AbstractQuadRule, d::Integer)

  Return the quadrature points and weights for the entities
  of topological dimension `d`.
"""
quadData(qr::AbstractQuadRule, d::Integer) = (qr.points[d], qr.weights[d])

#------------#
# Type Simp1 #
#------------#
type Simp1 <: AbstractQuadRule
end
typealias QuadRuleSimp1 QuadRule{Simp1}

function QuadRuleSimp1()
    order = 1
    points = ([0.0; 1.0],
              [0.0 0.0; 1.0 0.0; 0.0 1.0],
              [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    weights = ([1.0, 1.0]./6,
               [1.0, 1.0, 1.0]./6,
               [1.0, 1.0, 1.0, 1.0]./6)
    return QuadRule(Simp1, order, points, weights)
end


end # of module Quadrature
