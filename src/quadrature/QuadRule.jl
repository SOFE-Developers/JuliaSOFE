__precompile__()

module Quadrature

export AbstractQuadRule, QuadRule, QuadRuleSimp2
export qorder, qpoints, qweights, quadData

typealias Float AbstractFloat

abstract AbstractQuadRule

#---------------#
# Type QuadRule #
#---------------#
type QuadRule{T<:AbstractQuadRule} <: AbstractQuadRule
    order :: Int
    points :: Tuple
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
@inline qorder(qr::AbstractQuadRule) = qr.order

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
# Type Simp2 #
#------------#
type Simp2 <: AbstractQuadRule
end
typealias QuadRuleSimp2 QuadRule{Simp2}

function QuadRuleSimp2()
    order = 2
    points = ([0.11270167,  0.5       ,  0.88729833]'',
              [0.16666667  0.16666667;
               0.66666667  0.16666667;
               0.16666667  0.66666667],
              [0.1381966  0.1381966  0.1381966; 
               0.5854102  0.1381966  0.1381966;
               0.1381966  0.5854102  0.1381966;
               0.1381966  0.1381966  0.5854102])
    weights = ([0.27777778, 0.44444444, 0.27777778],
               [0.16666667, 0.16666667, 0.16666667],
               [0.04166667, 0.04166667, 0.04166667, 0.04166667])
    
    return QuadRule(Simp2, order, points, weights)
end


end # of module Quadrature
