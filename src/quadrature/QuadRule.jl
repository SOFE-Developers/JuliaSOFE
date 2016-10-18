__precompile__()

module Quadrature

abstract AbstractQuadRuleX

#----------------#
# Type QuadRuleX #
#----------------#
type QuadRuleX <: AbstractQuadRuleX
    order :: Integer
    points :: Array{AbstractFloat, 2}
    weights :: Array{AbstractFloat, 1}

    function QuadRule{T<:AbstractFloat}(dim::Integer, points::AbstractArray{T,2},
                                        weights::AbstractArray{T,1})
        return new(dim, points, weights)
    end
end

# Associated Methods
# -------------------

getQuadData(qr::AbstractQuadRuleX) = (qr.points, qr.weights)

################################################################################
# Concrete Types

type QuadTriNodal <: AbstractQuadRuleX
    super :: QuadRuleX
    
    function QuadTriNodal(order::Int)
        if order > 1
            warn("Quad order $order cannot be satisfied")
        end
        return new(QuadRule(1, [0.0 0.0; 1.0 0.0; 0.0 1.0], [1,1,1]/6))
    end
end

type QuadTriNodal2 <: AbstractQuadRuleX
    super :: QuadRuleX

    function QuadTriNodal2(order::Int)
        if order > 2
            warn("Quad order $order cannot be satisfied")
        end
        return new(QuadRule(2, [0.5 0.0; 0.5 0.5; 0.0 0.5], [1,1,1]/6))
    end
end

type QuadTriGauss6 <: AbstractQuadRuleX
    super :: QuadRuleX

    function QuadTriGauss6(order::Int)
        if order > 6
            warn("Quad order $order cannot be satisfied")
        end
        return new(QuadRule(6, [0.659027622374092  0.231933368553031;
                                0.659027622374092  0.109039009072877;
                                0.231933368553031  0.659027622374092;
                                0.231933368553031  0.109039009072877;
                                0.109039009072877  0.659027622374092;
                                0.109039009072877  0.231933368553031], [1,1,1,1,1,1]/12))
    end
end

type QuadIntNodal <: AbstractQuadRuleX
    super :: QuadRuleX

    function QuadIntNodal(order::Int)
        if order > 1
            warn("Quad order $order cannot be satisfied")
        end
        return new(QuadRule(1, reshape([0.0; 1.0],2,1), [1,1]/2))
    end
end
################################################################################

end # of module Quadrature
