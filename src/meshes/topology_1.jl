
#--------------------#
# Type MeshTopologyX #
#--------------------#
type MeshTopologyX <: AbstractMeshTopology
    dimension :: Integer
    nodes :: Array{AbstractFloat, 2}
    connectivities :: Dict{Tuple{Int, Int}, Array{Int,2}}

    function MeshTopologyX{Tn<:AbstractFloat, Te<:Integer}
        (nodes::AbstractArray{Tn,2}, elems::AbstractArray{Te,2})
        dim = size(nodes, 2)
        connect = Dict((dim,0) => elems)
        return new(dim, nodes, connect)
    end
end

# Associated Methods
# -------------------
function getConnectivity(mt::MeshTopologyX, d::Integer, dd::Integer)
    return mt.super.connectivities[(d,dd)]
end

function setConnectivity(mt::MeshTopologyX, d::Integer, dd::Integer, connect::Array{Int, 2})
    mt.super.connectivities[(d,dd)] = connect;
end

function getEntities(mt::MeshTopologyX, d::Integer)
    if dim == 0
        return mt.super.nodes
    else
        return mt.super.connectivities[(d,0)]
    end
end

function getNumber(mt::MeshTopologyX, d::Integer)
    return size(getEntities(mt, d), 1)
end

function isBoundary(mt::MeshTopologyX, f=x->x[:,1].>Inf)
    e2F = mt.super.connectivities[(2,1)][:]
    R = (sparse(e2F, Array{Int64}(ones(length(e2F))), Array{Int64}(ones(length(e2F)))).==1)[:];
    face = getEntities(mt, 1)
    center = permutedims(sum(mt.super.nodes[face,:],2)/size(face,2), [1 3 2])
    R = R & where(center)
    return R
end


getNodes(mt::MeshTopologyX) = mt.super.nodes
getDim(mt::MeshTopologyX) = mt.super.dimension

#----------------------#
# Type MeshTopologyTri #
#----------------------#
type MeshTopologyTri <: AbstractMeshTopology
    super :: MeshTopologyX

    function MeshTopologyTri{Tn<:AbstractFloat, Te::Integer}
        (nodes::AbstractArray{Tn,2}, elems::AbstractArray{Te,2})
        mt = new(MeshTopologyX(nodes, elems))
        updateConnectivity(mt)
        return mt
    end
end

# Associated Methods
# -------------------
function updateConnectivity(mt::MeshTopologyTri)
    D = getDim(mt)
    setConnectivity(mt, 0, 0, [1:getNumber(mt, 0);]);
    setConnectivity(mt, 1, 0, sort([getConnectivity(mt, 2, 0)[:,[1,2]];
                                    getConnectivity(mt, 2, 0)[:,[2,3]];
                                    getConnectivity(mt, 2, 0)[:,[1,3]]], 2));
    I, J = myUnique(getConnectivity(mt, 1, 0));
    setConnectivity(mt, 1, 0, I);
    setConnectivity(mt, 2, 1, J);
    setConnectivity(mt, 2, 1, reshape(getConnectivity(mt, 2, 1), getNumber(mt, D), 3));
    for d = 0:D
        setConnectivity(obj, d, d, [1:getNumber(mt, d);]);
    end
end

function getQuadRule(mt::MeshTopologyTri, order::Int, codim)
  if codim == 0
    return QuadTriGauss6(order)
  elseif codim == 1
    return QuadIntNodal(order)
  end
end

#-----------------------#
# Type MeshTopologyQuad #
#-----------------------#
type MeshTopologyQuad <: AbstractMeshTopology
    super :: MeshTopologyX

    function MeshTopologyQuad{Tn<:AbstractFloat, Te::Integer}
        (nodes::AbstractArray{Tn,2}, elems::AbstractArray{Te,2})
        mt = new(MeshTopologyX(nodes, elems))
        updateConnectivity(mt)
        return mt
    end
end

function updateConnectivity(mt::MeshTopologyQuad)
    D = getDim(mt)
    setConnectivity(mt, 0, 0, [1:getNumber(mt, 0;]);
    setConnectivity(mt, 1, 0, sort([getConnectivity(mt, 2, 0)[:,[1,2]];
                                    getConnectivity(mt, 2, 0)[:,[3,4]];
                                    getConnectivity(mt, 2, 0)[:,[1,3]];
                                    getConnectivity(mt, 2, 0)[:,[2,4]]], 2));
    I, J = myUnique(getConnectivity(mt, 1, 0));
    setConnectivity(mt, 1, 0, I);
    setConnectivity(mt, 2, 1, J);
    setConnectivity(mt, 2, 1, reshape(getConnectivity(mt, 2, 1), getNumber(mt, D), 4));
    for d = 0:D
        setConnectivity(mt, d, d, [1:getNumber(mt, d);]);
    end
end
