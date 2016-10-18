# Helper Function
function myUnique(a) # TODO
  ua = unique(a,1)
  N = size(a,1)
  M = size(ua,1)
  hasha = [hash(a[k,:]) for k = 1:size(a,1)]
  hashua = [hash(ua[k,:]) for k = 1:size(ua,1)]
  ind = zeros(Int,N)
  for n = 1:N
    for m = 1:M
      if hasha[n] == hashua[m]
         ind[n] = m
         break
      end
    end
  end
  return (ua, ind)
end


abstract AbstractMeshTopologyX <: AbstractMeshTopology

#--------------------#
# Type MeshTopologyX #
#--------------------#
type MeshTopologyX <: AbstractMeshTopologyX
    dimension :: Integer
    nodes :: Array{AbstractFloat, 2}
    connectivities :: Dict{Tuple{Int, Int}, Array{Int,2}}

    function MeshTopologyX{Tn<:AbstractFloat, Te<:Integer}(nodes::AbstractArray{Tn,2}, elems::AbstractArray{Te,2})
        dim = size(nodes, 2)
        connect = Dict((dim,0) => elems)
        return new(dim, nodes, connect)
    end
end

# Associated Methods
# -------------------
function getConnectivity(mt::AbstractMeshTopologyX, d::Integer, dd::Integer)
    return mt.super.connectivities[(d,dd)]
end

function setConnectivity(mt::AbstractMeshTopologyX, d::Integer, dd::Integer, connect::Array{Int, 2})
    mt.super.connectivities[(d,dd)] = connect;
end

function getEntities(mt::AbstractMeshTopologyX, d::Integer)
    if d == 0
        return mt.super.nodes
    else
        return mt.super.connectivities[(d,0)]
    end
end

function getNumber(mt::AbstractMeshTopologyX, d::Integer)
    return size(getEntities(mt, d), 1)
end

function getBoundary(mt::AbstractMeshTopologyX, f::Function=x->true)
    e2F = mt.super.connectivities[(2,1)][:]
    R = (sparse(e2F, Array{Int64}(ones(length(e2F))), Array{Int64}(ones(length(e2F)))).==1)[:];
    face = getEntities(mt, mt.dimension-1)
    center = permutedims(sum(mt.super.nodes[face,:],2)/size(face,2), [1 3 2])
    R = R & f(center)
    return R
end


getNodes(mt::AbstractMeshTopologyX) = mt.super.nodes
getDim(mt::AbstractMeshTopologyX) = mt.super.dimension

#----------------------#
# Type MeshTopologyTri #
#----------------------#
type MeshTopologyTri <: AbstractMeshTopologyX
    super :: MeshTopologyX

    function MeshTopologyTri{Tn<:AbstractFloat, Te<:Integer}(nodes::AbstractArray{Tn,2}, elems::AbstractArray{Te,2})
        mt = new(MeshTopologyX(nodes, elems))
        #updateConnectivity(mt)
        return mt
    end
end

# Associated Methods
# -------------------
function updateConnectivity(mt::MeshTopologyTri)
    D = getDim(mt)
    setConnectivity(mt, 0, 0, collect(1:getNumber(mt, 0))'');
    setConnectivity(mt, 1, 0, sort([getConnectivity(mt, 2, 0)[:,[1,2]];
                                    getConnectivity(mt, 2, 0)[:,[2,3]];
                                    getConnectivity(mt, 2, 0)[:,[1,3]]], 2));
    I, J = myUnique(getConnectivity(mt, 1, 0));
    setConnectivity(mt, 1, 0, I);
    setConnectivity(mt, 2, 1, reshape(J, getNumber(mt, 2), 3));
    for d = 0:D
        setConnectivity(mt, d, d, [1:getNumber(mt, d);]);
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
type MeshTopologyQuad <: AbstractMeshTopologyX
    super :: MeshTopologyX

    function MeshTopologyQuad{Tn<:AbstractFloat, Te<:Integer}(nodes::AbstractArray{Tn,2}, elems::AbstractArray{Te,2})
        mt = new(MeshTopologyX(nodes, elems))
        updateConnectivity(mt)
        return mt
    end
end

function updateConnectivity(mt::MeshTopologyQuad)
    D = getDim(mt)
    setConnectivity(mt, 0, 0, [1:getNumber(mt, 0)]);
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
