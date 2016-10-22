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

typealias Float AbstractFloat

#--------------------#
# Type MeshTopologyX #
#--------------------#
type MeshTopologyX{X, Tn<:Float, Te<:Integer} <: AbstractMeshTopology
    dimension :: Int
    nodes :: Array{Tn, 2}
    connectivities :: Dict{Tuple{Int, Int}, Array{Te,2}}
end

function MeshTopologyX{X<:AbstractMeshTopology, Tn<:Float, Te<:Integer}(::Type{X},
                                                                        nodes::AbstractArray{Tn,2},
                                                                        elems::AbstractArray{Te,2})
    dim = size(nodes, 2)
    connect = Dict((dim,0) => elems)
    mt = MeshTopologyX{X,Tn,Te}(dim, nodes, connect)
    updateConnectivity!(mt)
    return mt
end

# Associated Methods
# -------------------
function setConnectivity(mt::AbstractMeshTopology, d::Integer, dd::Integer, connect::AbstractArray{Int, 2})
    mt.connectivities[(d,dd)] = connect;
end

function getBoundary{X<:AbstractMeshTopology}(mt::MeshTopologyX{X}, f::Function=x->true)
    e2F = mt.connectivities[(2,1)][:]
    R = (sparse(e2F, Array{Int64}(ones(length(e2F))), Array{Int64}(ones(length(e2F)))).==1)[:];
    face = getEntities(mt, mt.dimension-1)
    center = permutedims(sum(mt.nodes[face,:],2)/size(face,2), [1 3 2])
    R = R & f(center)
    return R
end

#------------------------#
# Type MeshTopology{Tri} #
#------------------------#
type Tri <: AbstractMeshTopology
end

typealias MeshTopologyTri MeshTopologyX{Tri}
MeshTopologyTri(nodes, cells) = MeshTopologyX(Tri, nodes, cells)

# Associated Methods
# -------------------
function updateConnectivity!(mt::MeshTopologyX{Tri})
    D = getDim(mt)
    setConnectivity(mt, 0, 0, collect(1:size(getNodes(mt),1))'');
    setConnectivity(mt, 1, 0, sort([getConnectivity(mt, 2, 0)[:,[1,2]];
                                    getConnectivity(mt, 2, 0)[:,[2,3]];
                                    getConnectivity(mt, 2, 0)[:,[1,3]]], 2));
    I, J = myUnique(getConnectivity(mt, 1, 0));
    setConnectivity(mt, 1, 0, I);
    setConnectivity(mt, 2, 1, reshape(J, getNumber(mt, 2), 3));
    for d = 0:D
        nE = getNumber(mt, d)
        setConnectivity(mt, d, d, reshape(1:nE, nE, 1));
    end
end

function getQuadRule(mt::MeshTopologyX{Tri}, order::Int, codim)
  if codim == 0
    return QuadTriGauss6(order)
  elseif codim == 1
    return QuadIntNodal(order)
  end
end

#-------------------------#
# Type MeshTopology{Quad} #
#-------------------------#
type Quad <: AbstractMeshTopology
end

typealias MeshTopologyQuad MeshTopologyX{Quad}
MeshTopologyQuad(nodes, cells) = MeshTopologyX(Quad, nodes, cells)

function updateConnectivity!(mt::MeshTopologyX{Quad})
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
