module Topology

export MeshTopologyX, MeshTopologyGeneric
export getDim, getConnectivity

abstract AbstractMeshTopology

include("topology_1.jl")
include("topology_2.jl")

end # of module Topology
