module Topology

#export MeshTopologyX, MeshTopologyGeneric
#export getDim, getConnectivity

#-----------------------------#
# Abstract Mesh Topology Type #
#-----------------------------#
abstract AbstractMeshTopology

# Associated Methods
# -------------------

"""
Return the spatial dimension of the mesh topology,
i.e. he topological dimension of the cells.
"""
function getDim(mt::AbstractMeshTopology) end

"""

    getConnectivity(mt::AbstractMeshTopology, d::Integer, dd::Integer)

Return the incidence relation `d -> dd` computing it first if necessary.
"""
function getConnectivity(mt::AbstractMeshTopology, d::Integer, dd::Integer) end

"""

    getEntities(mt::AbstractMeshTopology, d::Integer)

Return the vertex index connectivity array for the mesh
entities of topological dimension `d`.
"""
function getEntities(mt::AbstractMeshTopology, d::Integer) end

"""

    getNumber(mt::AbstractMeshTopology, d::Integer)

Return the number of `d`-dimensional mesh entities.
"""
function getNumber(mt::AbstractMeshTopology, d::Integer) end

"""

    getBoundary(mt::AbstractMeshTopology)

Determine the boundary facets.
"""
function getBoundary(mt::AbstractMeshTopology) end
function getBoundary(mt::AbstractMeshTopology, f::Function) end

# Include Concrete Type Implementations
include("topology_1.jl")
include("topology_2.jl")

end # of module Topology
