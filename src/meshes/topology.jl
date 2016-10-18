module Topology

export AbstractMeshTopology, MeshTopologyTri, MeshTopologyQuad, MeshTopologyGeneric
export MeshTopology, getDim, getNodes, getConnectivity, getEntities, getNumber

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

#------------------------------#
# Concrete Mesh Topology Types #
#------------------------------#

USE_GENERIC = false

# Include Concrete Type Implementations
if !USE_GENERIC
    include("topology_1.jl")
else
    include("topology_2.jl")
end

"""
Return a mesh topology instance of appropriate type. 
"""
function MeshTopology(nodes, cells)
    if !USE_GENERIC
        if size(cells, 2) == 3
            return MeshTopologyTri(nodes, cells)
        elseif size(cells, 2) == 4
            return MeshTopologyQuad(nodes, cells)
        else
            error("Invalid cell connectivity!")
        end
    else
        return MeshTopologyGeneric(size(nodes, 2), nodes, cells)
    end
end


end # of module Topology
