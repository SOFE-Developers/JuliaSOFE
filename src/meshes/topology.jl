module Topology

export AbstractMeshTopology#, MeshTopologyTri, MeshTopologyQuad, MeshTopologyGeneric
export MeshTopology, getDim, getNodes, getConnectivity, getEntities, getNumber

#-----------------------------#
# Abstract Mesh Topology Type #
#-----------------------------#
abstract AbstractMeshTopology

# Associated Methods
# -------------------
"""

    getDim(mt::AbstractMeshTopology)

Return the spatial dimension of the mesh topology,
i.e. he topological dimension of the cells.
"""
function getDim(mt::AbstractMeshTopology)
    return mt.dimension
end

"""

    getNodes(mt::AbstractMeshTopology)

Return the coordinates of the mesh vertices.
"""
function getNodes(mt::AbstractMeshTopology)
    return mt.nodes
end

"""

    getConnectivity(mt::AbstractMeshTopology, d::Integer, dd::Integer)

  Return the incidence relation `d -> dd` computing it first if necessary.
"""
function getConnectivity(mt::AbstractMeshTopology, d::Integer, dd::Integer)
    return mt.connectivities[(d,dd)]
end

"""

    getEntities(mt::AbstractMeshTopology, d::Integer)

Return the vertex index connectivity array for the mesh
entities of topological dimension `d`.
"""
function getEntities(mt::AbstractMeshTopology, d::Integer)
    return getConnectivity(mt, d, 0)
end

"""

    getNumber(mt::AbstractMeshTopology, d::Integer)

Return the number of `d`-dimensional mesh entities.
"""
function getNumber(mt::AbstractMeshTopology, d::Integer)
    return size(getEntities(mt, d), 1)
end

"""

    getBoundary(mt::AbstractMeshTopology)

Determine the boundary facets.
"""
function getBoundary(mt::AbstractMeshTopology) end
function getBoundary(mt::AbstractMeshTopology, f::Function) end

#------------------------------#
# Concrete Mesh Topology Types #
#------------------------------#

USE_X = false

# Include Concrete Type Implementations
if USE_X
    include("topology_1.jl")
else
    include("topology_2.jl")
end

"""
Return a mesh topology instance of appropriate type. 
"""
function MeshTopology(nodes, cells)
    if USE_X
        if size(cells, 2) == 3
            # return MeshTopologyX(Tri, nodes, cells)
            return MeshTopologyTri(nodes, cells)
        elseif size(cells, 2) == 4
            # return MeshTopologyX(Quad, nodes, cells)
            return MeshTopologyQuad(nodes, cells)
        else
            error("Invalid cell connectivity!")
        end
    else
        return MeshTopologyGeneric(size(nodes, 2), nodes, cells)
    end
end


end # of module Topology
