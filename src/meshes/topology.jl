module Topology

export AbstractMeshTopology
export MeshTopology
export dimension, nodes, connectivity, entities
export getDim, getNodes, getConnectivity, getEntities, getNumber

import ...Helpers: dimension

#-----------------------------#
# Abstract Mesh Topology Type #
#-----------------------------#
abstract AbstractMeshTopology

# Associated Methods
# -------------------
"""

    dimension(mt::AbstractMeshTopology)

  Return the spatial dimension of the mesh topology,
  i.e. the topological dimension of the cells.
"""
@inline dimension(mt::AbstractMeshTopology) = getfield(mt, :dimension)
@inline getDim(mt::AbstractMeshTopology) = dimension(mt)

"""

    nodes(mt::AbstractMeshTopology)

  Return the coordinates of the mesh vertices.
"""
@inline nodes(mt::AbstractMeshTopology) = getfield(mt, :nodes)
@inline getNodes(mt::AbstractMeshTopology) = nodes(mt)

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
    nV = size(cells, 2)
    
    if USE_X
        if nV == 3
            # return MeshTopologyX(Tri, nodes, cells)
            return MeshTopologyTri(nodes, cells)
        elseif nV == 4
            # return MeshTopologyX(Quad, nodes, cells)
            return MeshTopologyQuad(nodes, cells)
        else
            error("Invalid cell connectivity!")
        end
    else
        if size(nodes, 2) == 1
            nV == 2 && return MeshTopology(Simp, nodes, cells)
        elseif size(nodes, 2) == 2
            nV == 3 && return MeshTopology(Simp, nodes, cells)
            nV == 4 && return MeshTopology(Orth, nodes, cells)
        elseif size(nodes, 2) == 3
            nV == 4 && return MeshTopology(Simp, nodes, cells)
            nV == 8 && return MeshTopology(Orth, nodes, cells)
        end

        error("Invalid case...")
    end
end


end # of module Topology
