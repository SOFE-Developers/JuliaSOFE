module Spaces

using ..Elements
using ..Meshes

abstract AbstractFESpace

#--------------#
# Type FESpace #
#--------------#
type FESpace{M<:AbstractMesh, E<:AbstractElement} <: AbstractFESpace
    mesh :: M
    element :: E
    freeDOF :: Array{Bool, 1}
end

function FESpace{M<:AbstractMesh, E<:AbstractElement}(m::M, e::E, bfnc::Function=x->trues(size(x,1)))
    freeDOF = !getBoundaryDOFs(m, bfnc)
    return FESpace(m, e, freeDOF)
end

# Associated Methods
# -------------------
"""

    dofMap(fes::FESpace, d::Integer)

Returns the degrees of freedom mapping that connects the global mesh
entities of topological dimension `d` to the local reference element.

Establishes the connection between the local and global degrees of freedom
via the connectivity array `C` where `C[i,j] = k` connects the `i`-th
global basis function on the `j`-th element to the `k`-th local basis 
function on the reference element.
"""
function dofMap(fes::FESpace, d::Integer)
    dofTuple = fes.element.dofTUple
    dofPerDim = fes.element.dofPerDim
    nEntities = [getNumber(fes.mesh.topology, dd) for dd = 0:d]
    dofsNeeded = [nEntities[dd+1] * dofTuple[dd+1] for dd = 0:d]
    ndofs = [0, cumsum(dofsNeeded)...]
    dofs = [reshape(ndofs[i]+1:ndofs[i+1], dofTuple[i], nEntities[i]) for i = 1:d+1]

    dofMap = [zeros(Int, dofPerDim[i], nEntities[d+1]) for i = 1:d+1]
    # first iterate over subdims
    for dd = 0:d-1
        for (i, d_dd_i) in enumerate(getConnectivity(fes.mesh.topology, d, dd))
            println(d_dd_i)
            for (j, d_dd_ij) in enumerate(d_dd_i)
                r = (j-1)*dofTuple[dd+1]+1 : j*dofTuple[dd+1]
                dofMap[dd+1][r,i] = dofs[dd+1][:,j]
            end
        end
    end

    # set dofs of query dim
    dofMap[d+1] = dofs[d+1]

    return vcat(dofMap...)
end

"""

    nDoF(fes::FESpace)

Return the total number of degrees of freedom for the
finite element space `fes`.
"""
function nDoF(fes::FESpace)
    return maxabs(getDOFMap(fes, fes.mesh.dimension))
end

"""

    dofIndices(fes::FESpace, d::Integer)

Return the indices of the degrees of freedom for the finite
element space `fes` associated with the mesh entities of 
topological dimension `d`.
"""
function dofIndices(fes::FESpace, d::Integer)
    dofMap = dofMap(fes, d)
    return sort!(unique(dofMap))
end

"""

    dofMask(fes::FESpace, d::Integer)

Return a boolean mask specifying the degrees of freedom
for the finite element space `fes` associated with the 
mesh entities of topological dimension `d`.
"""
function dofMask(fes::FESpace, d::Integer)
    mask = zeros(Bool, nDoF(fes))
    for i in dofIndices(fes, d)
        mask[i] = true
    end
    return mask
end
    
end # of module Spaces
