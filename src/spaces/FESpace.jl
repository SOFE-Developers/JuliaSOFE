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
function getDOFMap(fes::FESpace, d::Integer)
    dofTuple = fes.element.dofTUple
    dofPerDim = fes.element.dofPerDim
    nEntities = [getNumber(fes.mesh.topology, dd) for dd = 0:d]
    dofsNeeded = [nEntities[dd+1] * dofTuple[dd+1]
                  for dd = 0:d]
    ndofs = [0, cumsum(dofsNeeded)...]
    dofs = [reshape(ndofs[i]+1:ndofs[i+1], dofTuple[i], nEntities[i]) for i = 1:d+1]

    dofMap = [zeros(Int, dofPerDim[i], nEntities[d]) for i = 1:d+1]
    # first iterate over subdims
    for dd = 0:d-1
        for (i, d_dd_i) in enumerate(getConnectivity(fes.mesh.topology, d, dd))
            for (j, d_dd_ij) in enumerate(d_dd_i)
                dofMap[dd][:,i] = ...
            end
        end
    end

    
end

end # of module Spaces
