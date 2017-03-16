export dofMap, nDoF, dofIndices, dofMask, extractDoFs, fixedDoF, freeDoF

"""

    dofMap(fes::FESpace, d::Integer, mask::AbstractVector)

  Return the degrees of freedom mapping that connects the global mesh
  entities of topological dimension `d` to the local reference element.

  Establishes the connection between the local and global degrees of freedom
  via a connectivity array `C` where `C[i,j] = k` connects the `i`-th
  local basis function  on the reference element to the `k`-th global basis 
  function on the `j`-th element.

  # Arguments
  * `d::Integer`: The topological dimension of the entities
                    for which to compute the dof map
  * `mask::Vector{T<:Integer}`: A mask marking specific entities
                    for which to compute the dof map
"""
function dofMap{T<:Integer}(fes::FESpace, d::Integer, mask::AbstractArray{T,1})
    dofs = generateDoFs(element(fes), mesh(fes), d)
    ndof = nDoF(element(fes), d)

    #M = zeros(T, sum(ndof), number(mesh(fes), d))
    M = zeros(Int, sum(ndof), number(mesh(fes), d))

    rb = 0
    for dd = 0:d-1
        d_dd = array(connectivity(topology(mesh(fes)), d, dd))
        ra = rb + 1; rb += ndof[dd+1]
        fill_dofMap!(view(M, ra:rb, :), dofs[dd+1], d_dd)
    end

    for j = 1:size(dofs[d+1], 2)
        for i = 1:size(dofs[d+1], 1)
            M[rb+i,j] = dofs[d+1][i,j]
        end
    end

    return M[:,mask]
end
dofMap(fes::FESpace, d::Integer) = dofMap(fes, d, 1:number(mesh(fes), d))
dofMap(fes::FESpace) = dofMap(fes, dimension(element(fes)))

function fill_dofMap!{T<:Integer}(M::AbstractArray{T,2}, dofs::AbstractArray{T,2}, inc::AbstractArray{T,2})
    ptr = 1
    for i = 1:size(inc, 1)
        for j = 1:size(inc, 2)
            for k = 1:size(dofs, 1)
                M[ptr,i] = dofs[k, inc[i,j]]
                ptr += 1
            end
        end
        ptr = 1
    end
    return nothing
end

function generateDoFs{T<:AbstractElement}(elem::T, mesh::Mesh, d::Integer)
    doftuple = dofTuple(elem)[1:d+1]
    nentities = [number(mesh, dd) for dd = 0:d]
    ndofs = map(*, doftuple, nentities)
    dofrange = [1, 1 + cumsum(ndofs)...]
    dofs = [reshape(dofrange[i]:dofrange[i+1]-1, doftuple[i], nentities[i])
            for i = 1:d+1]

    return dofs
end

function assembleDoFs{T<:Integer}(dofs::AbstractArray{T,1}, d::Integer)
end

"""

    nDoF(fes::FESpace)

  Return the total number of degrees of freedom for the
  finite element space `fes`.
"""
function nDoF(fes::FESpace)
    return maxabs(dofMap(fes, dimension(mesh(fes))))
end

"""

    dofIndices(fes::FESpace, args...)

  Return the indices of the degrees of freedom for the finite
  element space `fes` associated with the mesh entities of 
  topological dimension `d`.

  # Arguments
  * `d::Integer`: The topological dimension of the entities
                    for which to compute the dof indices
  * `mask::Vector{T<:Integer}`: A mask marking specific entities
                    for which to compute the dof indices
"""
function dofIndices(fes::FESpace, args...)
    M = dofMap(fes, args...)
    return sort!(unique(M))
end

"""

    dofMask(fes::FESpace, args...)

  Return a boolean mask specifying the degrees of freedom
  for the finite element space `fes` associated with the 
  mesh entities of topological dimension `d`.

  # Arguments
  * `d::Integer`: The topological dimension of the entities
                    for which to compute the dof mask
  * `mask::Vector{T<:Integer}`: A mask marking specific entities
                    for which to compute the dof mask
"""
function dofMask(fes::FESpace, args...)
    mask = zeros(Bool, nDoF(fes))
    for i in dofIndices(fes, args...)
        mask[i] = true
    end
    return mask
end

extractDoF(fes::FESpace, args...) = dofMask(fes, args...)

"""

    fixedDoF(fes::FESpace)

  Return a boolean mask marking the constrained
  degrees of freedom of the finite element space.
"""
function fixedDoF(fes::FESpace)
    dim = dimension(mesh(fes)) - 1
    bmask = boundary(mesh(fes), domain(fes))
    return extractDoF(fes, dim, bmask)
end

"""

    freeDoF(fes::FESpace)

  Return a boolean mask marking the unconstrained
  degrees of freedom of the finite element space.
"""
freeDoF(fes::FESpace) = !fixedDoF(fes)

