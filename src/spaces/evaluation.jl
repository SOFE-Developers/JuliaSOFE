import ..CMeshExtensions: evaluate

export evaluate

"""

    evaluate{T<:AbstractFloat}(fes::FESpace, dofs::AbstractVector{T}, 
                               points::AbstractArray{T,2}, deriv::Integer=0)

  Evaluate the linear combination of the finite element
  space's basis functions or their derivatives in the 
  given local `points` w.r.t. the given `dof` values.
"""
function evaluate{T<:Float,S<:Float}(fes::FESpace, dofs::AbstractVector{S},
                                     points::AbstractArray{T,2}, deriv::Integer=0)
    nP, dimP = size(points)
    dofmap = dofMap(fes, dimP)
    nB, nE = size(dofmap)
    basis = evalBasis(element(fes), points, deriv)
    nB, nP, nC, nD = size(basis, 1:4...)

    dofs = convert(Vector{promote_type(T, S)}, dofs)
    
    if deriv == 0
        @assert nD == 1
        U = zeros(T,nE,nP,nC)
    elseif deriv == 1
        U = zeros(T,nE,nP,nC,nD)
    elseif deriv == 2
        U = zeros(T,nE,nP,nC,nD,nD)
    end

    fill_dofvalues!(U, basis, dofs, dofmap)

    return U
end

function fill_dofvalues!{T<:AbstractFloat}(U::AbstractArray{T,3}, basis::Array{T,3},
                                   dofs::AbstractVector{T}, dofmap::Array{Int,2})
    nE, nP, nC = size(U)
    nB, nP, nC = size(basis)
    nB, nE = size(dofmap)
    
    for ic = 1:nC
        for ip = 1:nP
            for ie = 1:nE
                for ib = 1:nB
                    U[ie,ip,ic] += dofs[dofmap[ib,ie]] * basis[ib,ip,ic]
                end
            end
        end
    end

    return nothing
end

function fill_dofvalues!{T<:AbstractFloat}(U::Array{T,4}, basis::Array{T,4},
                                   dofs::AbstractVector{T}, dofmap::Array{Int,2})
    nE, nP, nC, nD = size(U)
    nB, nP, nC, nD = size(basis)
    nB, nE = size(dofmap)
    
    for id = 1:nD
        fill_dofvalues!(view(U, :, :, :, id), basis[:,:,:,id], dofs, dofmap)
    end

    return nothing
end

function fill_dofvalues!{T<:AbstractFloat}(U::Array{T,5}, basis::Array{T,5},
                                   dofs::AbstractVector{T}, dofmap::Array{Int,2})
    nE, nP, nC, nDi, nDj = size(U)
    nB, nP, nC, nDi, nDj = size(basis)
    nB, nE = size(dofmap)

    for jd = 1:nDj
        for id = 1:nDi
            fill_dofvalues!(view(U, :, :, :, id, jd), basis[:,:,:,id,jd], dofs, dofmap)
        end
    end

    return nothing
end

