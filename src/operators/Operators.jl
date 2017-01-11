__precompile__()

module Operators

using CMesh.Meshes
using ..CMeshExtensions
using ..Elements
using ..Spaces
using ..Quadrature
using ..Helpers

import ..CMeshExtensions: evaluate

export AbstractOperator
export space, coeff, assemble, assemble!

typealias Float AbstractFloat

include("coefficients.jl")
include("basic.jl")

#------------------------#
# Variational Form Types #
#------------------------#
abstract AbstractVariationalOperator

include("bilinear.jl")
include("linear.jl")

#---------------------#
# Assembling Routines #
#---------------------#
assemble(op::AbstractVariationalOperator) = assemble(op, dimension(mesh(testspace(op))))
assemble!(op::AbstractVariationalOperator) = assemble!(op, dimension(mesh(testspace(op))))
assemble!(a::BilinearForm, d::Integer) = matrix!(a, assemble(a, d))
assemble!(l::LinearForm, d::Integer) = vector!(l, assemble(l, d))

function assemble(a::BilinearForm, d::Integer)
    qpoints, qweights = quadData(a.quadrule, d)
    
    fes = trialspace(a)
    dMap = dofMap(fes, d)
    ndof = nDoF(fes)
    
    nB, nE = size(dMap)
    nP, nD = size(qpoints)

    dofI = zeros(Int, nE, nB, nB)
    dofJ = zeros(Int, nE, nB, nB)

    for jb = 1:nB
        for ib = 1:nB
            for ie = 1:nE
                dofI[ie,ib,jb] = dMap[ib,ie]
                dofJ[ie,ib,jb] = dMap[jb,ie]
            end
        end
    end
    
    C, U, V = evaluate(a, d)
    Det = abs(evalJacobianDeterminat(mesh(fes), qpoints))
    
    entries = zeros(eltype(qpoints), nE, nB, nB)
    # fill_entries!(a, entries, C, U, V, qweights, Det)
    fill_entries!(entries, C, U, V, qweights, Det)

    A = sparse(dofI[:], dofJ[:], entries[:], ndof, ndof)
    return A
end

function assemble(l::LinearForm, d::Integer)
    qpoints, qweights = quadData(l.quadrule, d)

    fes = testspace(l)
    dMap = dofMap(fes, d)
    ndof = nDoF(fes)

    nB, nE = size(dMap)
    nP, nD = size(qpoints)

    dofI = zeros(Int, nE, nB)
    
    for ib = 1:nB
        for ie = 1:nE
            dofI[ie,ib] = dMap[ib,ie]
        end
    end
    
    C, V = evaluate(l, d)
    Det = abs(evalJacobianDeterminat(mesh(fes), qpoints))
    
    entries = zeros(eltype(qpoints), nE, nB)
    # fill_entries!(l, entries, C, V, qweights, Det)
    fill_entries!(entries, C, V, qweights, Det)
    
    L = sparsevec(dofI[:], entries[:], ndof)
    return full(L)
end

include("entries.jl")

end # of module Operators

