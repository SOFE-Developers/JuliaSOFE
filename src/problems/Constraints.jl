module Constraints

using ..Helpers
using ..Elements
using ..CMesh
using ..Meshes
using ..CMeshExtensions
using ..Spaces
using ..Operators

export AbstractConstraint, BoundaryCondition
export PeriodicBoundary, ZeroAverage
export transform, apply

abstract AbstractConstraint

abstract BoundaryCondition <: AbstractConstraint

#-----------------------#
# Type PeriodicBoundary #
#-----------------------#
type PeriodicBoundary <: BoundaryCondition
    fes :: FESpace
    master :: Function
    slave :: Function
end

# Associated Methods
# -------------------

"""

    transform(pbc::PeriodicBoundary)

Compute the transformation matrix associated
with the periodic boundary condition.
"""
function transform(pbc::PeriodicBoundary)
    sbmask = boundary(mesh(pbc.fes), pbc.slave)
    sdofind = dofIndices(pbc.fes, 1, sbmask)

    sdofnodes = evalReferenceMaps(mesh(pbc.fes),
                                  lagrangeNodesP(dimension(element(pbc.fes))-1,
                                                 order(element(pbc.fes))), 0)[sbmask,:,:]
    sdofnodes = unique(reshape(sdofnodes, (prod(size(sdofnodes, 1, 2)), size(sdofnodes, 3))), 1)

    mbmask = boundary(mesh(pbc.fes), pbc.master)
    mbfacets = facets(mesh(pbc.fes))[mbmask,:]

    ids = zeros(Int, size(sdofnodes, 1))
    dists = distance(sdofnodes, nodes(mesh(pbc.fes)), mbfacets; ids=ids)

    ns = normalize!(normals(nodes(mesh(pbc.fes)), mbfacets[ids,:]))
    mbfcells = cells(mesh(pbc.fes))[connectivity(topology(mesh(pbc.fes)), 1, 2)[mbmask,1],:]
    c = centroids(nodes(mesh(pbc.fes)), mbfcells[ids,:])
    p0 = nodes(mesh(pbc.fes))[mbfacets[ids,:],:][:,1,:]
    dot = sum(ns .* (c - p0), 2)

    msdofnodes = sdofnodes - copysign(dists, dot) .* ns
    hosts = connectivity(topology(mesh(pbc.fes)), 1, 2)[mbmask,1][ids,1]
    preimg = evalInverseMaps(mesh(pbc.fes), msdofnodes, hosts)

    mbasis = evalBasis(element(pbc.fes), preimg, 0)
    mdofs = dofMap(pbc.fes, dimension(element(pbc.fes)))[:,hosts]

    Is = repeat(sdofind, inner=(3,))
    Js = mdofs[:]
    Vs = mbasis[:]
    
    sdofmask = dofMask(pbc.fes, 1, sbmask)
    Ins = Jns = find(~sdofmask)
    Vns = ones(Ins)

    I = vcat(Ins, Is)
    J = vcat(Jns, Js)
    V = vcat(Vns, Vs)
    
    M = sparse(I, J, V, nDoF(pbc.fes), nDoF(pbc.fes))
    # M = M[:,Jns]

    return M
end

function apply(pbc::PeriodicBoundary, A::AbstractMatrix, b::AbstractVector)
    M = transform(pbc)

    A = M' * A * M
    b = M' * b

    sbmask = boundary(mesh(pbc.fes), pbc.slave)
    sdofind = dofIndices(pbc.fes, 1, sbmask)
    for sdof in sdofind
        A[sdof,sdof] = 1.
    end

    return A, b
end

#------------------#
# Type ZeroAverage #
#------------------#
type ZeroAverage <: AbstractConstraint
    fes :: FESpace
end

# Associated Methods
# -------------------

function apply(zac::ZeroAverage, A::AbstractMatrix, b::AbstractVector;
               idx::Integer=1)
    l = LinearForm(id, zac.fes, 1)
    assemble!(l)
    L = vector(l)
    
    A[idx,:] = L
    b[idx] = 0.

    return A, b
end



    
end # of module Constraints
