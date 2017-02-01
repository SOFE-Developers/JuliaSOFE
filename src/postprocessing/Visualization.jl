module Visualization

using PyCall

using CMesh
import CMesh: plot

using ..Elements
using ..CMeshExtensions
using ..Spaces

export plot

const pyvisofe = PyNULL()

function __init__()
    pygui_start(:qt)
    copy!(pyvisofe, pyimport("pyvisofe"))
end

visgrid{T<:PElement}(el::Element{T}, resolution::Integer) =
    UnitSimplex(Elements.dimension(el), resolution)
visgrid{T<:QElement}(el::Element{T}, resolution::Integer) =
    UnitOrthotope(Elements.dimension(el), resolution)

visgrid(m::Mesh, resolution::Integer) = visgrid(m, visgrid(m.element, resolution))
function visgrid(m::Mesh, refgrid::Mesh)
    refnodes = nodes(refgrid)
    refcells = entities(refgrid, dimension(refgrid))
    nE_ref, nV = size(refcells)
    nN_ref, nW = size(refnodes)

    visnodes = evalReferenceMaps(m, refnodes, 0)
    nE_vis, nN_ref, nW = size(visnodes)

    viscells = broadcast(+, reshape(refcells, (nE_ref,1,nV)),
                         reshape(nN_ref * collect(0:nE_vis-1), (1,nE_vis,1)))

    visnodes = permutedims(visnodes, [2,1,3])
    visnodes = reshape(visnodes, (nE_vis*nN_ref, nW))
    viscells = reshape(viscells, (nE_ref*nE_vis,nV))

    return Mesh(visnodes, viscells)
end

function plot(el::Element, i::Integer; resolution::Integer=100,
              size::Integer=2, alpha::AbstractFloat=1.)
    D = Elements.dimension(el)
    vismesh = visgrid(el, resolution)
    n = nodes(vismesh)
    basis = evalBasis(el, n, 0)

    if D == 1
        x = n[:,1]; y = basis[i,:]

        pyvisofe[:plot](x, y)
    elseif D == 2
        x = n[:,1]; y = n[:,2]; z = basis[i,:]
        faces = entities(vismesh, 2) - 1

        pyvisofe[:trisurface](x, y, z, faces)
    elseif D == 3
        x = n[:,1]; y = n[:,2]; z = n[:,3]
        c = basis[i,:]
        pyvisofe[:scatter3d](x, y, z, c=c,
                             size=size, alpha=alpha)
    end
end

function plot{T<:AbstractFloat}(fes::FESpace, uh::AbstractVector{T};
                                resolution::Integer=3, edge_color="black")
    refmesh = visgrid(element(fes), resolution)
    vismesh = visgrid(mesh(fes), refmesh)

    values = evaluate(fes, uh, nodes(refmesh), 0)
    nE, nP, nC = size(values)
    @assert nC == 1
    values = permutedims(values, [2,1,3])

    D = dimension(vismesh)
    coords = nodes(vismesh)
    if D == 1
        x = coords[:,1]; y = values[:]
        
        pyvisofe[:plot](x, y)
    elseif D == 2
        x = coords[:,1]; y = coords[:,2]; z = values[:]
        faces = entities(vismesh, 2) - 1
        edges = entities(vismesh, 1)[boundary(vismesh,1),:] - 1

        pyvisofe[:trisurface](x, y, z, faces,
                              edges=edges, edge_color=edge_color)
    elseif D == 3
        x = coords[:,1]; y = coords[:,2]; z = coords[:,3]
        c = values[:]

        pyvisofe[:scatter3d](x, y, z,
                             c=c, size=2., alpha=0.5)
    end
end

end # of module Visualization
