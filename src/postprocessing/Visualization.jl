module Visualization

using PyCall

using ..Elements
using ..Meshes
using ..Spaces

export show

const pyvisofe = PyNULL()

function __init__()
    pygui_start(:qt)
    copy!(pyvisofe, pyimport("pyvisofe"))
end

visgrid{T<:PElement}(el::Element{T}, resolution::Integer) =
    UnitSimplex(Elements.dimension(el), resolution)
visgrid{T<:QElement}(el::Element{T}, resolution::Integer) =
    UnitOrthotope(Elements.dimension(el), resolution)

function visgrid(m::Mesh, resolution::Integer)
    refmesh = visgrid(m.element, resolution)
    refnodes = nodes(refmesh)
    refcells = entities(refmesh, dimension(refmesh))
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

function show(m::Mesh)
    D = Meshes.dimension(m)
    
    if D == 1
        x = nodes(m)[:,1]; y = zeros(x)

        pyvisofe[:plot](x, y)
    elseif D == 2
        n = nodes(m)
        x = n[:,1]; y = n[:,2]
        faces = entities(m, 2) - 1
        
        pyvisofe[:triplot](x, y, faces)
    elseif D == 3
        n = nodes(m)
        x = n[:,1]; y = n[:,2]; z = n[:,3]
        faces = entities(m, 2) - 1

        pyvisofe[:wireframe](x, y, z, faces)
    else
        error("Invalid mesh dimension! ($D)")
    end
end

function show(el::Element, i::Integer; resolution::Integer=100,
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
        faces = entities(m, 2) - 1

        pyvisofe[:trisurface](x, y, z, faces)
    elseif D == 3
        x = n[:,1]; y = n[:,2]; z = n[:,3]
        c = basis[i,:]
        pyvisofe[:scatter3d](x, y, z, c=c,
                             size=size, alpha=alpha)
    end
end

function show{T<:AbstractFloat}(fes::FESpace, uh::AbstractVector{T};
                                resolution::Integer=3)
    vismesh_loc = visgrid(element(fes), resolution)
    
end

end # of module Visualization
