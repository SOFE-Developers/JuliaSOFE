module Visualization

using PyCall

using ..Elements
using ..Meshes

const pyvisofe = PyNULL()

function __init__()
    pygui_start(:qt)
    copy!(pyvisofe, pyimport("pyvisofe"))
end

visgrid{T<:PElement}(el::Element{T}, resolution::Integer) =
    UnitSimplex(Elements.dimension(el), resolution)
visgrid{T<:QElement}(el::Element{T}, resolution::Integer) =
    UnitOrthotope(Elements.dimension(el), resolution)

function showMesh(m::Mesh)
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

function showElement(el::Element, i::Integer, resolution::Integer=100)
    D = Elements.dimension(el)
    m = visgrid(el, resolution)
    n = nodes(m)
    basis = evalBasis(el, n)

    if D == 1
        x = n[:,1]; y = basis[i,:]

        pyvisofe[:plot](x, y)
    elseif D == 2
        x = n[:,1]; y = n[:,2]; z = basis[i,:]
        faces = entities(m, 2) - 1

        pyvisofe[:trisurface](x, y, z, faces)
    elseif D == 3
        error("Not Implemented!")
    end
end
    
end # of module Visualization
