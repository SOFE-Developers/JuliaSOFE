module Visualization

using PyCall

using ..Meshes

const pyvisofe = PyNULL()

function __init__()
    pygui_start(:qt)
    copy!(pyvisofe, pyimport("pyvisofe"))
end

function showMesh(m::Mesh)
    D = dimension(m)
    
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

end # of module Visualization
