using CMesh
using PyCall

pygui_start(:qt)
@pyimport pyvisofe
@pyimport scipy.ndimage as scipy_ndimage
convolve = scipy_ndimage.convolve
        
#--------------------------#
# Random Coefficient Field #
#--------------------------#
function checkerboard(ncells::Integer=2)
    n2 = Int(ceil(ncells/2))
    X = repeat([true  false;
                false true], outer=(n2, n2))
    if isodd(ncells)
        X = X[1:end-1,1:end-1]
    end

    return X
end

function random_checkerboard(ncells::Integer=2, conv::Bool=false;
                             beta::Real=1.)
    X = rand(Bool, ncells, ncells)

    if conv
        # Convolution
        kernel(x::Matrix; b::Real=1.) = b/π * exp(-b * sumabs2(x, 2))
        x = 1
        while kernel([x]'')[1] > 1e-8 && x < 1e4
            x += 1
        end
        nx = 2*x+1
        ls = linspace(-x, x, nx)
        grid = tensorprod(ls, ls)
        K = reshape(kernel(grid, b=beta), (nx, nx))
        
        C = convolve(float(X), K; mode="wrap")
        X = C .> mean(C)
    end

    return X
end

function show_checkerboard(epsilon::Real=(1/2);
                           rnd::Bool=false, cnv::Bool=false, beta::Real=1.)
    if rnd
        if cnv
            X = random_checkerboard(Int(1/epsilon), true, beta=beta)
        else
            X = random_checkerboard(Int(1/epsilon), false)
        end
    else
        X = checkerboard(Int(1/epsilon))
    end

    show_checkerboard(X)
end
function show_checkerboard(X::Union{Matrix{Bool}, BitArray{2}})
    epsilon = 1/sqrt(length(X))
    
    function inside{T<:Real}(xs::AbstractVector{T})
        cxs = Vector{Int}(ceil(xs/epsilon))
        #return all(iseven.(c)) | all(isodd.(c))
        return X[cxs[1], cxs[2]]
    end

    m = UnitSquare(Int(1/epsilon)+1)
    
    n = nodes(m)
    c = cells(m)

    cc = centroids(n , c)
    mask = [inside(cc[i,:]) for i = 1:size(cc, 1)]
    
    C = ones(size(cc, 1), 3)
    C[mask,:] = 0.
    
    pyvisofe.trisurface(n[:,1], n[:,2], zeros(size(n, 1)), c - 1,
                        face_colors=C)

    return X
end
