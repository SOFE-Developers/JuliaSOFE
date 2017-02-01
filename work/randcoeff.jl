using CMesh
using PyCall

pygui_start(:qt)
@pyimport pyvisofe

#--------------------------#
# Random Coefficient Field #
#--------------------------#
epsilon = 1/2^5

X = rand(Bool, Int(1/epsilon), Int(1/epsilon))

# Convolution
@pyimport scipy.ndimage as scipy
convolve = scipy.convolve

kernel(x::Matrix; b::Real=1.) = b/π * sumabs2(x, 2)


function inside{T<:Real}(xs::AbstractVector{T})
    c = Vector{Int}(ceil(xs/epsilon))
    #return all(iseven.(c)) | all(isodd.(c))
    return X[c[1], c[2]]
end

#------#
# Mesh #
#------#
m = UnitSquare(Int(1/epsilon)+1)

n = nodes(m)
c = cells(m)

cc = centroids(n , c)
mask = [inside(cc[i,:]) for i = 1:size(cc, 1)]

C = ones(size(c, 1), 3)
C[mask,:] = 0.

pyvisofe.trisurface(n[:,1], n[:,2], zeros(size(n, 1)), c - 1,
                    face_colors=C)
