using ..Helpers

export TensorProductMesh
export UnitOrthotope, UnitInterval, UnitSquare, UnitCube
export UnitSimplex, UnitTriangle

"""

    TensorProductMesh{T<:AbstractFloat}(grids::AbstractArray{T,1}...)

  Create a mesh where the nodes are generated from the tensor product
  of the given grid points.
"""
function TensorProductMesh{T<:AbstractFloat}(grids::AbstractArray{T,1}...)
    nodes = tensorprod(grids...)
    nnodes, dim = size(nodes)

    if dim == 1
        cells = transpose(reshape(repeat(1:nnodes, inner=2)[2:end-1], (2,nnodes-1)))
    elseif dim == 2
        m = length(grids[1])
        n = length(grids[2])

        rmask = ones(Bool, m*(n-1)-1)
        rmask[m:m:m*(n-1)-1] = false

        I = find(rmask)
        C = hcat(I, I+1, I+m, I+m+1)

        if false
            # with copy
            cells = vcat(C[:,[1,2,3]],
                         C[:,[4,3,2]])
        else
            # w/o copy
            cells = reshape(view(C, :, [1 2 3;
                                        4 3 2]),
                            (2*length(I),3))
        end
    elseif dim == 3
        m = length(grids[1])
        n = length(grids[2])
        p = length(grids[3])

        rmask = ones(Bool, m*n*(p-1)-m-1)
        r1 = m:m:m*n*(p-1)-m-1
        r2 = broadcast(+, m*(n-1)+1:m*n-1, collect((0:p-3)*m*n)')
        rmask[vcat(r1, r2[:])] = false

        I = find(rmask)
        C = hcat(I, I+1)
        C = hcat(C, C+m)
        C = hcat(C, C+m*n)

        if false
            # with copy
            # cells = vcat(C[:,[1,2,3,5]],
            #              C[:,[2,3,4,6]],
            #              C[:,[2,3,5,6]],
            #              C[:,[3,4,6,7]],
            #              C[:,[3,5,6,7]],
            #              C[:,[4,6,7,8]])
            cells = vcat(C[:,[1,2,3,5]],
                         C[:,[4,3,2,8]],
                         C[:,[7,5,8,3]],
                         C[:,[6,8,5,2]],
                         C[:,[2,3,5,8]])
        else
            #w/o copy
            # cells = reshape(view(C, :, [1 2 3 5;
            #                             2 3 4 6;
            #                             2 3 5 6;
            #                             3 4 6 7;
            #                             3 5 6 7;
            #                             4 6 7 8]),
            #                 (6*length(I), 4))
            cells = reshape(view(C, :, [1 2 3 5;
                                        4 3 2 8;
                                        7 5 8 3;
                                        6 8 5 2;
                                        2 3 5 8]),
                            (5*length(I), 4))
        end
    else
            error("Mesh dimension $dim not available!")
    end
    
    return Mesh(nodes, cells)
end

"""

    UnitInterval(n::Integer)

  Create a mesh discretizing the unit interval [0,1] using
  `n` nodes on the x-axis.
"""
UnitInterval(n::Integer) = TensorProductMesh(linspace(0,1,n))

#------------------#
# Orthotope Meshes #
#------------------#
# -> Hyperrectangles

"""

    UnitOrthotope(d::Integer, ns::Integer...)

  Create a mesh discretizing the unit orthotope [0,1]^`d` 
  using the `d` given number of nodes `ns` for each axis.
"""
function UnitOrthotope(d::Integer, ns::AbstractVector)
    @assert d == length(ns)
    return TensorProductMesh([linspace(0,1,n) for n in ns]...)
end
UnitOrthotope(d::Integer, n::Integer) = UnitOrthotope(d, repeat([n], inner=d))

"""

    UnitSquare(nx::Integer, ny::Integer)

  Create a mesh discretizing the unit square [0,1]² using
  `nx` nodes on the x-axis and `ny` nodes on the y-axis.
"""
UnitSquare(nx::Integer, ny::Integer) = UnitOrthotope(2, [nx, ny])
UnitSquare(n::Integer) = UnitSquare(n,n)

"""

    UnitCube(nx::Integer, ny::Integer, nz::Integer)

  Create a mesh discretizing the unit cube [0,1]³ using
  `nx` nodes on the x-axis, `ny` nodes on the y-axis
  and `nz` nodes on the z-axis.
"""
UnitCube(nx::Integer, ny::Integer, nz::Integer) = UnitOrthotope(3, [nx, ny, nz])
UnitCube(n::Integer) = UnitCube(n,n,n)

#----------------#
# Simplex Meshes #
#----------------#

"""

    UnitSimplex(d::Integer, n::Integer)

  Create a mesh discretizing the unit simplex in R^`d`
  using the `n` nodes for each axis.
"""
function UnitSimplex(d::Integer, n::Integer)
    m = UnitOrthotope(d, n)
    n = nodes(m)
    #c = entities(m, dimension(m))
    c = entities(m, d)

    mask = sum(n,2)[:] .< 1+1e-9
    off = cumsum(!mask)

    c = ndarray([c[i,:] for i = 1:size(c,1) if sum(mean(n[c[i,:],:],1)) < 1.])
    for j = 1:size(c,2)
        for i = 1:size(c,1)
            c[i,j] -= off[c[i,j]]
        end
    end

    n = n[mask,:]

    return Mesh(n,c)
end

"""

    UnitTriangle(n::Integer)

  Create a mesh discretizing the unit triangle 
  defined by the nodes (0,0),(1,0),(0,1) 
  using `n` nodes on the x- and y-axis.
"""
UnitTriangle(n::Integer) = UnitSimplex(2, n)

"""

    UnitTetrahedron(n::Integer)

  Create a mesh discretizing the unit tetrahedron 
  defined by the nodes (0,0,0),(1,0,0),(0,1,0),(0,0,1) 
  using `n` nodes on the x-, y- and z-axis.
"""
UnitTetrahedron(n::Integer) = UnitSimplex(3, n)

