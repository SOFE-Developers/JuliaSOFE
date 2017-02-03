using CMesh
using JuliaSOFE

include("randcoeff.jl")

#--------------------------#
# Random Coefficient Field #
#--------------------------#
epsilon = 1/2^6
beta = 0.2

# X = checkerboard(Int(1/epsilon))
X = random_checkerboard(Int(1/epsilon), false)
# X = random_checkerboard(Int(1/epsilon), true, beta=beta)

function inside{T<:Real}(xs::AbstractVector{T})
    # c = Vector{Int}(ceil(mod(xs, 1)/epsilon))
    c = Vector{Int}(ceil(xs))
    #return all(iseven.(c)) | all(isodd.(c))
    return X[c[1], c[2]]
end

# function a{T<:Real}(xs::AbstractArray{T,2})
#     A = zeros(T, size(xs, 1))
#     x = zeros(T, size(xs, 2))
#     for ip = 1:size(xs,1)
#         for id = 1:size(xs, 2)
#             x[id] = xs[ip,id]
#         end
#         A[ip] = inside(x) ? convert(T, 100) : convert(T, 1)
#     end
#     return A
# end
function a{T<:Real}(xs::AbstractArray{T,2})
    A = zeros(T, (size(xs, 1), 2, 2))
    x = zeros(T, size(xs, 2))
    for ip = 1:size(xs,1)
        for id = 1:size(xs, 2)
            x[id] = xs[ip,id]
        end
        A[ip,1,1] = A[ip,2,2] = inside(x) ? convert(T, 100) : convert(T, 1)
    end
    return A
end
function ae{T<:Real}(xs::AbstractArray{T,2}; ɛ::Real=epsilon)
    return a(xs/ɛ)
end

if true
    #---------#
    # Problem #
    #---------#

    # Mesh #
    #------#
    N = max(Int(1/epsilon)+1, 256)
    m = UnitSquare(N)

    # Function Space #
    #----------------#
    e = P1(2)
    fes = FESpace(m, e, x->true, x->zero(x))

    lhs = BilinearForm(Grad, Grad, fes, ae)
    rhs = LinearForm(id, fes, 1)
    pde = PDE(lhs, rhs)

    Ue = solve(pde)
end

if true
    #---------------#
    # Cell Problems #
    #---------------#
    
    # Mesh #
    #------#
    N = max(Int(1/epsilon)+1, 256)
    m = UnitSquare(N)
    
    # Function Space #
    #----------------#
    e = P1(2)
    fes = FESpace(m, e, x->false, x->zero(x))
    
    lhs = BilinearForm(Grad, Grad, fes, ae)
    rhs1 = LinearForm(Grad, fes, x -> -ae(x)[:,:,1])
    rhs2 = LinearForm(Grad, fes, x -> -ae(x)[:,:,2])
    pde1 = PDE(lhs, rhs1)
    pde2 = PDE(lhs, rhs2)
    
    #-------------------#
    # Periodic Boundary #
    #-------------------#
    master1{T<:Real}(x::Vector{T}) = x[1] < zero(T) + sqrt(eps(T))
    slave1{T<:Real}(x::Vector{T})  = x[1] > one(T) - sqrt(eps(T))
    master2{T<:Real}(x::Vector{T}) = x[2] < zero(T) + sqrt(eps(T))
    slave2{T<:Real}(x::Vector{T})  = x[2] > one(T) - sqrt(eps(T))
    
    pbc1 = PeriodicBoundary(fes, master1, slave1)
    pbc2 = PeriodicBoundary(fes, master2, slave2)

    #--------------#
    # Zero Average #
    #--------------#
    zac = ZeroAverageConstraint(fes)

    #----------#
    # Solution #
    #----------#
    assemble!(pde1)
    A1, b1 = system(pde1)
    A1, b1 = apply(pbc1, A1, b1)
    A1, b1 = apply(pbc2, A1, b1)
    A1, b1 = apply(zac, A1, b1)

    phi1 = A1 \ b1

    assemble!(pde2)
    A2, b2 = system(pde2)
    A2, b2 = apply(pbc1, A2, b2)
    A2, b2 = apply(pbc2, A2, b2)
    A2, b2 = apply(zac, A2, b2)

    phi2 = A2 \ b2

    M1 = transform(pbc1)
    M2 = transform(pbc2)

    phi1 = M2 * (M1 * phi1)
    phi2 = M2 * (M1 * phi2)

    #------------------------#
    # Effective Coefficients #
    #------------------------#
    qp, qw = quadData(lhs.quadrule, 2)
    jac_dets = evalJacobianDeterminat(m, qp)
    
    dphi1 = evaluate(fes, phi1, qp, 1)
    dphi2 = evaluate(fes, phi2, qp, 1)
    dphi = cat(3, dphi1, dphi2)

    A = evaluate(ae, qp, m)
    Ah = zeros((2,2))
    I = reshape(eye(2), (1,1,2,2))

    for i = 1:2
        for j = 1:2
            e = I[:,:,:,i]
            v = broadcast(+, e, dphi[:,:,i,:])
            v = sum(broadcast(*, A, reshape(v, size(v)..., 1)), 4)[:,:,:,1]
            dx = broadcast(*, jac_dets, reshape(qw, 1, length(qw)))
            Ah[i,j] = sum(v[:,:,j] .* dx)
        end
    end

    #---------------------#
    # Homogenized Problem #
    #---------------------#
    eh = P1(2)
    fesh = FESpace(m, eh, x->true, x->zero(x))
    lhsh = BilinearForm(Grad, Grad, fesh, Ah)
    rhsh = LinearForm(id, fesh, 1)
    pdeh = PDE(lhsh, rhsh)

    Uh = solve(pdeh)

end
