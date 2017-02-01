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

function a{T<:Real}(xs::AbstractArray{T,2})
    A = zeros(T, size(xs, 1))
    x = zeros(T, size(xs, 2))
    for ip = 1:size(xs,1)
        for id = 1:size(xs, 2)
            x[id] = xs[ip,id]
        end
        A[ip] = inside(x) ? convert(T, 100) : convert(T, 1)
    end
    return A
end
function ae{T<:Real}(xs::AbstractArray{T,2}; ɛ::Real=epsilon)
    return a(xs/ɛ)
end

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

Uh = solve(pde)


# #---------------#
# # Cell Problems #
# #---------------#

# # Mesh #
# #------#
# m = UnitSquare(10*Int(1/epsilon)+1)

# # Function Space #
# #----------------#
# e = P1(2)
# fes = FESpace(m, e, x->false, x->zero(x))

# lhs = BilinearForm(Grad, Grad, fes, a)
# rhs1 = LinearForm(Grad, fes, x -> a(x)[:,:,1])
# rhs2 = LinearForm(Grad, fes, x -> a(x)[:,:,2])
# pde1 = PDE(lhs, rhs1)
# pde2 = PDE(lhs, rhs2)

# #-------------------#
# # Periodic Boundary #
# #-------------------#
# master1{T<:Real}(x::Vector{T}) = x[1] < zero(T) + sqrt(eps(T))
# slave1{T<:Real}(x::Vector{T})  = x[1] > one(T) - sqrt(eps(T))
# master2{T<:Real}(x::Vector{T}) = x[2] < zero(T) + sqrt(eps(T))
# slave2{T<:Real}(x::Vector{T})  = x[2] > one(T) - sqrt(eps(T))

# pbc1 = PeriodicBoundary(fes, master1, slave1)
# pbc2 = PeriodicBoundary(fes, master2, slave2)

# #--------------#
# # Zero Average #
# #--------------#
# zac = ZeroAverageConstraint(fes)

# #----------#
# # Solution #
# #----------#
# assemble!(pde1)
# A1, b1 = system(pde1)
# A1, b1 = apply(pbc1, A1, b1)
# A1, b1 = apply(pbc2, A1, b1)
# A1, b1 = apply(zac, A1, b1)

# phi1 = A1 \ b1

# assemble!(pde2)
# A2, b2 = system(pde2)
# A2, b2 = apply(pbc1, A2, b2)
# A2, b2 = apply(pbc2, A2, b2)
# A2, b2 = apply(zac, A2, b2)

# phi2 = A2 \ b2

# M1 = transform(pbc1)
# M2 = transform(pbc2)

# phi1 = M2 * (M1 * phi1)
# phi2 = M2 * (M1 * phi2)
