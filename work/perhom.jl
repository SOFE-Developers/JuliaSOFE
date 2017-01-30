using CMesh
using JuliaSOFE

#--------------------------#
# Random Coefficient Field #
#--------------------------#
nsamples = 20
rlo = 0.02; rhi = 0.1
cs = rhi + (1 - 2rhi) * rand(nsamples, 2)
rs = rlo + (rhi - rlo) * rand(nsamples)

sdfcs = [DCircle(cs[i,:], rs[i]) for i = 1:nsamples]
sdf = DUnion(sdfcs...)

inside{T<:Real}(xs::AbstractArray{T,2}) = sdf(xs) .< zero(T)

function a{T<:Real}(xs::AbstractArray{T,2})
    # A = zeros(T, size(xs, 1))
    A = zeros(T, size(xs, 1), 2, 2)
    x = zeros(T, size(xs, 2))
    for ip = 1:size(xs,1)
        for id = 1:size(xs, 2)
            x[id] = xs[ip,id]
        end
        A[ip,1,1] = A[ip,2,2] = sdf(x) < zero(T) ? convert(T, 100) : convert(T, 1)
    end
    return A
end

#------#
# Mesh #
#------#
m = UnitSquare(Int(2/rlo))

#----------------#
# Function Space #
#----------------#
e = P1(2)
fes = FESpace(m, e, x->false, x->zero(x))

#--------------#
# Cell Problem #
#--------------#
lhs = BilinearForm(Grad, Grad, fes, a)
rhs = LinearForm(Grad, fes, x -> a(x)[:,:,1])
pde = PDE(lhs, rhs)

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
assemble!(pde)
A, b = system(pde)
A, b = apply(pbc1, A, b)
A, b = apply(pbc2, A, b)
A, b = apply(zac, A, b)

#Uh = solve(pde)
Uh = A \ b

M1 = transform(pb1)
M2 = transform(pb2)

Uh = M2 * (M1 * Uh)

#-----------------#
# Post Processing #
#-----------------#
plot(fes, Uh)
