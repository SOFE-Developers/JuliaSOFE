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

#---------------#
# Cell Problems #
#---------------#
lhs = BilinearForm(Grad, Grad, fes, a)
rhs1 = LinearForm(Grad, fes, x -> a(x)[:,:,1])
rhs2 = LinearForm(Grad, fes, x -> a(x)[:,:,2])
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

φ1 = A1 \ b1

assemble!(pde2)
A2, b2 = system(pde2)
A2, b2 = apply(pbc1, A2, b2)
A2, b2 = apply(pbc2, A2, b2)
A2, b2 = apply(zac, A2, b2)

φ2 = A2 \ b2

M1 = transform(pbc1)
M2 = transform(pbc2)

φ1 = M2 * (M1 * φ1)
φ2 = M2 * (M1 * φ2)

#-----------------#
# Post Processing #
#-----------------#
#plot(fes, Uh)
