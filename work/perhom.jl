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
    A = zeros(T, size(xs, 1))
    x = zeros(T, size(xs, 2))
    for ip = 1:size(xs,1)
        for id = 1:size(xs, 2)
            x[id] = xs[ip,id]
        end
        A[ip] = sdf(x) < zero(T) ? convert(T, 100) : convert(T, 1)
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
fes = FESpace(m, e, x->true, x->zero(x))

#---------#
# Problem #
#---------#
op = BilinearForm(Grad, Grad, fes, a)
fnc = LinearForm(id, fes, 1)
pde = PDE(op, fnc)

#----------#
# Solution #
#----------#
Uh = solve(pde)

#-----------------#
# Post Processing #
#-----------------#
plot(fes, Uh)
