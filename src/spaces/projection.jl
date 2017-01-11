using ..Quadrature

export project

# need this
include(joinpath("..", "operators", "entries.jl"))

function project(fes::FESpace, f::Function)
    if issimp(element(fes))
        qrule = GaussQuadSimp(2*order(fes.element), dimension(fes.mesh))
    else
        qrule = GaussQuadOrth(2*order(fes.element), dimension(fes.mesh))
    end

    qp = qnodes(qrule, dimension(element(fes)))
    qw = qweights(qrule, dimension(element(fes)))

    dMap = dofMap(fes, dimension(element(fes)))
    ndof = nDoF(fes)

    nB, nE = size(dMap)
    nP, nD = size(qp)

    dofI = zeros(Int, nE, nB, nB)
    dofJ = zeros(Int, nE, nB, nB)
    dofK = zeros(Int, nE, nB)
    for jb = 1:nB
        for ib = 1:nB
            for ie = 1:nE
                dofI[ie,ib,jb] = dMap[ib,ie]
                dofJ[ie,ib,jb] = dMap[jb,ie]
                jb == 1 && (dofK[ie,ib] = dMap[ib,ie])
            end
        end
    end

    C = ones(eltype(qp), nE, nP)
    F = evaluate(f, qp, mesh(fes))
    U = V = evalBasis(element(fes), qp, 0)
    D = abs(evalJacobianDeterminat(mesh(fes), qp))

    Ea = zeros(eltype(qp), nE, nB, nB); fill_entries!(Ea, C, U, V, qw, D)
    El = zeros(eltype(qp), nE, nB);     fill_entries!(El, F, V, qw, D)

    A = sparse(dofI[:], dofJ[:], Ea[:], ndof, ndof)
    L = full(sparsevec(dofK[:], El[:], ndof))
    P = A\L
    return P
end
