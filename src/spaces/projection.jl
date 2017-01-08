using ..Quadrature

export project



function project(fes::FESpace, f::Function)
    if issimp(element(fes))
        qrule = GaussQuadSimp(2*order(fes.element), dimension(fes.mesh))
    else
        qrule = GaussQuadOrth(2*order(fes.element), dimension(fes.mesh))
    end

    qp = qnodes(qrule, dimension(element(fes)))
    basis = evalBasis(element(fes), qp, 0)
    
end
