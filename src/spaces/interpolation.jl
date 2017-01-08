export interpolate

function interpolate(m::Mesh, el::Element, f::Function)
    @assert isnodal(el)

    d = dimension(el)
    p = order(el)

    v = nodes(m)
    e = p > 1 ? evalReferenceMaps(m, linspace(0,1,p+1)[2:end-1]) : zeros(0,d)
    i = (p > 2 && d > 1) ? evalReferenceMaps(m, lagrangeNodesP(d,p)[p*(d+1):end,:]) : zeros(0,d)

    n = vcat(v, e, i)

    return f(n)
end
interpolate(fes::FESpace, f::Function) = interpolate(fes.mesh, fes.element, f)

