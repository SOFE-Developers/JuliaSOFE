using Jacobi
using Polynomials

fac = factorial

"""

    poly_lobatto(k::Integer)

Return the `k`-th Lobatto polynomial.

"""
# """
# If ``k \geq 2`` the coefficients are computed using the relation

# ```math
# L_{k-1}(x) = \frac{1}{2k-1} \frac{d}{dx} (L_{k}(x) - L_{k-2})
# ```

# which after integration gives the following computable recurrence relation

# ```math
# l_k(x) = \frac{1}{\sqrt{4k-2}} (L_{k}(x) - L_{k-2}(x))
# ```

# The Lobatto polynomials are defined by

# ```math
# \begin{eqnarray}
#     l_0(x) = -\frac{1}{2} x + \frac{1}{2} \\
#     l_1(x) = \frac{1}{2} x + \frac{1}{2}  \\
#     l_k(x) = \frac{1}{\| L_{k-1} \|_2} \int_{-1}^{x} L_{k-1}(s) ds,\ k \geq 2
# \end{eqnarray}
# ```

# with ``L_k`` being `k`-th order Legendre polynomial and
# ``\| L_{k-1} \|_2 = \sqrt{\frac{2}{2 k - 1}}``.
# """
function poly_lobatto(k::Integer)
    k < 0 && error("Polynomial order has to be non-negative (", k, ")")

    if k == 0
        return Poly([0.5, -0.5])
    elseif k == 1
        return Poly([0.5, 0.5])
    else
        l_k = poly_legendre(k)
        l_km2 = poly_legendre(k-2)
        return (l_k - l_km2) / sqrt(4k-2)
    end
end

"""

    poly_dlobatto(k::Integer, n::Integer=1)

Return the `n`-th derivative of the `k`-th Lobatto polynomial
"""
# using the relation

# """```math
# \frac{d^n}{dx^n} l_k(x) = \frac{1}{2^n \sqrt{4k-2}} \left(\frac{(k+n)!}{k!} J_{k-n}^{(n,n)}(x) - \frac{(k+n-2)!}{(k-2)!} J_{k-n-2}^{(n,n)}(x) \right)
# ```

# with ``J_k^{(n,n)}`` being the `k`-th order Jacobi polynomial.
# """
function poly_dlobatto(k::Integer, n::Integer=1)
    if k == 0 && n == 0
        return Poly([0.5, -0.5])
    elseif k == 0 && n == 1
        return Poly([-0.5])
    elseif k == 1 && n == 0
        return Poly([0.5, 0.5])
    elseif k == 1 && n == 1
        return Poly([0.5])
    elseif k in (0, 1) && n > 1
        return Poly([0.0])
    elseif k >= 1 && k >= n
        if k-n-2 > 0
            J_kn = poly_jacobi(k-n, n, n)
            J_kn2 = poly_jacobi(k-n-2, n, n)

            return (fac(k+n)/fac(k) * J_kn - fac(k+n-2)/fac(k-2) * J_kn2) / (2^n * sqrt(4k-2))
        else
            dn_l_k = poly_lobatto(k)
            for i = 1:n
                dn_l_k = polyder(dn_l_k)
            end
            return dn_l_k
        end
    else
        error("Invalid pair of polynomial and derivation order (", k, n, ")")
    end
end

"""

    poly_lobatto_kernel(k::Integer)

Return the `k`-th Lobatto kernel polynomial.

"""
# """
# The kernel polynomials are obtained by decomposing 
# the Lobatto polynomials for `k \geq 2` into products
# of the form

# ```math
# l_{k}(x) = l_0(x) l_1(x) \varphi_{k-2}(x)
# ```
# """
function poly_lobatto_kernel(k::Integer)
    p = poly_lobatto(k+2)
    q = Poly([0.25, 0.0, -0.25])
    return div(p, q)
end
