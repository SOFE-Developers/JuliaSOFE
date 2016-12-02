#-------------------#
# Computing Entries #
#-------------------#

typealias op ScalarOperator
typealias Op VectorOperator
typealias OP MatrixOperator

typealias data Union{ScalarCoefficient, FunctionCoefficient}
typealias Data Union{VectorCoefficient, FunctionCoefficient}
typealias DATA Union{MatrixCoefficient, FunctionCoefficient}

# Bilinear Operators
# -------------------

# scalar coeff * scalar oper (local) * scalar oper (local)
function fill_entries!{T<:Real,Tc<:data,Tu<:op,Tv<:op}(::BilinearForm{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3}, # nE x nBi x nBj
                                                       C::AbstractArray{T,2}, # nE x nP
                                                       U::AbstractArray{T,3}, # nBj x nP x nC
                                                       V::AbstractArray{T,3}, # nBi x nP x nC
                                                       w::AbstractArray{T,1}, # nP
                                                       D::AbstractArray{T,2}) # nE x nP
    nE, nBi, nBj = size(E)
    nBj, nP, nC = size(U); @assert nC == 1
    nBi, nP, nC = size(V); @assert nC == 1

    for ip = 1:nP
        for jb = 1:nBj
            for ib = 1:nBi
                for ie = 1:nE
                    E[ie,ib,jb] += C[ie,ip] * U[jb,ip,1] * V[ib,ip,1] * w[ip] * D[ie,ip]
                end
            end
        end
    end

    return nothing
end

# scalar coeff * vector oper (global) * vector oper (global)
function fill_entries!{T<:Real,Tc<:data,Tu<:Op,Tv<:Op}(::BilinearForm{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3}, # nE x nBi x nBj
                                                       C::AbstractArray{T,2}, # nE x nP
                                                       U::AbstractArray{T,4}, # nE x nBj x nP x nC
                                                       V::AbstractArray{T,4}, # nE x nBi x nP x nC
                                                       w::AbstractArray{T,1}, # nP
                                                       D::AbstractArray{T,2}) # nE x nP
    nE, nBi, nBj = size(E)
    nE, nB, nP, nC = size(U)

    for ic = 1:nC
        for ip = 1:nP
            for jb = 1:nBj
                for ib = 1:nBi
                    for ie = 1:nE
                        E[ie,ib,jb] += C[ie,ip] * U[ie,jb,ip,ic] * V[ie,ib,ip,ic] * w[ip] * D[ie,ip]
                    end
                end
            end
        end
    end
    return nothing
end

# vector coeff * vector oper (global) * scalar oper (local)
function fill_entries!{T<:Real,Tc<:Data,Tu<:Op,Tv<:op}(::BilinearForm{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3}, # nE x nBi x nBj
                                                       C::AbstractArray{T,3}, # nE x nP x nD
                                                       U::AbstractArray{T,4}, # nE x nBj x nP x nD
                                                       V::AbstractArray{T,2}, # nBi x nP
                                                       w::AbstractArray{T,1}, # nP
                                                       D::AbstractArray{T,2}) # nE x nP
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for id = 1:nW
                        E[ie,ib,jb] += C[ie,ip,id] * U[ie,jb,ip,id] * V[ib,ip] * w[ip] * D[ie,ip]
                    end
                end
            end
        end
    end
    return nothing
end

# matrix coeff * vector oper (global) * vector oper (global)
function fill_entries!{T<:Real,Tc<:DATA,Tu<:Op,Tv<:Op}(::BilinearForm{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3}, # nE x nBi x nBj
                                                       C::AbstractArray{T,4}, # nE x nP x nCi x nCj
                                                       U::AbstractArray{T,4}, # nE x nBj x nP x nCj
                                                       V::AbstractArray{T,4}, # nE x nBi x nP x nCi
                                                       w::AbstractArray{T,1}, # nP
                                                       D::AbstractArray{T,2}) # nE x nP
    nE, nBi, nBj = size(E)
    nE, nP, nCi, nCj = size(C)
    nE, nB, nP, nCj = size(U)
    nE, nB, nP, nCi = size(V)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for jc = 1:nC
                        for ic = 1:nC
                            E[ie,ib,jb] += C[ie,ip,ic,jc] * U[ie,jb,ip,jc] * V[ie,ib,ip,ic] * w[ip] * D[ie,ip]
                        end
                    end
                end
            end
        end
    end
    return nothing
end

# Linear Operators
# -----------------

# scalar coeff * scalar operator (local)
function fill_entries!{T<:Real,Tc<:data,Tv<:op}(::LinearForm{Tc,Tv},
                                                E::AbstractArray{T,2}, # nE x nB
                                                C::AbstractArray{T,2}, # nE x nP
                                                V::AbstractArray{T,3}, # nB x nP x nC
                                                w::AbstractArray{T,1}, # nP
                                                D::AbstractArray{T,2}) # nE x nP
    nE, nB = size(E)
    nB, nP, nC = size(V); @assert nC == 1
    for ib = 1:nB
        for ie = 1:nE
            for ip = 1:nP
                E[ie,ib] += C[ie,ip] * V[ib,ip] * w[ip] * D[ie,ip]
            end
        end
    end

    return nothing
end

# vector coeff * vector operator (local)
function fill_entries!{T<:Real,Tc<:Data,Tv<:Op}(::LinearForm{Tc,Tv},
                                                E::AbstractArray{T,2}, # nE x nB
                                                C::AbstractArray{T,3}, # nE x nP x nC
                                                V::AbstractArray{T,3}, # nB x nP x nC
                                                w::AbstractArray{T,1}, # nP
                                                D::AbstractArray{T,2}) # nE x nP
    nE, nB = size(E)
    nB, nP, nC = size(V)
    # nW, = size(C)
    # nP = size(w, 1)
    for ib = 1:nB
        for ie = 1:nE
            for ip = 1:nP
                for ic = 1:nC
                    E[ie,ib] += C[ie,ip,ic] * V[ib,ip,ic] * w[ip] * D[ie,ip]
                end
            end
        end
    end

    return nothing
end

# vector coeff * vector operator (global)
function fill_entries!{T<:Real,Tc<:Data,Tv<:Op}(::LinearForm{Tc,Tv},
                                                E::AbstractArray{T,2}, # nE x nB
                                                C::AbstractArray{T,3}, # nE x nP x nC
                                                V::AbstractArray{T,4}, # nE x nB x nP x nC
                                                w::AbstractArray{T,1}, # nP
                                                D::AbstractArray{T,2}) # nE x nP
    nE, nB = size(E)
    nE, nB, nP, nC = size(V)
    # nW, = size(C)
    # nP = size(w, 1)
    for ib = 1:nB
        for ie = 1:nE
            for ip = 1:nP
                for ic = 1:nC
                    E[ie,ib] += C[ie,ip,ic] * V[ie,ib,ip,ic] * w[ip] * D[ie,ip]
                end
            end
        end
    end

    return nothing
end

