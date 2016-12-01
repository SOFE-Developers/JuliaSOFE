#-------------------#
# Computing Entries #
#-------------------#

typealias op ScalarOperator
typealias Op VectorOperator
typealias OP MatrixOperator

typealias data Union{ScalarCoefficient, FunctionCoefficient}
typealias Data Union{VectorCoefficient, FunctionCoefficient}
typealias DATA Union{MatrixCoefficient, FunctionCoefficient}
#typealias fdat FunctionCoefficient

# Bilinear Operators
# -------------------

# scalar coeff * scalar oper (local) * scalar oper (local)
function fill_entries!{T<:Real,Tc<:data,Tu<:op,Tv<:op}(::BilinearForm{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3}, # nE x nBi x nBj
                                                       C::AbstractArray{T,2}, # nE x nP
                                                       U::AbstractArray{T,2}, # nBj x nP
                                                       V::AbstractArray{T,2}, # nBi x nP
                                                       w::AbstractArray{T,1}, # nP
                                                       D::AbstractArray{T,2}) # nE x nP
    nE, nBi, nBj = size(E)
    nBj, nP = size(U)
    nBi, nP = size(V)
    nP = size(w, 1)

    for ip = 1:nP
        for jb = 1:nBj
            for ib = 1:nBi
                for ie = 1:nE
                    E[ie,ib,jb] += C[ie,ip] * U[jb,ip] * V[ib,ip] * w[ip] * D[ie,ip]
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
                                                       U::AbstractArray{T,4}, # nE x nBj x nP x nD
                                                       V::AbstractArray{T,4}, # nE x nBi x nP x nD
                                                       w::AbstractArray{T,1}, # nP
                                                       D::AbstractArray{T,2}) # nE x nP
    nE, nBi, nBj = size(E)
    nE, nB, nP, nD = size(U)

    for id = 1:nD
        for ip = 1:nP
            for jb = 1:nBj
                for ib = 1:nBi
                    for ie = 1:nE
                        E[ie,ib,jb] += C[ie,ip] * U[ie,jb,ip,id] * V[ie,ib,ip,id] * w[ip] * D[ie,ip]
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
                                                       C::AbstractArray{T,4}, # nE x nP x nD x nD
                                                       U::AbstractArray{T,4}, # nE x nBj x nP x nD
                                                       V::AbstractArray{T,4}, # nE x nBi x nP x nD
                                                       w::AbstractArray{T,1}, # nP
                                                       D::AbstractArray{T,2}) # nE x nP
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for jd = 1:nW
                        for id = 1:nW
                            E[ie,ib,jb] += C[ie,ip,id,jd] * U[ie,jb,ip,jd] * V[ie,ib,ip,id] * w[ip] * D[ie,ip]
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
                                                E::AbstractArray{T,2},
                                                C::AbstractArray{T,2},
                                                V::AbstractArray{T,2},
                                                w::AbstractArray{T,1},
                                                D::AbstractArray{T,2})
    nE, nB = size(E)
    nP = size(w, 1)
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
                                                E::AbstractArray{T,2},
                                                C::AbstractArray{T,3},
                                                V::AbstractArray{T,3},
                                                w::AbstractArray{T,1},
                                                D::AbstractArray{T,2})
    nE, nB = size(E)
    nB, nP, nW = size(V)
    # nW, = size(C)
    # nP = size(w, 1)
    for ib = 1:nB
        for ie = 1:nE
            for ip = 1:nP
                for id = 1:nW
                    E[ie,ib] += C[ie,ip,id] * V[ib,ip,id] * w[ip] * D[ie,ip]
                end
            end
        end
    end

    return nothing
end

# vector coeff * vector operator (global)
function fill_entries!{T<:Real,Tc<:Data,Tv<:Op}(::LinearForm{Tc,Tv},
                                                E::AbstractArray{T,2},
                                                C::AbstractArray{T,3},
                                                V::AbstractArray{T,4},
                                                w::AbstractArray{T,1},
                                                D::AbstractArray{T,2})
    nE, nB = size(E)
    nE, nB, nP, nW = size(V)
    # nW, = size(C)
    # nP = size(w, 1)
    for ib = 1:nB
        for ie = 1:nE
            for ip = 1:nP
                for id = 1:nW
                    E[ie,ib] += C[ie,ip,id] * V[ie,ib,ip,id] * w[ip] * D[ie,ip]
                end
            end
        end
    end

    return nothing
end

