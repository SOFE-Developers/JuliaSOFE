#-------------------#
# Computing Entries #
#-------------------#

typealias op ScalarOperator
typealias Op VectorOperator
typealias OP MatrixOperator

typealias data ScalarCoefficient
typealias Data VectorCoefficient
typealias DATA MatrixCoefficient


# Bilinear Operators
# ===================

# Constant Coefficients
#-----------------------

function fill_entries!{T<:Real,Tc<:data,Tu<:op,Tv<:op}(::Operator{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3},
                                                       C::T,
                                                       U::AbstractArray{T,2},
                                                       V::AbstractArray{T,2},
                                                       w::AbstractArray{T,1},
                                                       D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nP = size(w, 1)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    E[ie,ib,jb] += C * U[jb,ip] * V[ib,ip] * w[ip] * D[ie,ip]
                end
            end
        end
    end
    return nothing
end

function fill_entries!{T<:Real,Tc<:data,Tu<:Op,Tv<:Op}(::Operator{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3},
                                                       C::T,
                                                       U::AbstractArray{T,4},
                                                       V::AbstractArray{T,4},
                                                       w::AbstractArray{T,1},
                                                       D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)

    for ip = 1:nP
        for id = 1:nW
            for jb = 1:nBj
                for ib = 1:nBi
                    for ie = 1:nE
                        E[ie,ib,jb] += C * U[ie,jb,ip,id] * V[ie,ib,ip,id] * w[ip] * D[ie,ip]
                    end
                end
            end
        end
    end
    return nothing
end

function fill_entries!{T<:Real,Tc<:Data,Tu<:Op,Tv<:op}(::Operator{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3},
                                                       C::AbstractArray{T,1},
                                                       U::AbstractArray{T,4},
                                                       V::AbstractArray{T,2},
                                                       w::AbstractArray{T,1},
                                                       D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for id = 1:nW
                        E[ie,ib,jb] += C[id] * U[ie,jb,ip,id] * V[ib,ip] * w[ip] * D[ie,ip]
                    end
                end
            end
        end
    end
    return nothing
end

function fill_entries!{T<:Real,Tc<:DATA,Tu<:Op,Tv<:Op}(::Operator{Tc,Tu,Tv},
                                                       E::AbstractArray{T,3},
                                                       C::AbstractArray{T,2},
                                                       U::AbstractArray{T,4},
                                                       V::AbstractArray{T,4},
                                                       w::AbstractArray{T,1},
                                                       D::AbstractArray{T,2})
    nE, nBi, nBj = size(E)
    nE, nB, nP, nW = size(U)
    for jb = 1:nBj
        for ib = 1:nBi
            for ie = 1:nE
                for ip = 1:nP
                    for jd = 1:nW
                        for id = 1:nW
                            E[ie,ib,jb] += C[id,jd] * U[ie,jb,ip,jd] * V[ie,ib,ip,id] * w[ip] * D[ie,ip]
                        end
                    end
                end
            end
        end
    end
    return nothing
end

# Non-Constant Coefficients
# --------------------------


# Linear Operators
# =================

# Constant Coefficients
# ----------------------

function fill_entries!{T<:Real,Tc<:data,Tv<:op}(::Functional{Tc,Tv},
                                                E::AbstractArray{T,2},
                                                C::T,
                                                V::AbstractArray{T,2},
                                                w::AbstractArray{T,1},
                                                D::AbstractArray{T,2})
    nE, nB = size(E)
    nP = size(w, 1)
    for ib = 1:nB
        for ie = 1:nE
            for ip = 1:nP
                E[ie,ib] += C * V[ib,ip] * w[ip] * D[ie,ip]
            end
        end
    end

    return nothing
end

