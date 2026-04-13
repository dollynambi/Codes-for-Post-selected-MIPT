module HaarRandomCircuit

using Combinatorics, LinearAlgebra, Random

function initialize_state(dim::Integer)
    ψ = randn(ComplexF64,dim)
    normalize!(ψ)
end

function randomHaar!(U::Array{ComplexF64,2})
    # Look into how to make Haar from QR decomposition
    dim, dim = size(U)
    randn!(U)
    Q, R = qr(U)
    L = broadcast( x -> abs(x)>0 ? sign(x) : 0, Diagonal(R))
    U .= Matrix(Q)*L
end

# function gate_and_measure!(occupied_Fock_space::Array{Int64},ψ::Array{ComplexF64},U::Array{ComplexF64},i::Integer,
#     temp_vec1::Vector{ComplexF64},temp_vec2::Vector{ComplexF64},L::Integer)
#
#
#     for i=eachindex(occupied_Fock_space)
#         nR = mod(n,2^i)
#         nL = div(n,2^(i+1))
#
#         n1 = nR + nL*2^(i+1)
#         n2 = nR + nL*2^(i+1) + 2^(i-1)
#         n3 = nR + nL*2^(i+1) + 2^(i)
#         n4 = nR + nL*2^(i+1) + 2^(i) + 2^(i-1)
#
#         ind1 = findfirst(x -> x==n1,occupied_Fock_space)
#         ind2 = findfirst(x -> x==n2,occupied_Fock_space)
#         ind3 = findfirst(x -> x==n3,occupied_Fock_space)
#         ind4 = findfirst(x -> x==n4,occupied_Fock_space)
#
#         if ind1>=i && ind2>=i && ind3>=i && ind4 >=i
#
#
#
#
#     for nL=0:2^(L-1-i)-1, nR=0:2^(i-1)-1
#         ind1 = nR + nL*2^(i+1) + 1
#         ind2 = nR + nL*2^(i+1) + 2^(i-1) + 1
#         ind3 = nR + nL*2^(i+1) + 2^(i) + 1
#         ind4 = nR + nL*2^(i+1) + 2^(i) + 2^(i-1) + 1
#
#         temp_vec1[1] = ψ[ind1]
#         temp_vec1[2] = ψ[ind2]
#         temp_vec1[3] = ψ[ind3]
#         temp_vec1[4] = ψ[ind4]
#
#
#
# end


function gate!(ψ::Vector{ComplexF64},U::Array{ComplexF64,2},i::Integer,
    temp_vec1::Vector{ComplexF64},temp_vec2::Vector{ComplexF64},L::Integer)

    randomHaar!(U)
    for nL=0:2^(L-1-i)-1, nR=0:2^(i-1)-1
        ind1 = nR + nL*2^(i+1) + 1
        ind2 = nR + nL*2^(i+1) + 2^(i-1) + 1
        ind3 = nR + nL*2^(i+1) + 2^(i) + 1
        ind4 = nR + nL*2^(i+1) + 2^(i) + 2^(i-1) + 1

        temp_vec1[1] = ψ[ind1]
        temp_vec1[2] = ψ[ind2]
        temp_vec1[3] = ψ[ind3]
        temp_vec1[4] = ψ[ind4]

        mul!(temp_vec2,U,temp_vec1)

        ψ[ind1] = temp_vec2[1]
        ψ[ind2] = temp_vec2[2]
        ψ[ind3] = temp_vec2[3]
        ψ[ind4] = temp_vec2[4]

    end
end

function gate_periodic!(ψ::Vector{ComplexF64},U::Array{ComplexF64,2},
    temp_vec1::Vector{ComplexF64},temp_vec2::Vector{ComplexF64},L::Integer)

    randomHaar!(U)
    for n = 0:(2^(L-2)-1)
        ind1 = 2n + 1
        ind2 = 2n + 1 + 1
        ind3 = 2n + 2^(L-1) + 1
        ind4 = 2n + 2^(L-1) + 1 + 1

        temp_vec1[1] = ψ[ind1]
        temp_vec1[2] = ψ[ind2]
        temp_vec1[3] = ψ[ind3]
        temp_vec1[4] = ψ[ind4]

        mul!(temp_vec2,U,temp_vec1)

        ψ[ind1] = temp_vec2[1]
        ψ[ind2] = temp_vec2[2]
        ψ[ind3] = temp_vec2[3]
        ψ[ind4] = temp_vec2[4]

    end
end

# Always Projecting the qubit to the |0> state
function post_select_measure_0!(ψ::Vector{ComplexF64},i::Integer,L::Integer)

    for nL=0:(2^(L-i)-1), nR=(0:2^(i-1)-1)
        n = nR+2^(i)*nL+2^(i-1)+1
        ψ[n] = 0.0
    end
    normalize!(ψ)
end


# Always Projecting the qubit to the |0> state
function post_select_measure_1!(ψ::Vector{ComplexF64},i::Integer,L::Integer)

    for nL=0:(2^(L-i)-1), nR=(0:2^(i-1)-1)
        n = nR+2^(i)*nL +1
        ψ[n] = 0.0
    end 
    normalize!(ψ)
end

function post_select_measure_all_0!(ψ::Vector{ComplexF64},L::Integer,p::Float64)
    for i=1:L
        if rand()<p
            post_select_measure!(ψ,i,L)
        end
    end
end

function post_select_measure_all_1!(ψ::Vector{ComplexF64},L::Integer,p::Float64)
    for i=1:L
        if rand()<p
            post_select_measure!(ψ,i,L)
        end
    end
end

function post_select_odd_measure_all_0!(ψ::Vector{ComplexF64},L::Integer,p::Float64)
    for i=2:(L-1)
        if rand()<p
           post_select_measure_0!(ψ,i,L)
        end
    end
end

function post_select_odd_measure_all_1!(ψ::Vector{ComplexF64},L::Integer,p::Float64)
    for i=2:(L-1)
        if rand()<p
           post_select_measure_0!(ψ,i,L)
        end
    end
end


function fullcircuit!(ψ::Vector{ComplexF64},
    L::Integer,p::Float64,T::Integer;BCs="open")

    U = Array{ComplexF64}(undef,4,4)
    temp_vec1 = Array{ComplexF64}(undef,4)
    temp_vec2 = Array{ComplexF64}(undef,4)

    for t=1:T
        one_time_step!(ψ,L,p,U,temp_vec1,
            temp_vec2; BCs=BCs)
    end
end

function fullcircuit_OE!(ψ::Vector{ComplexF64},
    L::Integer,p::Float64,T::Integer;BCs="open")

    U = Array{ComplexF64}(undef,4,4)
    temp_vec1 = Array{ComplexF64}(undef,4)
    temp_vec2 = Array{ComplexF64}(undef,4)

    for t=1:T
        one_time_step_OE!(ψ,L,p,U,temp_vec1,
            temp_vec2; BCs=BCs)
    end
end


function one_time_step!(ψ::Vector{ComplexF64},L::Integer,p::Float64,
    U::Array{ComplexF64,2},temp_vec1::Vector{ComplexF64},temp_vec2::Vector{ComplexF64};
    BCs="open")

    even_gates!(ψ,L,U,temp_vec1,temp_vec2)
    post_select_measure_all_0!(ψ,L,p)
    odd_gates!(ψ,L,U,temp_vec1,temp_vec2;BCs=BCs)
    if BCs=="open"
        post_select_odd_measure_all_0!(ψ,L,p)
    else
        post_select_measure_all_0!(ψ,L,p)
    end

end

function one_time_step_OE!(ψ::Vector{ComplexF64},L::Integer,p::Float64,
    U::Array{ComplexF64,2},temp_vec1::Vector{ComplexF64},temp_vec2::Vector{ComplexF64};
    BCs="open")

    odd_gates!(ψ,L,U,temp_vec1,temp_vec2;BCs=BCs)
    if BCs=="open"
        post_select_odd_measure_all_0!(ψ,L,p)
    else
        post_select_measure_all_0!(ψ,L,p)
    end
    even_gates!(ψ,L,U,temp_vec1,temp_vec2)
    post_select_measure_all_0!(ψ,L,p)

end

function even_gates!(ψ::Vector{ComplexF64},L::Integer,
    U::Array{ComplexF64,2},temp_vec1::Vector{ComplexF64},
    temp_vec2::Vector{ComplexF64})

    for i=1:2:L
        gate!(ψ,U,i,temp_vec1,temp_vec2,L)
    end
end

function odd_gates!(ψ::Vector{ComplexF64},L::Integer,
    U::Array{ComplexF64,2},temp_vec1::Vector{ComplexF64},
    temp_vec2::Vector{ComplexF64}; BCs="open")

    for i=2:2:(L-1)
        gate!(ψ,U,i,temp_vec1,temp_vec2,L)
    end

    if BCs=="periodic"
        gate_periodic!(ψ,U,temp_vec1,temp_vec2,L)
    end
end


function svd_full(ψ::Vector{ComplexF64},A::Vector{Int64},L::Int64)

    Schmidt_values = Array{Float64}(undef,2^length(A))
    B = [j for j=1:L if !(j in A)]

    ΨSVD = Array{ComplexF64}(undef,2^length(A),2^length(B))

    for RA=CartesianIndices(Tuple(0:1 for j=1:length(A)))
        for RB=CartesianIndices(Tuple(0:1 for j=1:length(B)))
            n = dot(Tuple(RA),2 .^(A .- 1)) + dot(Tuple(RB),2 .^(B .- 1) ) + 1
            ii = dot(Tuple(RA),2 .^(0:(length(A)-1))) + 1
            jj = dot(Tuple(RB),2 .^(0:(length(B)-1))) + 1
            ΨSVD[ii,jj] = ψ[n]
        end
    end

    return svdvals(ΨSVD)
end

function Renyi_entropy(n::Number,ps)
    S = 0
    if n==0
        S = log2(sum(ps .> eps(eltype(ps))))
    elseif n==1
        S = -sum([p*log2(p) for p=ps if 0<p<=1])
    elseif n==Inf
        S = -log2(maximum(ps))
    else
        for p in ps
            S += p^n
        end
        S = log2(S)/(1-n)
    end
    return S
end

function mutual_information(ψ::Vector{ComplexF64},
    A::Vector{Int64},B::Vector{Int64},L::Integer,ns::Vector{Float64})

    pA = abs2.(svd_full(ψ,A,L))
    pB = abs2.(svd_full(ψ,B,L))
    pAB = abs2.(svd_full(ψ,unique([A;B]),L))

    SA = Array{Float64}(undef,length(ns))
    SB = Array{Float64}(undef,length(ns))
    SAB = Array{Float64}(undef,length(ns))
    for i = eachindex(ns)
        n = ns[i]
        SA[i] = Renyi_entropy(n,pA)
        SB[i] = Renyi_entropy(n,pB)
        SAB[i] = Renyi_entropy(n,pAB)
    end

    MI = SA .+ SB .- SAB
    return MI
end

function top_mutual_information(ψ::Vector{ComplexF64},
    A::Vector{Int64},B::Vector{Int64},C::Vector{Int64},L::Integer,ns::Vector{Float64})

    pA = abs2.(svd_full(ψ,A,L))
    pB = abs2.(svd_full(ψ,B,L))
    pC = abs2.(svd_full(ψ,C,L))
    pAB = abs2.(svd_full(ψ,unique([A;B]),L))
    pAC = abs2.(svd_full(ψ,unique([A;C]),L))
    pBC = abs2.(svd_full(ψ,unique([B;C]),L))
    pABC = abs2.(svd_full(ψ,unique([A;B;C]),L))

    SA = Array{Float64}(undef,length(ns))
    SB = Array{Float64}(undef,length(ns))
    SC = Array{Float64}(undef,length(ns))
    SAB = Array{Float64}(undef,length(ns))
    SAC = Array{Float64}(undef,length(ns))
    SBC = Array{Float64}(undef,length(ns))
    SABC = Array{Float64}(undef,length(ns))

    for i = eachindex(ns)
        n = ns[i]
        SA[i] = Renyi_entropy(n,pA)
        SB[i] = Renyi_entropy(n,pB)
        SC[i] = Renyi_entropy(n,pC)
        SAB[i] = Renyi_entropy(n,pAB)
        SAC[i] = Renyi_entropy(n,pAC)
        SBC[i] = Renyi_entropy(n,pBC)
        SABC[i] = Renyi_entropy(n,pABC)
    end 

    MI = SA .+ SB .+ SC .- SAB .- SAC .- SBC .+ SABC
    return MI, SA, SB, SC, SAB, SAC, SBC, SABC

end

mutual_information(ψ::Vector{ComplexF64}, A::Vector{Int64}, B::Vector{Int64}, L::Integer) = mutual_information(ψ, A, B, L, [1])[1]


end
