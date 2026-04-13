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
            post_select_measure_0!(ψ,i,L)
        end
    end
end

function post_select_measure_all_1!(ψ::Vector{ComplexF64},L::Integer,p::Float64)
    for i=1:L
        if rand()<p
            post_select_measure_1!(ψ,i,L)
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

function average_entropies(ψ::Vector{ComplexF64},L::Int64,ns::Vector{Float64})
    MI = zeros(Float64,5)
    SA = zeros(Float64,5)
    SB = zeros(Float64,5)
    SC = zeros(Float64,5)
    SAB = zeros(Float64,5)
    SAC = zeros(Float64,5)
    SBC = zeros(Float64,5)
    SABC = zeros(Float64,5)
    for jj=1:(div(L,2))
        A = collect(jj:(jj+div(L,4)-1))
        B = mod1.(collect((jj+div(L,4)):(jj+2*div(L,4)-1)),L)
        C = mod1.(collect((jj+2*div(L,4)):(jj+3*div(L,4)-1)),L)
        r = HaarRandomCircuit.top_mutual_information(ψ,A,B,C,L,ns)
        MI .+= r[1]
        SA .+= r[2]
        SB .+= r[3]
        SC .+= r[4]
        SAB .+= r[5]
        SAC .+= r[6]
        SBC .+= r[7]
        SABC .+= r[8]
    end
    MI ./= div(L,2)
    SA ./= div(L,2)
    SB ./= div(L,2)
    SC ./= div(L,2)
    SAB ./= div(L,2)
    SAC ./= div(L,2)
    SBC ./= div(L,2)
    SABC ./= div(L,2)
    return MI, SA, SB, SC, SAB, SAC, SBC, SABC
end

function fullcircuit!(ψ::Vector{ComplexF64},
    L::Integer,p::Float64,T::Integer;BCs="open",MI,SA,SB,SC,SAB,SAC,SBC,SABC)

    U = Array{ComplexF64}(undef,4,4)
    temp_vec1 = Array{ComplexF64}(undef,4)
    temp_vec2 = Array{ComplexF64}(undef,4)

    r = average_entropies(ψ,L,[0,1,2,3,Inf])
    MI[1,:] = r[1]
    SA[1,:] = r[2]
    SB[1,:] = r[3]
    SC[1,:] = r[4]
    SAB[1,:] = r[5]
    SAC[1,:] = r[6]
    SBC[1,:] = r[7]
    SABC[1,:] = r[8]

    for t=1:T
        even_gates!(ψ,L,U,temp_vec1,temp_vec2)
        r1 = average_entropies(ψ,L,[0,1,2,3,Inf])
        post_select_measure_all_0!(ψ,L,p)
        odd_gates!(ψ,L,U,temp_vec1,temp_vec2;BCs=BCs)
        r2 = average_entropies(ψ,L,[0,1,2,3,Inf])
        if BCs=="open"
            post_select_odd_measure_all_0!(ψ,L,p)
        else
            post_select_measure_all_0!(ψ,L,p)
        end

        MI[t+1,:] = (r1[1] + r2[1])/2
        SA[t+1,:] = (r1[2] + r2[2])/2
        SB[t+1,:] = (r1[3] + r2[3])/2
        SC[t+1,:] = (r1[4] + r2[4])/2
        SAB[t+1,:] = (r1[5] + r2[5])/2
        SAC[t+1,:] = (r1[6] + r2[6])/2
        SBC[t+1,:] = (r1[7] + r2[7])/2
        SABC[t+1,:] = (r1[8] + r2[8])/2
    end
end




function fullcircuit_OE!(ψ::Vector{ComplexF64},
    L::Integer,p::Float64,T::Integer;BCs="open",MI,SA,SB,SC,SAB,SAC,SBC,SABC)

    U = Array{ComplexF64}(undef,4,4)
    temp_vec1 = Array{ComplexF64}(undef,4)
    temp_vec2 = Array{ComplexF64}(undef,4)

    r = average_entropies(ψ,L,[0,1,2,3,Inf])
    MI[1,:] = r[1]
    SA[1,:] = r[2]
    SB[1,:] = r[3]
    SC[1,:] = r[4]
    SAB[1,:] = r[5]
    SAC[1,:] = r[6]
    SBC[1,:] = r[7]
    SABC[1,:] = r[8]

    

    for t=1:T
        odd_gates!(ψ,L,U,temp_vec1,temp_vec2;BCs=BCs)
        r1 = average_entropies(ψ,L,[0,1,2,3,Inf])
        if BCs=="open"
            post_select_odd_measure_all_0!(ψ,L,p)
        else
            post_select_measure_all_0!(ψ,L,p)
        end
        even_gates!(ψ,L,U,temp_vec1,temp_vec2)
        r2 = average_entropies(ψ,L,[0,1,2,3,Inf])
        post_select_measure_all_0!(ψ,L,p)
        MI[t+1,:] = (r1[1] + r2[1])/2
        SA[t+1,:] = (r1[2] + r2[2])/2
        SB[t+1,:] = (r1[3] + r2[3])/2
        SC[t+1,:] = (r1[4] + r2[4])/2
        SAB[t+1,:] = (r1[5] + r2[5])/2
        SAC[t+1,:] = (r1[6] + r2[6])/2
        SBC[t+1,:] = (r1[7] + r2[7])/2
        SABC[t+1,:] = (r1[8] + r2[8])/2
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
        S = 1/(1-n)*log2(sum([p^n for p=ps]))
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
