include("HaarRandomCircuitAncilla.jl")
using .HaarRandomCircuitAncilla
using LinearAlgebra
using CSV, DataFrames, Random

function main(p::Float64,L::Int64,r::Int64)
    Random.seed!(r)
    s = run_periodic(p,L)
    return s
    df = DataFrame(s,:auto)
    CSV.write("spatial_correlation_$(p)_$(L)_$(r).csv",df,header=["s1", "s2", "s3","sinf"] )
end


function entangle_ancillas!(ψ::Vector{ComplexF64},L::Int64)

    ψmat = reshape(ψ,(2^(L),4))
    i = 1
    j = div(L,2)+1
    
    dim = length(ψmat[:,1])
    dimL = dim ÷ 2^j 
    dimM = 2^(j-1-i)
    dimR = 2^(i-1)

    for nL=0:(dimL-1),nM= (0:dimM-1),nR=(0:dimR-1)
        n00 = nR + 2^(i)*nM + 2^(j)*nL + 1
        n01 = n00 + 2^(i-1)
        n10 = n00 + 2^(j-1)
        n11 = n00 + 2^(i-1) + 2^(j-1)

        ψmat[n01,1],ψmat[n10,1],ψmat[n11,1] = 0, 0, 0
        ψmat[n01,2], ψmat[n10,3], ψmat[n11,4] = ψmat[n00,1], ψmat[n00,1], ψmat[n00,1]
    end
    
    ψ .= reshape(ψmat,(2^(L+2),1))
    normalize!(ψ)   
end


function entanglePBC!(ψ::Vector{ComplexF64},L::Integer)

    temp_vec1 = Array{ComplexF64,1}(undef,4)
    temp_vec2 = Array{ComplexF64,1}(undef,4)

    # What is this matrix?
    U = 1/sqrt(2)*[1 0 1 0; 0 1 0 1; 0 1 0 -1; 1 0 -1 0]
    
    # Entangling the first qubit with the (L-1)th qubit
    i=1
    for nL=0:1,nM=0:2^(L-2-i)-1
    	ind1 = nL*2^(L-1)+ nM*2^(i) + 1
        ind2 = nL*2^(L-1)+ nM*2^(i) + 2^(i-1) + 1
        ind3 = nL*2^(L-1)+ nM*2^(i) + 2^(L-2) + 1
        ind4 = nL*2^(L-1)+ nM*2^(i) + 2^(L-2) + 2^(i-1) + 1

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

    # Entangling the L/2 +1 th qubit with the Lth qubit
    i = div(L-2,2)
    for nM=0:2^(L-i-1)-1, nR=0:(2^(i-1)-1)
    	ind1 = nM*2^(i)  +nR + 1
        ind2 = nM*2^(i) +nR + 2^(i-1) + 1
        ind3 = nM*2^(i) +nR + 2^(L-1) + 1
        ind4 = nM*2^(i)+nR + 2^(L-1) + 2^(i-1) + 1

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


function run_periodic(p::Float64,L::Int64)
    
    T = 2L
    s = zeros(Float64,8T+1,4)
    num_q = L+2
    # Defining the initial state
    ψ = zeros(ComplexF64,2^(num_q))
    ψ[1] = 1

    HaarRandomCircuitAncilla.full_evolution!(ψ,T,p,num_q)
    entanglePBC!(ψ,num_q)
    HaarRandomCircuitAncilla.full_evolution_ancilla!(ψ,8T,p,num_q,s)  
    return s
end


function multiple_runs(p::Float64,L::Int64,iter::Int64)
    T = 2L
    MI_t = zeros(Float64,8T+1,4)

    for r in 1:iter
        MI_t += main(p,L,r)
        println(r)
    end
    MI_t ./= iter
    df = DataFrame(MI_t,:auto)
    CSV.write("spatial_correlation_$(p)_$(L)_dolly.csv",df,header=["s1", "s2", "s3","sinf"])
    
end
# pwd()
pvalues = [0.0,0.05,0.1,0.12,0.14,0.16,0.18,0.20,0.22,0.24]
for p in pvalues
    println("p = ",p)
    multiple_runs(p,8,1200)
end
main(0.258,1)
# Checking Product of 2 Bell states
ψ = zeros(ComplexF64, 2^4)

ψ[1],ψ[1+5],ψ[1+10],ψ[1+15] = 0.5exp(im*2),0.5exp(im*2),0.5exp(im*2),0.5exp(im*2)
HaarRandomCircuitAncilla.spatial_correlation(ψ)


# Checking a Bell state
ψ = zeros(ComplexF64, 2^4)

ψ[1],ψ[1+15] = sqrt(0.5),sqrt(0.5)
HaarRandomCircuitAncilla.spatial_correlation(ψ)

#Checking only 1 ancilla entangled to the system
# Ancilla 2

ψ = zeros(ComplexF64, 2^4)
ψ[1],ψ[7] = 1/sqrt(2),1/sqrt(2)
HaarRandomCircuitAncilla.spatial_correlation(ψ)

# Ancilla 1

ψ = zeros(ComplexF64, 2^4)
ψ[1],ψ[11] = 1/sqrt(2),1/sqrt(2)
HaarRandomCircuitAncilla.spatial_correlation(ψ)
