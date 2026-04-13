include("HaarRandomCircuitAncilla.jl")
using .HaarRandomCircuitAncilla
using LinearAlgebra
using CSV, DataFrames, Random

function main(p::Float64,L::Int64,r::Int64,t0::Int64)
    Random.seed!(r)
    s = run_periodic(p,L,t0)
    df = DataFrame(s,:auto)
    CSV.write("temporal_correlation_$(p)_$(L)_$(r).csv",df,header=["s1", "s2", "s3","sinf"] )
end

function entangle_ancilla2!(ψ::Vector{ComplexF64},L::Int64)
    ψmat = reshape(ψ,2^L,4)
    i = 1
    dim = length(ψmat[:,1])
    dimL = dim ÷ 2^i 
    dimR = 2^(i-1)

    for nL=0:(dimL-1),nR=(0:dimR-1)
        n0 = nR + 2^(i)*nL + 1
        n1 = n0 + 2^(i-1)
        
        ψmat[n1,1], ψmat[n0,2], ψmat[n1,3], ψmat[n0,4] = 0,0,0,0
        ψmat[n1,2] = ψmat[n0,1]
        ψmat[n1,4] = ψmat[n0,3]
    end
    
    ψ .= reshape(ψmat,2^(L+2),1)
    normalize!(ψ)   
end


function entangle_ancilla1!(ψ::Vector{ComplexF64},L::Int64)
    ψmat = reshape(ψ,(2^(L+1),2))
    i = 1
    dim = length(ψmat[:,1])
    dimL = dim ÷ 2^i 
    dimR = 2^(i-1)

    for nL=0:(dimL-1),nR=(0:dimR-1)
        n0 = nR + 2^(i)*nL + 1
        n1 = n0 + 2^(i-1)
        
        ψmat[n1,1], ψmat[n0,2] = 0,0
        ψmat[n1,2] = ψmat[n0,1]
    end
    
    ψ .= reshape(ψmat,(2^(L+2),1))
 
    normalize!(ψ)   
end

function run_periodic(p::Float64,L::Int64,t0::Int64)
    
    T = 2L
    s = zeros(Float64,8T+1,4)
    num_q = L+2
    # Defining the initial state
    ψ = zeros(ComplexF64,2^(num_q))
    ψ[1] = 1

    HaarRandomCircuitAncilla.full_evolution!(ψ,T,p,num_q)
    entangle_ancilla1!(ψ,L)
    HaarRandomCircuitAncilla.full_evolution!(ψ,t0,p,num_q)
    U = Array{ComplexF64}(undef, 4, 4)
    temp_vec1 = Array{ComplexF64}(undef, 4)
    temp_vec2 = Array{ComplexF64}(undef, 4)
    HaarRandomCircuitAncilla.layer_ancilla!(ψ, U, temp_vec1, temp_vec2, num_q; even_layer = true)
    HaarRandomCircuitAncilla.post_select_measure_set_0!(ψ, 1:L, p)
    entangle_ancilla2!(ψ,L) 
    HaarRandomCircuitAncilla.full_evolution_ancilla!(ψ,4,p,num_q,s;start_layer_even = false)  
    return s
end



main(p,L,r,t0)