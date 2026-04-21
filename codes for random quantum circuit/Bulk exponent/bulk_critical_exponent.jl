include("HaarRandomCircuitAncilla.jl")
using .HaarRandomCircuitAncilla
using LinearAlgebra
using CSV, DataFrames, Random

function main(p::Float64,L::Int64,r::Int64)
    Random.seed!(r)
    s = run_periodic(p,L)
    df = DataFrame(s,:auto)
    CSV.write("bulk_ancilla_$(p)_$(L)_$(r).csv",df,header=["s1", "s2", "s3","sinf"] )
end


function entangle_ancilla!(ψ::Vector{ComplexF64},L::Int64)
    ψmat = reshape(ψ,(2^(L),2))
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
    
    ψ .= reshape(ψmat,(2^(L+1),1))
    normalize!(ψ)   
end


function run_periodic(p::Float64,L::Int64)
    
    T = 2L
    s = zeros(Float64,8T+1,4)
    num_q = L+1
    # Defining the initial state
    ψ = zeros(ComplexF64,2^(num_q))
    ψ[1] = 1

    HaarRandomCircuitAncilla.full_evolution!(ψ,T,p,num_q)
    entangle_ancilla!(ψ,L)
    HaarRandomCircuitAncilla.full_evolution_ancilla!(ψ,8T,p,num_q,s)  
    return s
end
p = parse(Float64,ARGS[1])
L = parse(Int,ARGS[2])
r = parse(Int,ARGS[3])
main(p,L,r)
