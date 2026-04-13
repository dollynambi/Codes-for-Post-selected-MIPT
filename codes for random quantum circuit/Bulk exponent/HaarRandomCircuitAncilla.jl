module HaarRandomCircuitAncilla

using Combinatorics, LinearAlgebra, Random


function randomHaar!(U::Array{ComplexF64,2})

    dim, dim = size(U)
    randn!(U)
    Q, R = qr(U)
    L = broadcast( x -> abs(x)>0 ? sign(x) : 0, Diagonal(R))
    U .= Matrix(Q)*L
end


function randomHaar!(U::Array{ComplexF64,2})
    # Look into how to make Haar from QR decomposition
    dim, dim = size(U)
    randn!(U)
    Q, R = qr(U)
    L = broadcast( x -> abs(x)>0 ? sign(x) : 0, Diagonal(R))
    U .= Matrix(Q)*L
end

function gate!(ψ::Vector{ComplexF64}, U::Matrix{ComplexF64}, i::Integer, j::Integer, 
               temp_vec1::Vector{ComplexF64}, temp_vec2::Vector{ComplexF64})

    i, j = i < j ? (i, j) : (j, i)
    randomHaar!(U)

    ## States have the form (nL, zi, nM, zj, nR) 
    ## where nL has i-1 qubits in it so dim = 2^(i - 1)
    ##       nM has j-i-1 qubits in it so dim = 2^(j - i - 1)
    ##       nR has L-j-1 qubits in it so dim = 2^(L - j)
    ## Example: L = 20, i = 5, j = 13
    ## [1 2 3 4] 5 [6 7 8 9 10 11 12] 13 [14 15 16 17 18 19 20]
    ## ----4---- X ---------7-------- X  ----------7-----------
    ## 2^(i - 1) = 2^4
    ## 2^(j - i - 1) = 2^7
    ## 2^(20 - j) = 2^7
    ## To represent states as integers we use sum( b[n]*2^(n-1) ) + 1 for b[n] = 0,1
    ## In this case we have:
    ## 1 + nL + b[i]*2^(i-1) + nM * 2^i + b[j]*2^(j - 1) + nR * 2^j
    ## nR = 0:2^(L - j) - 1
    ## nM = 0:2^(j - i - 1) - 1
    ## nL = 0:2^(i - 1) - 1
    ## Test extremes:
    ## (*) 1 + (2^(i-1) - 1) + 2^(i - 1) + (2^(j - 1) - 2^i) + 2^(j-1) + 2^L - 2^j
    ## 1 + (2^(i-1) - 1) + 2^(i - 1) = 2^i  
    ## (2^(j - 1) - 2^i) + 2^(j - 1) = 2^j - 2^i
    ## => (*) is equal to 2^L 

    ibit = 2^(i-1)
    jbit = 2^(j-1)
    ijbit = ibit+jbit

    dim = length(ψ)
    dimL = 2^(i - 1)
    dimM = 2^(j - i - 1)
    dimR = dim ÷ 2^j

    for nR = 0:dimR - 1, nM = 0:dimM - 1, nL = 0:dimL - 1
        ind1 = 1 + nL + nM * 2^i + nR * 2^j
        ind2 = ind1 + ibit 
        ind3 = ind1 + jbit
        ind4 = ind1 + ijbit 

        temp_vec1[1] = ψ[ind1]
        temp_vec1[2] = ψ[ind2]
        temp_vec1[3] = ψ[ind3]
        temp_vec1[4] = ψ[ind4]

        mul!(temp_vec2, U, temp_vec1)

        ψ[ind1] = temp_vec2[1]
        ψ[ind2] = temp_vec2[2]
        ψ[ind3] = temp_vec2[3]
        ψ[ind4] = temp_vec2[4]

    end

end

function layer_ancilla!(ψ::Vector{ComplexF64}, U::Matrix{ComplexF64}, temp_vec1::Vector{ComplexF64}, 
    temp_vec2::Vector{ComplexF64}, N::Integer; even_layer::Bool = true)

    L = N - 1 
    for i = (2 - even_layer):2:(L - 1)
        gate!(ψ, U, i, i + 1, temp_vec1, temp_vec2)
    end
    if ~even_layer
        gate!(ψ, U, 1, L, temp_vec1, temp_vec2)
    end

end

# Always Projecting the qubit to the |0> state
function post_select_measure_0!(ψ::Vector{ComplexF64},i::Integer)
    dim = length(ψ)
    dimL = dim ÷ 2^i 
    dimR = 2^(i-1)

    for nL=0:(dimL-1), nR=(0:dimR-1)
        n = nR+2^(i)*nL+2^(i-1)+1
        ψ[n] = 0.0
    end
    normalize!(ψ)
end

# Always Projecting the qubit to the |1> state
function post_select_measure_1!(ψ::Vector{ComplexF64},i::Integer)

    dim = length(ψ)
    dimL = dim ÷ 2^i 
    dimR = 2^(i-1)

    for nL=0:(dimL-1), nR=(0:dimR-1)
        n = nR+2^(i)*nL +1
        ψ[n] = 0.0
    end 
    normalize!(ψ)
end

function post_select_measure_set_0!(ψ::Vector{ComplexF64},sites::AbstractVector,p::Float64; use_ancilla::Bool = false)
    for i=sites
        if rand()<p
            post_select_measure_0!(ψ,i)
        end
    end
end

function brick_layer_ancilla!(ψ::Vector{ComplexF64}, U::Matrix{ComplexF64}, temp_vec1::Vector{ComplexF64}, 
    temp_vec2::Vector{ComplexF64}, p::Float64, N::Integer; start_layer_even::Bool = true)

    L = N - 1
    layer_ancilla!(ψ, U, temp_vec1, temp_vec2, N; even_layer = start_layer_even)
    post_select_measure_set_0!(ψ, 1:L, p)
    layer_ancilla!(ψ, U, temp_vec1, temp_vec2, N; even_layer = ~start_layer_even)
    post_select_measure_set_0!(ψ, 1:L, p)
end

function full_evolution!(ψ::Vector{ComplexF64}, T::Integer, p::Float64, N::Integer; start_layer_even::Bool = true)
    @assert 2^N == length(ψ) "Must have N qubits in ψ"
    U = Array{ComplexF64}(undef, 4, 4)
    temp_vec1 = Array{ComplexF64}(undef, 4)
    temp_vec2 = Array{ComplexF64}(undef, 4)

    for t = 1:T
        brick_layer_ancilla!(ψ, U, temp_vec1, temp_vec2, p, N; start_layer_even = start_layer_even)
    end

end

function full_evolution_ancilla!(ψ::Vector{ComplexF64}, T::Integer, p::Float64, N::Integer,s::Array{Float64}; start_layer_even::Bool = true)
    @assert 2^N == length(ψ) "Must have N qubits in ψ"
    U = Array{ComplexF64}(undef, 4, 4)
    temp_vec1 = Array{ComplexF64}(undef, 4)
    temp_vec2 = Array{ComplexF64}(undef, 4)
    s[1,:] = ancilla_entropy(ψ)

    for t = 1:T
        brick_layer_ancilla!(ψ, U, temp_vec1, temp_vec2, p, N; start_layer_even = start_layer_even)
        s[t+1,:] = ancilla_entropy(ψ)
    end

end

function ancilla_entropy(ψ::Vector{ComplexF64})
    L = Int(log2(length(ψ))) - 1
    rho = zeros(ComplexF64,2,2)
    ψmat =  reshape(ψ,(2^L,2))
    for i in 1:2, j in 1:2
        rho[i,j] = dot((ψmat[:,i]),ψmat[:,j])
    end

    lambda = eigvals(rho)

    s1,s2,s3 = 0,0,0
    
    for i in lambda
        p = abs(i)
        if p != 0
            s1 = s1 - p*log2(p)
        end
        s2 = s2 + p^2
        s3 = s3 + p^3 
    end
    sinf = log2(1/maximum(abs.(lambda)))
    return [s1,-log2(s2),-log2(s3)/2,sinf] 
end



end

