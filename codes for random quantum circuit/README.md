
# How to run the codes for the post-selected random quantum circuits

This folder contains sub-folders corresponding to computing different entanglement measures for computing various critical exponents. The sub folders are as follows:

## `Bulk exponent`
This folder contains the files:
 -  `HaarRandomCircuitAncilla.jl` :
     Contains all the functions relevant for running the circuit and computing the ancilla entropy.
 -  `bulk_critical_exponent.jl` :
     Code to run the simulation and compute the local order parameter (entropy of the ancilla entangled into the bulk) for each time step till t = 16L for different system sizes and measurement rates.
    
### Running the code
To run the code all the files in `Bulk exponent` should be in the same directory, then it can be run by 
```bash
julia bulk_critical_exponent.jl Arg1 Arg2 Arg3
```
## `Free energy`

## `tmi_p`
This folder contains the files:
 -  `HaarRandomCircuit.jl` :
     Contains all the functions relevant for running the circuit and computing the half-cut entropy and the tripartite mutual information.
 -  `mutualInformationHaar.jl` :
     Script to run the simulation and compute the half-cut entropy and the tripartite mutual information at time t = 2L,  and t =4L for different system sizes and measurement rates.
    
### Running the code
To run the code all the files in `tmi_p` should be in the same directory, then it can be run by 
```bash
julia mutualInformationHaar.jl Arg1 Arg2 Arg3
```
**Arguments:**
- `Arg1` (`p`, Float64): Measurement rate within 0 to 1.
- `Arg2` (`L`, Int): Number of qubits (system size).
- `Arg3` (`r`, Int): Random seed for reproducibility.

### Output files
The code generates CSV files of the form:

- `tmi_$(p)_$(L)_$(r)_<tag>.csv` - Tripartite mutual information
- `sab_$(p)_$(L)_$(r)_<tag>.csv` - Half-cut entropy
where:
- `$(p)` = measurement rate  
- `$(L)` = system size (number of qubits)  
- `$(r)` = random seed  
The `<tag>` indicates the time at which the ancilla entropy is computed:

- `L`     → \(t = 2L\)  
- `2L`    → \(t = 4L\)

The headers for these files are  S0 (Hartley entropy), S1 (von Neumann entropy) , S2 (2nd order Renyi entropy), S3 (3nd order Renyi entropy) and Sinf (infinite-order (min) Renyi entropy).

## `tmi_t`
This folder contains the files:
 -  `HaarRandomCircuit.jl` :
     Contains all the functions relevant for running the circuit and computing the half-cut entropy and the tripartite mutual information.
 -  `mutualInformationHaar.jl` :
     Code to run the simulation and compute the entropies of the different partitions and the tripartite mutual information at time step till t = 4L for different system sizes and measurement rates.
    
### Running the code
To run the code all the files in `tmi_t` should be in the same directory, then it can be run by 
```bash
julia mutualInformationHaar.jl Arg1 Arg2 Arg3
```
**Arguments:**
- `Arg1` (`p`, Float64): Measurement rate within 0 to 1.
- `Arg2` (`L`, Int): Number of qubits (system size).
- `Arg3` (`r`, Int): Random seed for reproducibility.

### Output files
The code generates CSV files of the form:

- `tmi_$(p)_$(L)_$(r).csv` - Tripartite mutual information
- `s<tag>_$(p)_$(L)_$(r).csv` - von Neumann entropy of the density matrix corresponding part `<tag>`
where:
- `$(p)` = measurement rate  
- `$(L)` = system size (number of qubits)  
- `$(r)` = random seed  
The `<tag>` indicates the partition

- `a`     → A and the rest (BCD)  
- `b`    → B and the rest (ACD)
- `c`    → C and the rest (ABD)
- `abc`    → D and the rest (ABC)
- `ab`    → AB and the rest (CD)
- `bc`    → BC and the rest (AD)
- `ac`    → AC and the rest (BD)

The headers for these files are  S0 (Hartley entropy), S1 (von Neumann entropy) , S2 (2nd order Renyi entropy), S3 (3nd order Renyi entropy) and Sinf (infinite-order (min) Renyi entropy).

## `ancilla_p`
This folder contains the files:
 -  `HaarRandomCircuitAncilla.jl` :
     Contains all the functions relevant for running the circuit and computing the ancilla entropy.
 -  `AncillaEntropy.jl` :
     Script to run the simulation and compute the ancilla entropy at time t = L, t = 2L, t = 3L and t =4L for different system sizes and measurement rates.
    
### Running the code
To run the code all the files in `ancilla_p` should be in the same directory, then it can be run by 
```bash
julia AncillaEntropy.jl Arg1 Arg2 Arg3
```
**Arguments:**
- `Arg1` (`p`, Float64): Measurement rate within 0 to 1.
- `Arg2` (`L`, Int): Number of qubits (system size).
- `Arg3` (`r`, Int): Random seed for reproducibility.

### Output files
The code generates CSV files of the form:

`ancilla_entropy(<tag>)_$(p)_$(L)_$(r).csv`
where:
- `$(p)` = measurement rate  
- `$(L)` = system size (number of qubits)  
- `$(r)` = random seed  
The `<tag>` indicates the time at which the ancilla entropy is computed:

- `L_2`   → \(t = L\)  
- `L`     → \(t = 2L\)  
- `3L_2`  → \(t = 3L\)  
- `2L`    → \(t = 4L\)

The headers for these files are s1 (von Neumann entropy of the ancilla) , s2 (2nd order Renyi entropy of the ancilla), s3 (3nd order Renyi entropy of the ancilla) and sinf (infinite-order (min) Renyi entropy) 


## `ancilla_t`
This folder contains the files:
 -  `HaarRandomCircuitAncilla.jl` :
     Contains all the functions relevant for running the circuit and computing the ancilla entropy.
 -  `AncillaEntropy.jl` :
     Script to run the simulation and compute the ancilla entropy each time step for different system sizes and measurement rates.
    
### Running the code
To run the code all the files in `ancilla_p` should be in the same directory, then it can be run by 
```bash
julia AncillaEntropy.jl Arg1 Arg2 Arg3
```
**Arguments:**
- `Arg1` (`p`, Float64): Measurement rate within 0 to 1.
- `Arg2` (`L`, Int): Number of qubits (system size).
- `Arg3` (`seed`, Int): Random seed for reproducibility.

### Output files
The code generates a CSV file of the form:

`ancilla_entropy_$(p)_$(L)_$(r).csv`
where:
- `$(p)` = measurement rate  
- `$(L)` = system size (number of qubits)  
- `$(r)` = random seed

The headers for these files are s1 (von Neumann entropy of the ancilla) , s2 (2nd order Renyi entropy of the ancilla), s3 (3nd order Renyi entropy of the ancilla) and sinf (infinite-order (min) Renyi entropy).





 
