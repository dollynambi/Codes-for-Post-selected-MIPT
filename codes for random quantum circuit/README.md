
# How to run the codes for the post-selected random quantum circuits ?

**Environment**

The code was developed and tested with:

- Julia version: 1.11.1

Required packages

- LinearAlgebra
- Random 
- Combinatorics
- CSV
- DataFrames

This folder contains sub-folders corresponding to computing different entanglement measures for computing various critical exponents. The sub folders are as follows:

## `Bulk exponent`-  local order parameter exponent $\beta$
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

**Arguments:**
- `Arg1` (`p`, Float64): Measurement rate within 0 to 1.
- `Arg2` (`L`, Int): Number of qubits (system size).
- `Arg3` (`r`, Int): Random seed for reproducibility.

### Output files
The code generates a CSV file of the form:

`bulk_ancilla_$(p)_$(L)_$(r).csv`
where:
- `$(p)` = measurement rate  
- `$(L)` = system size (number of qubits)  
- `$(r)` = random seed

The headers for these files are s1 (von Neumann entropy of the ancilla) , s2 (2nd order Renyi entropy of the ancilla), s3 (3nd order Renyi entropy of the ancilla) and sinf (infinite-order (min) Renyi entropy).


## `tmi_p` - Critical point $p_c$ and correlation length exponent $\nu$ of the entanglement transition
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

## `tmi_t` - Computing the coeffecients $\alpha_n$ ($S_n \sim \alpha_n log L$) of the Renyi entropy $S_n$ 
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

## `ancilla_p` - Critical point $p_c$ and correlation length exponent $\nu$ of the purification transition
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


## `ancilla_t`- Computing the dynamical exponent $z$ of the purification transition
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

## spatial_correlation - Computing the spatial correlation length exponent, $\eta$
This folder contains the files:
- `HaarRandomCircuitAncilla.jl` : Has the functions to run the quantum circuit entangled with two ancillas.
- `spatial_correlation.jl` : Runs the circuit using the functions in `HaarRandomCircuitAncilla.jl` where the ancillas entangled to bulk of the circuit to two different qubits at the same time.

### Running the code
To compute the spatial correlation all the files in `spatial_correlation` should be in the same directory.

```bash
julia spatial_correlation.jl Arg1 Arg2 Arg3
```  
**Arguments:**
- `Arg1` (`p`, Float64): Measurement rate within 0 to 1.
- `Arg2` (`L`, Int): Number of qubits (system size).
- `Arg3` (`seed`, Int): Random seed for reproducibility.

### Output files

The `spatial_correlation.jl` generates a CSV file of the form:

`spatial_correlation_$(p)_$(L)_$(r).csv`

where:
- `$(p)` = measurement rate  
- `$(L)` = system size (number of qubits)  
- `$(r)` = random seed

The headers are s1 (von Neumann entropy of the ancilla) , s2 (2nd order Renyi entropy of the ancilla), s3 (3nd order Renyi entropy of the ancilla) and sinf (infinite-order (min) Renyi entropy).

## `Free energy` - Computing the effective central charge $c_{eff}$ and the anisotropy factor $v$
This folder contains the files:
- `HaarRandomCircuit.jl` : Has the functions to run the quantum circuit.
- `FreeEnergy.jl` : Runs the circuit using the functions in  `HaarRandomCircuit.jl` and computes the free energy.
- `HaarRandomCircuitAncilla.jl` : Has the functions to run the quantum circuit entangled with two ancillas.
- `temporal_correlation.jl` : Runs the circuit using the functions in `HaarRandomCircuitAncilla.jl` where the ancillas entangled to bulk of the circuit at two different times.

### Running the code 
- To compute the free energy, both `HaarRandomCircuit.jl` and `FreeEnergy.jl` should be in the same directory, then it can be run by 
```bash
julia FreeEnergy.jl Arg1 Arg2 Arg3
```  
- To compute the temporal correlation, both `HaarRandomCircuitAncilla.jl` and `temporal_correlation.jl` should be in the sam directory, then it can be run by

```bash
julia temporal_correlation.jl Arg1 Arg2 Arg3
```  
**Arguments:**
- `Arg1` (`p`, Float64): Measurement rate within 0 to 1.
- `Arg2` (`L`, Int): Number of qubits (system size).
- `Arg3` (`seed`, Int): Random seed for reproducibility.

### Output files
The `FreeEnergy.jl` generates a CSV file of the form:

`logp_$(p)_$(L)_$(r).csv`

The `temporal_correlation.jl` generates a CSV file of the form:

`temporal_correlation_$(p)_$(L)_$(r).csv`

where:
- `$(p)` = measurement rate  
- `$(L)` = system size (number of qubits)  
- `$(r)` = random seed

The headers for `temporal_correlation_$(p)_$(L)_$(r).csv` are s1 (von Neumann entropy of the ancilla) , s2 (2nd order Renyi entropy of the ancilla), s3 (3nd order Renyi entropy of the ancilla) and sinf (infinite-order (min) Renyi entropy).






 
