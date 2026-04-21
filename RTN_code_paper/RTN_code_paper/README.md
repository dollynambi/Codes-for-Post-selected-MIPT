# RTN Code Paper

This folder contains the standalone Python scripts used for extracting critical exponents for the entanglement transition in random tensor newtorks:

- `rtn_circ_beta.py`
- `rtn_circ_dyn_exp.py`
- `rtn_circ_ee.py`
- `rtn_circ_free_energy.py`
- `rtn_circ_mi.py`
- `rtn_circuit_eta.py`
- `rtn_circuit_style_ancilla_op.py`
- `rtn_circuit_style_tmi.py`
- `rtn_uniform_circ_style_tmi.py`

## Environment

The dependency versions used for this folder are recorded in `pyproject.toml`.

Core dependencies:

- `numpy==1.26.4`
- `scipy==1.11.4`
- `matplotlib==3.8.4`

Python version:

- `3.12.3`

## Running

Each script is standalone and can be run directly, for example:

```bash
python3 rtn_circ_ee.py --help
```
