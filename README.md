# polarization-ionization-gating-hhg

Fortran implementation of Lewenstein strong-field approximation model to simulate
high-harmonic generation (HHG) and isolated attosecond pulse generation for He atom
driven by orthogonally polarized combined fields.  (Demo version)

# polarization-ionization-gating-hhg

Two-step Fortran workflow:
1) **Program 1 (HHG)** generates high-harmonic spectrum data.
2) **Program 2 (IAP)** synthesizes an isolated attosecond pulse by superposing selected harmonics.

Related publication:
Yu, W. et al., *Attosecond pulse generation isolated with a polarization-ionization gating scheme*,
Eur. Phys. J. D 73, 236 (2019).

---

## Quick start

### Requirements
- Fortran compiler: `gfortran` (recommended) or Intel `ifort`

---

### Step 1 — HHG spectrum generation (Program 1)
**Source:** `src/source1.f90`

Compile:
```bash
gfortran -O2 src/source1.f90 -o hhg_demo
Run:

./hhg_demo


How to change parameters
This project does not use a separate input file. To modify simulation parameters, edit them directly in src/source1.f90, then recompile and rerun.

Output:

HHG spectrum data files are written to the output folder defined in the code.

Example outputs: output_example/
