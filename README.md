# polarization-ionization-gating-hhg

Fortran implementation of Lewenstein strong-field approximation model to simulate
high-harmonic generation (HHG) and isolated attosecond pulse generation for He atom
driven by orthogonally polarized combined fields.  (Demo version)

## Related publication
Yu, W. et al., "Attosecond pulse generation isolated with a polarization-ionization gating scheme",
Eur. Phys. J. D 73, 236 (2019).

## Quick start (demo)
### Requirements
- Fortran compiler: gfortran (recommended) or ifort
- OS: Linux/macOS/Windows (with compiler)

### Compile
```bash
gfortran -O1 src/<main_file>.f90 -o hhg_demo
gfortran -O2 src/<main_file_2>.f90 -o isolated_attosecond_pulse_demo
