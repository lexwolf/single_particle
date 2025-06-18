# Quasi-Static Metal Nanoparticle Simulation

This repository contains a minimal and self-contained implementation of the quasi-static model for a metal nanoparticle embedded in an infinite gain medium, as described in the open-access publication:

**Reference**: ["Multipolar, time-dynamical model for the loss compensation and lasing of a spherical plasmonic nanoparticle spaser immersed in an active gain medium"](https://www.nature.com/articles/srep33018), *Scientific Reports* (2016).

---

## Overview

The code models the optical response of a spherical nanoparticle in the quasi-static limit, using a time-domain and frequency-domain approach. The implementation supports:

- Analytical steady-state polarizability calculation
- Full time-dependent emission dynamics in the weak-excitation regime
- Material parametrization via Lorentz gain and Drude/spline metal models

The code is kept minimal and does not rely on Mie theory or nanoshell extensions. This is the core starting point for more advanced spaser models.

---

## Compilation

### Requirements

- C++ compiler (e.g., `g++`)
- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
- [Armadillo](http://arma.sourceforge.net/)

### Steps

1. Run the `configure` script to verify library dependencies:
   ```bash
   ./configure
   ```

2. Compile the project:
   ```bash
   make
   ```

3. Clean compiled binaries:
   ```bash
   make clean
   ```

---

## Input Files

### `data/input/nanosphere_eV.dat`

| Parameter               | Description                                                          | Example |
|-------------------------|----------------------------------------------------------------------|---------|
| Radius (nm)             | Nanoparticle radius                                                  | 10      |
| Emission width (eV)     | \(\hbar \Delta\), gain linewidth                                     | 0.15    |
| Emission frequency (eV) | \(\hbar \omega_0\), gain emission frequency                          | 2.81    |
| Gain level              | \(G\), complex permittivity amplification factor                     | 0.033   |
| Spectrum range (eV)     | \(\omega_{\text{min}}, \omega_{\text{max}}\)                         | 2.0–4.5 |
| Metal type              | e.g., `silver`, `gold`                                               | silver  |
| Metal model             | `drude` or `spline`                                                  | drude   |
| Gain model              | `lorentz`, `flat`, etc.                                              | lorentz |
| Host medium             | e.g., `water`                                                        | water   |
| Excitation amplitude    | \(E_0\) (arbitrary units)                                            | 1.e-30  |

### `data/input/time.dat`

| Parameter           | Description                         | Example |
|---------------------|-------------------------------------|---------|
| Total time (ps)     | Duration of time-domain simulation  | 100     |
| Pump-on time (ps)   | Time when gain becomes active       | 10      |

---

## Executables

| Executable | Source File       | Description                                                   |
|------------|-------------------|---------------------------------------------------------------|
| `anl`      | `anl_time.cxx`    | Computes analytical polarizability over time and frequency    |
| `sgl`      | `single.cxx`      | Full simulation with steady-state and dynamic results         |
| `num`      | `num_time.cxx`    | Numerical integration of time-dependent dipole amplitude      |
| `fro`      | `frohlich.cxx`    | Computes Frohlich condition and response spectrum             |

Each binary reads input from `data/input/` and writes output to `data/output/` or standard output.

---

## License

This project is released under the GNU GENERAL PUBLIC LICENSE. See the `LICENSE` file for details.

---

## Contact

For questions or collaboration, contact:

**Alessandro Veltri**  
📧 alessandro.veltri@gmail.com  
