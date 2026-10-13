# Active learning Bi--Sb dislocation


This repository implements an **automated force-following active learning workflow** for developing a **Gaussian Approximation Potential (GAP)** for **dislocations in the Bi–Sb topological insulator system**.  
The workflow couples **LAMMPS**, **QUIP/GAP**, and **Density Functional Theory (DFT)** calculations using **Quantum ESPRESSO**, and is designed for execution on HPC systems using **SLURM, MPI, and OpenMP**.

---

## Overview

The goal of this workflow is to efficiently construct a **dislocation-aware interatomic potential** by iteratively identifying and correcting regions of high model uncertainty.  
Rather than uniformly sampling configuration space, the method focuses DFT calculations on atomic environments where the GAP model exhibits large force variance.

At each iteration, the workflow:
1. Minimizes atomic configurations using the current GAP model.
2. Evaluates GAP uncertainty on the minimized structure.
3. Identifies the atomic environment with the maximum uncertainty.
4. Performs a DFT calculation for that environment.
5. Updates the GAP training database and refits the potential.
6. Repeats until uncertainty falls below a prescribed threshold.

---

## GAP Potential

The GAP potential files used for the compression calculations are archived on Zenodo:

- [GAP potential files (Zenodo record 23346774)](https://zenodo.org/records/23346774)

Download the XML file together with all associated `sparseX` files.

---

## Main Script

- **`main.sh`**
  - SLURM job script that controls the full active learning loop
  - Handles evaluation, minimization, DFT execution, and GAP retraining
  - Designed for multi-node CPU execution

---

## Software Requirements

The following software must be available on the target system:

- **LAMMPS** (compiled with QUIP support)
- **QUIP / GAP**
- **Quantum ESPRESSO** (v7.x or compatible)
- **Python 3**
- **MPI (OpenMPI)**
- **BLAS/LAPACK or Intel MKL**

The script assumes that required Python utilities are available in the active environment.

---

## HPC and Parallelization

The workflow is configured for large-scale parallel execution:

- **Parallel components**
  - LAMMPS: MPI
  - QUIP/GAP fitting: MPI + OpenMP
  - Quantum ESPRESSO: MPI


