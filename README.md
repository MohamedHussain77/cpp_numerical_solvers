# CFD Solver Development in C++

A from-scratch CFD solver project written in C++, developed as a hands-on journey toward **scientific computing, HPC, parallel programming, GPU computing, and OpenFOAM solver development**.

The goal is not to reproduce a production CFD package, but to understand the fundamental building blocks behind one by implementing them progressively from scratch.

## Current Status

The project currently contains:

* 1D linear advection solver
* 1D diffusion solver
* 1D Euler solver — **in development**
* Structured mesh representation
* Field-based data structures
* Initial conditions and boundary conditions
* Finite-difference / finite-volume numerical methods
* Explicit time integration
* RK2 time stepping

## Project Structure

```text
.
├── common/
│   └── Reusable CFD infrastructure
│
├── advection1d/
│   └── 1D linear advection equation
│
├── diffusion1d/
│   └── 1D diffusion equation
│
└── euler1d/
    └── 1D compressible Euler equations
```

The separation is intentional: physics-specific components live inside their respective solver directories, while genuinely reusable infrastructure is placed in `common/`.

## Numerical Problems

### 1D Linear Advection

$$
\frac{\partial u}{\partial t}
+
a\frac{\partial u}{\partial x}=0
$$

Implemented numerical methods include:

* First-order upwind
* Lax-Friedrichs
* Lax-Wendroff

### 1D Diffusion

The diffusion equation is used to develop the basic spatial discretization and time-integration framework.

### 1D Euler — In Progress

The next major step is a compressible-flow solver based on the conservative formulation:

$$
\frac{\partial \mathbf U}{\partial t}
+
\frac{\partial \mathbf F}{\partial x}=0
$$

with

$$
\mathbf U =
\begin{bmatrix}
\rho\\
\rho u\\
\rho E
\end{bmatrix}.
$$

The Euler solver will introduce a multi-variable state representation and eventually be extended to 2D.

## Learning Objectives

This project is being developed to build practical understanding of:

* Modern C++ for scientific computing
* Object-oriented and data-oriented solver design
* Mesh and field abstractions
* Numerical discretization
* Stability and CFL conditions
* Time integration
* Conservation laws
* Compressible-flow solvers
* Parallel CFD
* HPC programming
* GPU programming

## Roadmap

```text
1D Advection
     ↓
1D Diffusion
     ↓
1D Euler
     ↓
2D Euler
     ↓
OpenMP
     ↓
MPI / Domain Decomposition
     ↓
MPI + OpenMP
     ↓
GPU / CUDA
     ↓
OpenFOAM Source-Code Development
     ↓
High-Performance CFD Solvers
```

The roadmap will evolve as the implementations become more sophisticated.

## Philosophy

The main objective is to understand **how a CFD solver is built**, rather than simply using an existing CFD package.

Each new feature is introduced only after understanding its numerical and computational role. The project therefore progresses from simple scalar equations toward coupled compressible-flow equations and eventually toward parallel and heterogeneous computing.

## Tools

* C++
* CMake / compiler toolchain
* Python for post-processing and visualization
* Git / GitHub
* OpenFOAM — planned for source-code study and solver development

## Future Goals

The long-term goal is to develop the skills required to work on:

**CFD algorithms + scientific C++ + HPC + parallel computing + GPU computing**

with particular emphasis on understanding and developing high-performance CFD solvers.
