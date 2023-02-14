# ITVOLT: An Iterative Solver for the Time-Dependent Schrödinger Equation
In atomic units, the time-dependent Schrödinger equation (TDSE) takes the form

$$ i \frac{\partial}{\partial t} \psi (t) = {\bf H}(t) \psi(t) $$

for a Hamiltonian ${\bf H}$ and a corresponding wave function $\psi$. ITVOLT (short for Iterative Volterra Propagator) is a novel method for solving the TDSE that propagates a solution on intervals $[\tau_j, \tau_{j+1}]$ by solving the equivalent Volterra integral equation

$$ \psi(t) = e^{-i {\bf H}_j(t - \tau_j)} \psi(\tau_j) - i \int_{\tau_j}^t e^{-i {\bf H}_j(t-t')} {\bf V}_j(t') \psi(t')dt', $$ 

where ${\bf H}_j = {\bf H}_0 + {\bf W}(\tau_j + \frac{\Delta \tau}{2})$ and ${\bf V}_j(t) = {\bf W}(t) - {\bf W}\left( \tau_j + \frac{\Delta \tau}{2} \right)$ for ${\bf H}_0$ and ${\bf W}(t)$ the time-independent and time-dependent parts of ${\bf H}$ respectively and $\Delta \tau$ the size of the interval. ITVOLT proceeds by applying a Lagrange interpolation to the integrand, reducing the Volterra integral equation to a linear system that is then solved iteratively. The numerical details of the method are presented in a [forthcoming paper.](https://arxiv.org/abs/2210.15677)

This repository contains Fortran 90 code that applies ITVOLT to the TDSE for the following quantum systems:
 * A two level atom exposed to a laser.
 * A driven harmonic oscillator.
 
Both of these examples demonstrate how the repository can be used to solve the TDSE for any problem satisfying the following criteria:
1. The Hamiltonian can be represented by a symmetric banded matrix.
2. The time dependent piece of the Hamiltonian can be written as $E(t)*{\bf V}$ for some function of time $E$ and some fixed, symmetric banded matrix ${\bf V}$.

Because ITVOLT applies more broadly to any problem that can be written as a Volterra integral equation, the repository contains additional code to solve two non-physical problems:
 * A simple two channel problem of Wang and Wang [1].
 * The following one dimensional ODE (with initial condition $\psi(0) = 1$):
 
 $$ \left[ i \frac{ \partial}{\partial t} - t \right] \psi(t) = 0. $$

## Usage
After downloading the repository, begin by creating a build directory within Iterative_Volterra_Propagator.
```
mkdir Build
```
Next, run CMAKE to create a Makefile in Build.
```
cmake -B Build
```
You can now move to Build and compile the code.
```
cd Build
make
```
This creates a separate executable in Source (the one inside Build) for each example problem. To run the harmonic oscillator, for example, input the following.

```
cd Source
./harmonic_oscil.exe
```

Each problem has its own input file in the Input directory off of Iterative_Volterra_Propagator where you can change things like the number of quadrature points used, the maximum number of iterations allowed, and even what method is used to treat matrix exponentials (if applicable). To edit the inputs for the harmonic oscillator without leaving the Source directory, input the following (with your favorite editor in place of emacs).

```
emacs ../../Input/harmonic_oscillator_input.in
```
While not necessary, you can also specify the input file in the command line when executing one of the problems. For example:

```
./model_ode.exe ../../Input/model_ode.input.in
```

Numerical results for each problem print out in the terminal when run. Once finished, the Build directory can be safely deleted without affecting the rest of the code.

## Sample Input/Output
As mentioned above, input files for each example problem are provided in the Input directory. If you run the two level atom with the default inputs by entering the command
```
./two_level.exe
```
the following outputs are displayed (with error and run times varying slightly based on hardware).

```
************************************
Beginning computations!
************************************
Problem: Driven Two-Level Atom
Iteration type: gmres
Pulse amplitude: 0.34906500000000001
Total propagation time:   9000.0000000000000
Propagation step size:   1000.0000000000000
Quadrature type: gauss
Number of quadrature points:           36
************************************
Results:
Maximum ground state error:   8.5321516518632734E-010
Maximum excited state error:   6.0210891827949808E-010
Maximum number of iterations:          43
**************************************************
System Time Elapsed:   6.1000000685453415E-002
Cpu Time Elapsed:   6.0284994542598724E-002
```

## Source Files
Here we give a brief summary of each module in the Source directory.
1. **banded_matrices.90** - Defines a banded symmetric matrix type and contains routines for multiplying with them. Throughout, we assume that the Hamiltonian of the TDSE being solved can be represented as a banded symmetric matrix. Accordingly, this module also contains a routine for constructing a Hamiltonian on a grid via finite difference methods.
2. **fconfig.f90** - Sets up key-value read from files. [Originally due to Kyle Gerheiser.](https://github.com/kgerheiser/fconfig)
3. **grid.f90** - Constructs a spatial grid.
4. **harmonic_oscillator.f90** - Sets up and solves the driven harmonic oscillator TDSE.
5. **integral_method.f90** - Contains main subroutines for ITVOLT: iterative_loop runs the Jacobi or Gauss/Seidel versions while linear_solve runs the GMRES version. Each is to be called once on $[\tau_j, \tau_{j+1}]$.
6. **model_ode.f90** - Sets up and solves the one-dimensional ODE mentioned above.
7. **parameter_read.f90** - Reads problem and method parameters from input files.
8. **parameters.f90** - Defines all common parameters.
9. **potential.f90** - Defines various spatial potentials.
10. **propagator.f90** - Contains code for various methods of computing matrix exponentials. Each method consists of an initialization routine (to be called once on each interval) and an actual propagator routine, which applies a matrix exponential to a vector.
11. **pulse_module.f90** - Defines various time-dependent pulses. More generally, the pulse is the time-dependent function that multiplies the fixed matrix ${\bf V}$ to obtain the time-dependent part of ${\bf H}$.
12. **two_channel.f90** - Sets up and solves the two-channel problem of Wang and Wang.
13. **two_level_atom.f90** - Sets up and solves the TDSE for a two-level atom exposed to a laser. 

Routines for computing Gauss-Lobatto quadrature points and the corresponding set of Lagrange weights as well as a GMRES routine due to Frayssé et al. [2] can be found in the Library directory.

## Contact
Thanks for taking the time to look through our work! If you have questions or feedback about using this repository, reach out to Ryan Schneider at ryschnei@ucsd.edu.

## References
<a id="1">[1]</a> 
W. Wang and X. Wang. *A generalized block-by-block method for the system of linear Volterra equations of the second kind.* Wuhan University Journal of Natural Sciences (2011). 


<a id="2">[2]</a> 
V. Frayssé, L. Giraud, S. Gratton, and J. Langou. *Algorithm 842: A set of GMRES routines for real and complex arithmetics on high performance computers.* ACM Transactions on Mathematical Software  (2005).


