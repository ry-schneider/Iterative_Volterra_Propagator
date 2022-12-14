# ITVOLT: An Iterative Solver for Volterra Integral Equations
Volterra integral equations of the second kind take the form

$$ {\bf f}(t) = {\bf g}(t) + \int_{t_0}^t {\bf K}(t,t') {\bf f}(t') dt' $$

on an interval $t_0 \le t \leq t_f$ for two vector-valued functions ${\bf f}$ and ${\bf g}$ and a matrix integral kernel ${\bf K}$. Given an inhomogeneity ${\bf g}$ and a kernel ${\bf K}$, the equations are solved for the unknown function ${\bf f}$. 
 
ITVOLT, short for Iterative Volterra Propagator, is a novel method for solving these equations via global Lagrange interpolation. Given an approximation ${\bf f}^{(k)}$ of ${\bf f}$, the method proceeds by choosing a set of quadrature points in $[t_0, t_f]$ and expanding ${\bf K}(t,t'){\bf f}^{(k)}(t')$ in Lagrange polynomials, finding ${\bf f}^{(k+1)}$ by evaluating the Volterra equation via quadrature. The numerical details of the method are presented in a [forthcoming paper.](https://arxiv.org/abs/2210.15677)

This repository contains Fortran 90 code that applies ITVOLT to four example problems:
 * A simple two channel problem of Wang and Wang [2].
 * The following one dimensional ODE (with initial condition $\psi(0) = 1$):
 
 $$ \left[ i \frac{ \partial}{\partial t} - t \right] \psi(t) = 0. $$
 
 * The time dependent Schrödinger equation (TDSE) for a two level atom exposed to a laser.
 * The TDSE for a driven harmonic oscillator.
 
The latter two problems, as well as the conversion of the TDSE into a Volterra integral equation, are covered in work of Ndong et al. [1].
 
Of particular interest for other researchers, the harmonic oscillator example demonstrates how the repository can be used to solve the TDSE for any problem that satisfies the following criteria:
1. The Hamiltonian can be represented by a symmetric banded matrix.
2. The time dependent piece of the Hamiltonian can be written as $E(t)*{\bf V}$ for some function of time $E$ and some fixed, symmetric banded matrix ${\bf V}$.

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

## Contact
Thanks for taking the time to look through our work! If you have questions or feedback about using this repository, reach out to Ryan Schneider at ryschnei@ucsd.edu.

## References
<a id="1">[1]</a> 
M. Ndong, H. Tal-Ezer, R. Kosloff, and C.P. Koch. *A Chebyshev propagator with iterative time ordering for explicitly time-dependent hamiltonians.* The Journal of Chemical Physics (2010).

<a id="2">[2]</a> 
W. Wang and X. Wang. *A generalized block-by-block method for the system of linear Volterra equations of the second kind.* Wuhan University Journal of Natural Sciences (2011). 
