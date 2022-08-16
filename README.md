# Iterative_Volterra_Propagator
Volterra integral equations of the second kind take the form
f (t) = g(t) + ∫ K(t, s)f(s)ds, t_0 ≤ t ≤ t_f (1.1)
for two functions f and g and an integral kernel K. Studied first by Liouville [8] and later Volterra [14],
these equations arise in a variety of disciplines, including material science (linear viscoelasticity equations)
and probability (the renewal equation). More generally, any first order differential equation can be converted
to a Volterra integral equation by integrating. Given an inhomogeneity, g and and integral kernel, K, the
equations are solved for the unknown, f .
Numerical approaches to these problems have been developed for several decades and are the subject
of numerous papers. Early approaches focused on applying a quadrature to the integral in the equation,
beginning with the simple trapezoidal rule and generalizing to more complicated quadrature rules. For an in-
depth survey of these techniques as well as a discussion of analytic solutions to (1.1), see work of Brunner and
van der Houwen [1]. Recent approaches make use of more complicated mathematical machinery, including
interpolation [13], block methods [15], and homotopy theory [2].
In this paper, we present a new iterative method for solving Volterra equations based on Lagrange
interpolation. In sections two and three, we discuss the numerical details of the method and demonstrate its
utility on a handful of toy problems. In the remainder of the paper, we consider in depth the application of
our method to the time dependent Schr ̈odinger equation, illustrating its ability to achieve highly accurate
results on a number of physical problems.
