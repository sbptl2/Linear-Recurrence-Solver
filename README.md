# Linear-Recurrence-Solver
A simple tool powered by Sage to solve non-homogeneous linear recurrences

This program reduces non-homogeneous linear recurrences into homogeneous ones by using the following trick:
  The recurence a_(n+k)*x_k + a_(n+k-1)*x_(n+k-1) + ... a_0*x_n = p_1(n)*e_1(n) + ... + p_z(n)*e_z(n)
    where a_i are complex coeifficients, x_i is the sequence, p_i are polynomial functions, and e_i are exponential functions
  can be reduced to a homogeneous system as the righthand side satisfies a linear homogeneous recurrence equation by means of 
  adding, scaling, and shifting the original equation. This trick that I came up with, makes these equations fairly easy to solve
  deterministically and helps to avoid the normal guessing that is required to solve these problems.
