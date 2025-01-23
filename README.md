
# **Finite Difference Poisson — C++/Eigen + Python**

This repository provides an instructive example of a somewhat-simplistic solver for the 2D Poisson problem. We use C++ for the solver itself, and Python 3 for the postprocessing step. This is the first part of an extensive series to the generalisation to a SEEC (Spectral Element Exterior Calculus) solver. We begin by constructing a simple framework for simple problems, and extending it to suit increasingly-complex problems. The first few repositories are geared towards exposing interested people to lower-level languages. The latter repositories will be about exposing interested people to (more-)modern methods and techniques to solve problems.

### Motivation
Currently in the TUDelft Aerospace Faculty, most numerics is taught at the bachelor-level and masters-level via Python/Matlab. This is perfectly fine as most examples are simple enough to be done in these languages. However, very little exposure to lower-level languages and computing codes leaves most numericist-to-be at a disadvantage in the future. We wish to remedy this. 

### Goal
By the end of this repository, you should know the following:
* Why we are setting up a framework for such a simple problem.
* Why/how we setup typedefs and logging methods (optional).
* How we setup the problem to be solved.
* How we solved the problem using iterative methods (optional).
* Why we save the solution the way it is, and why we use Python to postprocess it.

## Requirements

### LaTeX

### Cmake

### COMPILERS ADD INFO HERE

### Python

## Framework