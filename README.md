# Lightweight Linear Algebra Library FORTRAN

A relatively simple, lightweight module of FORTRAN subroutines that perform some common linear algebra procedures. Particularly useful for coding user-defined material models (e.g., plasticity) within the finite element method. 

Sometimes you just want a simple subroutine that will multiply two matrices together, or solve a linear system of equations. While optimized, established libraries do exist freely online (e.g., LINPACK), the user may not need something that advanced for a prototype code. This repository contains some common linear alebra procedures in FORTRAN that I have required in the past for developing new user-material models within the finite element method, such as the UMAT subroutine in Abaqus Standard/Explicit. 

The include FORTRAN subroutines have been verified through a series of simple use cases. Additionally, they are thread safe. The code has been written in FORTRAN 95, and is suitable for modern FORTRAN compilers, including both the GNU Fortran compiler ([link](https://gcc.gnu.org/fortran/)) and the Intel Fortran compiler ([link](https://software.intel.com/en-us/fortran-compilers)).

In the included FORTRAN module, you can find codes that do the following:

*  Matrix transpose
*  Inner dot product between matrices
*  A horizontal vector dotted with a matrix
*  A matrix dotted with a vertical vector
*  Calculating the Frobenius norm of a matrix
*  Calculating the vector norm (i.e., magnitude)
*  Getting the outer product between two vectors, which results in a matrix
*  Inner dot product between vectors
*  Getting the three invariants of a 3 X 3 tensor
*  Calculating the principal values of 3 X 3 symmetric tensor
*  Getting the eigenvalues and eigenvectors of an n X n symmetric tensor. This one was actually quite tricky, and is based on using the Householder reduction of a matrix to make it tridiagonal, and then using the QL algorithm with implicit shifts to determine the eigenvalues.
*  Transforming a matrix into its LU-decomposition (often needed for other routines, like calculating the determinant, solving a linear system of equations, and calculating the inverse)
*  Calculating the determinant for an n X n matrix
*  Solving for the inverse of an n X n matrix
*  Solving a linear system of equations
*  And few miscellaneous routines that transform Voigt notation into Mandel notation, which I find to be more consistent and easier to avoid programming errors with.

All of the subroutines can be found in "nmoser_linear_algebra.f". All of the subroutines do have a comment block. For some examples on how to properly call these functions, refer to main file that drives the verification of these subroutines, "main_verification.f". This one is a bit more messy and less commented (my apologies), but it essentially tests each subroutine and writes a summary of the tests in an output text file. An example Makefile is also given to compile the two FORTRAN files using GNU gfortran.  