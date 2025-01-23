#pragma once

#include "CoreIncludes.hpp"

/************************************************************************************************************************ 
 *  @brief All Krylov iterative solvers are stored underneath this namespace.
 * 
 *  @details
 *  This currently includes CG, BiCGstab(l).
 ************************************************************************************************************************/ 
namespace KrylovSolver{
/************************************************************************************************************************ 
 *  @brief A hybrid Petrov-Galerkin iterative solver. Used with any type of A matrix.
 * 
 *  @details
 *  The idea behind iterative solvers is that we want better estimates of our exact solution:
 *  x1 <-- x0
 *  x2 <-- x1
 *  x3 <-- x2
 *  and so forth until x*, our exact solution. If we define r_k = x_k+1 - x_k, then we can show that
 *  x* = x0 + r0 + r1 + r2 + ...
 * 
 *  where r0 + r1 + r2 + ... can be shown to be a linear combination of vectors in span(r0, Ar0, A^2r0, ... A^n r0) = 
 *  K^n+1(A, r0). This set is called a Krylov subspace. We normally rewrite the basis of the Krylov subspace than the 
 *  ones spanned by powers of the sparse matrix A because the basis vectors tend to become linearly-dependent on each 
 *  other in finite-precision. You need a suitable set of basis in order to explore/design algorithms. 
 * 
 *  Ritz-Galerkin approaches seek to construct a new residual estimate r_{k+1} that is orthogonal to the rest of the 
 *  subspace spanned by v_1, v_2, ..., v_k, where v_k are the basis of the Krylov subspace K^{k+1}(A; r_0) determined by 
 *  the Arnoldi algorithm (see Chapter 3.3 of Henk van der Vorst 2003, around Fig.~3.1). Since the new r_{k+1} is 
 *  orthogonal to the rest of the subspace, the residual can act as the v_{k+1} basis.
 * 
 *  It turns out that based on a certain relation when using symmetric A matrix, you can design the next basis vector 
 *  of the Krylov subspace by only using all previous basis vectors, with the equations established by noting that we
 *  want a Ritz-Galerkin approach. As it turns out, storing all basis vectors is undesirable given that you can have
 *  thousands and thousands of iterations. CG makes a further assumption on the A matrix to solve this issue, which is 
 *  that it is positive definite. This allows for a very simple recursion relation that allows us to not have to store
 *  all basis vectors anymore.
 * 
 *  As it turns out, the 'search' vectors that lead to orthogonal residual vectors are conjugate (A-orthogonal), hence
 *  the name of the method.
 * 
 *  * see "Iterative Krylov Methods for Large Linear Systems" by Henk van der Vorst 2003
 *  * see "A Brief Introduction to Krylov Space Methods for Solving Linear Systems" by Martin H. Gutknecht 2007
 *  * see Section 3.1 https://homepage.tudelft.nl/d2b4e/burgers/lin_notes.pdf
 *  * see Section 4.1 https://homepage.tudelft.nl/d2b4e/burgers/lin_notes.pdf
 *  * see https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
 ************************************************************************************************************************/ 
template<u32 level> class BiCGstab{

    public:
        // ---------------- //
        // member functions //
        // ---------------- // 

        /**< Default construction takes a reference to the sparse A matrix, and resizes all internal vectors to the appropriate shape */
        BiCGstab(Eigen::SparseMatrix<f64> &A);    

        /**< Disabled construction using another BiCGstab solver */
        BiCGstab(const BiCGstab&) = delete;             

        /**< Disabled construction by equating to another BiCGstab solver */
        BiCGstab& operator =(const BiCGstab&) = delete; 



        /************************************************************************************************************************ 
         *  @brief Runs through the conjugate-gradient algorithm to find the solution to Au = b.
         * 
         *  @details
         *  The details of this methodology can be found in the main class descriptor, and algorithms in the references provided.
         * 
         *  @param u       reference to the solution vector of the system Au = b.
         *  @param b       reference to the forcing vector of the system Au = b.
         *  @param tol     tolerance for convergence, default 1e-15.
         *  @param maxiter maximum number of iterations for convergence, default 5000.
         * 
         *  @return None
         ************************************************************************************************************************/ 
        void solve(EigenDefs::Vector<f64> &u,
                   EigenDefs::Vector<f64> &b,
                   f64 tol = 1e-15,
                   u32 maxiter = 5000);



    private:
        // ---------------- //
        // member functions //
        // ---------------- // 
        
        // ---------------- //
        // member variables //
        // ---------------- // 
        const u8 l = level;

        EigenDefs::Vector<f64> rk;            /**< residual vector */
        std::vector< EigenDefs::Vector<f64> > hu;
        std::vector< EigenDefs::Vector<f64> > hr;

        EigenDefs::Vector<f64> tr0;
        
        f64 alpha, beta, omega; /**< update coefficients */
        f64 rho0, rho1;
        f64 tau[level+1][level+1];  // ignore index 0, hence the +1;
        f64 sigma[level+1];     // ignore index 0;
        f64 gam;
        f64 gamma[level+1];
        f64 gammap[level+1];
        f64 gammapp[level+1];
        
        // Eigen::SparseMatrix<f64> A;      /**< Internal reference of the sparse matrix*/
        // EigenDefs::Vector<f64> rk, rkp1; /**< residual vector */
        // EigenDefs::Vector<f64> zk, zkp1; /**< preconditioned residual vector */
        // EigenDefs::Vector<f64> pk;       /**< search/conjugate direction vector */
        // EigenDefs::Vector<f64> qk;       /**< search/conjugate direction vector */
        // f64 alphak, betak;               
    
};

} // end KrylovSolver



