#include "CoreIncludes.hpp"
#include "CG.hpp"

namespace KrylovSolver{

CG::CG(Eigen::SparseMatrix<f64> &A) : A(A) {

    // Set new vectors
    u32 n = A.rows();
    u32 m = A.cols();
    CHECK_FATAL_ASSERT(n==m, "Number of rows and columns of sparse matrix A do not match.")

    rk.setZero(m);
    rkp1.setZero(m);
    zk.setZero(m);
    zkp1.setZero(m);
    pk.setZero(m);
    qk.setZero(m);  
}

void CG::solve(EigenDefs::Vector<f64> &u,
               EigenDefs::Vector<f64> &b,
               f64 tol, u32 iterMax){

    // Initialization
    u32 iter = 0;     /**< Iterate count */
    f64 err = 1./0.;  /**< residual error */
    rk = b - A*u;     // Initial guess
    // zk = Mm1*rk;
    zk = rk;
    pk = zk;

    // N.B. We write it this way to skip the if-else statement in Figure 5.2 of Henk van der Vorst 2003
    do {
        // Update iterate
        qk     = A*pk;
        alphak = rk.dot(zk) / pk.dot(qk);
        u      = u + alphak*pk;

        // Update residual
        rkp1   = rk - alphak*qk;
        err    = std::sqrt( rk.dot(rk)/rk.size() );

        // Termination criteria
        CHECK_FATAL_ITERERROR(iter, err);
        INFO_MSG("iter = %-5u err = %1.4e", iter, err); 
        if (tol > err) break;

        // Calculate preconditioning residual vector
        // zkp1 = Mm1*rkp1;                
        zkp1   = rkp1;         

        // Update search direction
        betak  = rkp1.dot(zkp1) / rk.dot(zk);
        pk     = zkp1 + betak*pk;    

        // Update iteration
        iter++;
        rk     = rkp1; 
        zk     = zkp1;

    } while (iter < iterMax); 
}

} // end KrylovSolver




