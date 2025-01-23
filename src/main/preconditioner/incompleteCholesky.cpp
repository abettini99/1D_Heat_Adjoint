#include "CoreIncludes.hpp"
#include "preconditioners.hpp"

Eigen::SparseMatrix<f64> Preconditioner::incompleteCholesky(Eigen::SparseMatrix<f64> A){

    // Get size of matrix. 
    u32 n = A.rows(); /**< the #rows and #cols of the preconditioner M(n,n), should be equal to A.cols(). */

    // Get Jacobi Preconditioner
    Eigen::SparseMatrix<f64> L(n,n); /**< Preconditioner L*/
    // for (u32 i=0; i<n; i++){


    // }

    return L;

}