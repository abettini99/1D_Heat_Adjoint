#include "CoreIncludes.hpp"
#include "preconditioners.hpp"

Eigen::SparseMatrix<f64> Preconditioner::Jacobi(Eigen::SparseMatrix<f64> A){

    // Get size of matrix. 
    u32 n = A.rows(); /**< the #rows and #cols of the preconditioner M(n,n), should be equal to A.cols(). */

    // Get Jacobi Preconditioner
    Eigen::SparseMatrix<f64> M(n,n); /**< Preconditioner M*/
    M = A.diagonal().asDiagonal();

    return M;

}
