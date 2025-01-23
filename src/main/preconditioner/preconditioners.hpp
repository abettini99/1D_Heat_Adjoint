#pragma once

#include "CoreIncludes.hpp"

namespace Preconditioner {

    // Diagonal Preconditioner
    // - see https://diamhomes.ewi.tudelft.nl/~mvangijzen/PhDCourse_DTU/LES5/TRANSPARANTEN/les5.pdf 
    // - see Section 4.1 https://homepage.tudelft.nl/d2b4e/burgers/lin_notes.pdf
    Eigen::SparseMatrix<f64> Jacobi(Eigen::SparseMatrix<f64> A);


    Eigen::SparseMatrix<f64> incompleteCholesky(Eigen::SparseMatrix<f64> A);
    
}