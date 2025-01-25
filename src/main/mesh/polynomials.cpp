#include "CoreIncludes.hpp"
#include "Polynomials.hpp"
#include <math.h>

namespace Polynomials{

f64   Legendre(u8 n, f64 xi) {
    if (n == 0) {
        f64 tmp = 1.;
        return tmp;
    }
    else if (n == 1) {
        return xi;
    }
    else {
        f64 fP = 1.;
        f64 sP = xi;
        f64 nP = 0.;

        for (u8 i=2; i<n+1; i++){
            nP = ((2*i-1)*xi*sP-(i-1)*fP)/i;
            fP = sP; sP = nP;
        }
        return nP;
    }
}

f64 d1Legendre(u8 n, f64 xi) {
    f64 tmp = n*(Polynomials::Legendre(n-1, xi) - xi*Polynomials::Legendre(n, xi)) \
                                  / (1 - xi*xi);
    return tmp;
}

f64 d2Legendre(u8 n, f64 xi) {
    f64 tmp = (2*xi*Polynomials::d1Legendre(n, xi) - n*(n+1) \
               * Polynomials::Legendre(n, xi)) / (1 - xi*xi);
    return tmp;
}

f64 d3Legendre(u8 n, f64 xi) {
    f64 tmp = (4*xi*Polynomials::d2Legendre(n, xi) - (n*(n+1)-2) \
               * Polynomials::d1Legendre(n, xi)) / (1 - xi*xi);
    return tmp;
}

PolyInterp1D::PolyInterp1D(EigenDefs::Array1D<f64> X, EigenDefs::Array1D<f64> Y) {

    CHECK_FATAL_ASSERT(X.rows() == Y.rows(), "inputs should have matching dimensions.")
    CHECK_FATAL_ASSERT(X.rows() < 256, "Number of interpolating values too high")

    EigenDefs::Matrix<f64> A = EigenDefs::Matrix<f64>::Zero(X.rows(), X.rows());
    Eigen::ColPivHouseholderQR<EigenDefs::Matrix<f64>> solver;
    
    // Vandermonde matrix
    for (u8 i=0; i<X.rows(); i++) {
        A.col(i) = X.pow(i);
    }
    //INFO_MSG("%-4f", A[0,0])
    solver.compute(A);
    coeffs = solver.solve(Y.matrix()); // converts array to matrix (basically a vector)
}

PolyInterp1D::PolyInterp1D(EigenDefs::Array1D<f64> coeffs_) : coeffs(coeffs_) {}

EigenDefs::Array1D<f64> PolyInterp1D::operator()(EigenDefs::Array1D<f64> X) {
    
    EigenDefs::Array1D<f64> out(X.rows());
    
    for (u8 i=0; i<coeffs.rows(); i++){
        out += coeffs[i]*X.pow(i);
    }

    return out;
}

f64 PolyInterp1D::operator()(f64 X) {
    
    f64 out;
    
    for (u8 i=0; i<coeffs.rows(); i++){
        out += coeffs[i]*std::pow(X,i);
    }
    
    return out;
}

PolyInterp1D PolyInterp1D::derivative() {

    if (coeffs.rows() == 1) {
        return PolyInterp1D(EigenDefs::Array1D<f64>::Zero(1));
    }
    else {
        // The constant drops out, hence the -1
        EigenDefs::Array1D<f64> derivCoeffs(coeffs.rows()-1);
        
        // Basic differentiation of a*x^n = n*a*x^{n-1}. the coefficient is n*a
        for (u8 i=1; i<coeffs.rows(); i++){
            derivCoeffs[i-1] = coeffs[i]*i;
        };
        return PolyInterp1D(derivCoeffs);
    }
}

} // end Polynomials
