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
    CHECK_FATAL_ASSERT(X.rows() > 255, "Number of interpolating values too high")

    EigenDefs::Matrix<f64> A = EigenDefs::Matrix<f64>::Ones(X.rows(), X.rows());

    for (u8 i=0; i<X.rows(); i++){
        A.col(i) = X.pow(X.rows()-(i+1));
    }

}

PolyInterp1D::PolyInterp1D(EigenDefs::Array1D<f64> coeffs_) : coeffs(coeffs_) {}

} // end Polynomials
