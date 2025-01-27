#include "CoreIncludes.hpp"
#include "polynomials.hpp"
#include <math.h>

namespace Polynomials{

// TODO: Add debug + trace messages
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

// TODO: Add debug + trace messages
f64 d1Legendre(u8 n, f64 xi) {
    f64 tmp = n*(Polynomials::Legendre(n-1, xi) - xi*Polynomials::Legendre(n, xi)) \
                                  / (1 - xi*xi);
    return tmp;
}

// TODO: Add debug + trace messages
f64 d2Legendre(u8 n, f64 xi) {
    f64 tmp = (2*xi*Polynomials::d1Legendre(n, xi) - n*(n+1) \
               * Polynomials::Legendre(n, xi)) / (1 - xi*xi);
    return tmp;
}

// TODO: Add debug + trace messages
f64 d3Legendre(u8 n, f64 xi) {
    f64 tmp = (4*xi*Polynomials::d2Legendre(n, xi) - (n*(n+1)-2) \
               * Polynomials::d1Legendre(n, xi)) / (1 - xi*xi);
    return tmp;
}

PolyInterp1D::PolyInterp1D(EigenDefs::Array1D<f64> X, EigenDefs::Array1D<f64> Y) {

    CHECK_FATAL_ASSERT(X.rows() == Y.rows(), "inputs should have matching dimensions.")
    CHECK_FATAL_ASSERT(X.rows() < 256, "Number of interpolating values too high")
    if (X.rows() > 7) WARN_MSG("PolyInterp1D(X,Y) : Unknown whether Vandermonde matrix will have issues due to repeated exponentiation with X.size() = %i elements.", X.rows())

    EigenDefs::Matrix<f64> A = EigenDefs::Matrix<f64>::Zero(X.rows(), X.rows());
    Eigen::ColPivHouseholderQR<EigenDefs::Matrix<f64>> solver;
    TRACE_MSG("PolyInterp1D(X,Y) : Passed variable declarations / initialisation")
    
    // Vandermonde matrix
    for (u8 i=0; i<X.rows(); i++) {
        A.col(i) = X.pow(i);
    }
    TRACE_MSG("PolyInterp1D(X,Y) : Passed matrix construction") 

    solver.compute(A);
    coeffs = solver.solve(Y.matrix()); // converts array to matrix (basically a vector)
    TRACE_MSG("PolyInterp1D(X,Y) : Passed solver step") 
    
    // DEBUG SUMMARY
    #if RELEASE==0
        std::stringstream printArr;
        printArr << std::fixed << std::setprecision( 4 );
        // print A
	printArr << "PolyInterp1D(X,Y) : A = \n[";
        for (u64 i=0; i<A.rows(); i++) {
	    if (i != 0) printArr << " "; // Extra padding so that all [ line up
            printArr << " [ ";
	    for (u64 j=0; j<A.cols(); j++) {
		printArr << A(i,j) << " ";
            }
	    printArr << "]";
	    if (i != A.rows()-1) printArr << "\n"; // Nicer formatting so that the next ] does not appear on next line
	}
	printArr << " ]";
        DEBUG_MSG("%s", printArr.str().c_str())
        
	printArr.str(std::string());
	// print coeffs 
        printArr << std::fixed << std::setprecision( 4 );
        printArr << "PolyInterp1D(X,Y) : coeffs = [ ";
        for (u64 i=0; i<coeffs.size(); i++) printArr << coeffs[i] << " ";
        printArr << "]";
        DEBUG_MSG("%s", printArr.str().c_str())
    #endif
}

PolyInterp1D::PolyInterp1D(EigenDefs::Array1D<f64> coeffs_) : coeffs(coeffs_) {
    #if RELEASE == 0
    std::stringstream printArr;
    printArr << std::fixed << std::setprecision( 4 );
    printArr << "PolyInterp1D(coeffs) : coeffs = [ ";
    for (u64 i=0; i<coeffs.size(); i++) printArr << coeffs[i] << " ";
    printArr << "]";
    DEBUG_MSG("%s", printArr.str().c_str())
    #endif
}

// TODO: Add assertions, debug + trace messages
EigenDefs::Array1D<f64> PolyInterp1D::operator()(EigenDefs::Array1D<f64> X) {
    
    EigenDefs::Array1D<f64> out = EigenDefs::Array1D<f64>::Zero(X.rows());
    
    for (u8 i=0; i<coeffs.rows(); i++){
        out += coeffs[i]*X.pow(i);
    }

    return out;
}

// TODO: Add assertions, debug + trace messages
f64 PolyInterp1D::operator()(f64 X) {
    
    f64 out = 0;
    
    for (u8 i=0; i<coeffs.rows(); i++){
        out += coeffs[i]*std::pow(X,i);
    }
    
    return out;
}

PolyInterp1D PolyInterp1D::derivative() {

    if (coeffs.rows() == 1) {
	TRACE_MSG("PolyInterp1D.derivative : return 0")
        return PolyInterp1D(EigenDefs::Array1D<f64>::Zero(1));
    }
    else {
        // The constant drops out, hence the -1
        EigenDefs::Array1D<f64> derivCoeffs(coeffs.rows()-1);
        TRACE_MSG("PolyInterp1D.derivative : passed variable declaration")

        // Basic differentiation of a*x^n = n*a*x^{n-1}. the coefficient is n*a
        for (u8 i=1; i<coeffs.rows(); i++){
            derivCoeffs[i-1] = coeffs[i]*i;
        };
        TRACE_MSG("PolyInterp1D.derivative : passed coefficient determination") 
        return PolyInterp1D(derivCoeffs);
    }
}

} // end Polynomials
