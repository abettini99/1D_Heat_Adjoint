#include "CoreIncludes.hpp"
#include "mesh.hpp"
#include "polynomials.hpp"
#include <math.h>
#include <stdarg.h>

namespace Mesh{

MasterElement::MasterElement(u8 nDims_) : nDims(nDims_) {

    CHECK_FATAL_ASSERT(nDims > 0, "Number of dimensions must be bigger than 0")
}

void MasterElement::setnVars(u8 nVars_) {

    CHECK_FATAL_ASSERT(nDims  > 0, "nDims must be set first before calling upon this function")
    CHECK_FATAL_ASSERT(nVars_ > 0, "Number of variables must be bigger than 0")
    nVars = nVars_;
}

void MasterElement::setLGLOrder(u8 Var, ...){

    DEBUG_MSG("MasterElement.setLGLOrder : ===========")
    DEBUG_MSG("MasterElement.setLGLOrder : Var %i", Var)
    DEBUG_MSG("MasterElement.setLGLOrder : ===========")

    CHECK_FATAL_ASSERT(nDims > 0, "nDims must be set first before calling upon this function")
    CHECK_FATAL_ASSERT(nVars > 0, "nVars must be set first before calling upon this function")
    CHECK_FATAL_ASSERT(nVars > Var, "Variable number accessed too large")
    f64 epsilon = 1e-15;
    std::vector<u8> tmp; 
    i32 polyOrder;
    va_list  argPtr;
    va_start(argPtr, Var);
    for (u8 i=0; i<nDims; i++) {
        polyOrder = va_arg(argPtr, i32);
	if (!polyOrder) break;
	
	INFO_MSG("%i", polyOrder)
	CHECK_FATAL_ASSERT(polyOrder > 1,   "Polynomial order must be larger than 1")
	tmp.push_back( polyOrder);	
	
    }
    va_end(argPtr);
    INFO_MSG("%i %i", tmp.size(), nDims)
    CHECK_FATAL_ASSERT(tmp.size() == nDims, "Number of polyOrder inputs does not match nDim")

    // CHECK_FATAL_ASSERT(polyOrder > 1, "Polynomial Order must be larger than 1")
    // // --------------------- //
    // // LGL Weights and Nodes //
    // // --------------------- // 

    // u8  n       = polyOrder + 1; // Number of gridpoints 
    // EigenDefs::Array1D<f64> x  = EigenDefs::Array1D<f64>::Zero(n);
    // EigenDefs::Array1D<f64> w  = EigenDefs::Array1D<f64>::Zero(n);
    // x[0] = -1;              x[n-1] = 1;
    // w[0] = 2.0/(n*(n-1));   w[n-1] = w[0];

    // u8 n_2 = n/2; // Floor division if n is odd
    // f64 error;
    // f64 xi, dxi;
    // f64 y1, y2, y3;
    // TRACE_MSG("MasterElement.setLGLOrder : Passed variable declaration / initialisation") 

    // for (u8 i=1; i<n_2; i++) {
    //     xi = (1 - 3*(n-2) / (8 * (n-1)*(n-1)*(n-1))) \
    //          * std::cos((4*i+1)*EIGEN_PI/(4*(n-1)+1));

    //     error = 1.;

    //     do{
    //         y1 = Polynomials::d1Legendre(n-1, xi);
    //         y2 = Polynomials::d2Legendre(n-1, xi);
    //         y3 = Polynomials::d3Legendre(n-1, xi);

    //         dxi = 2*y1*y2 / (2*y2*y2-y1*y3);
    //         xi -= dxi;
    //         error = std::abs(dxi);

    //     } while (error > epsilon);

    //     x[i]     = -xi;
    //     x[n-i-1] =  xi;

    //     w[i]     = 2/(n*(n-1)*std::pow(Polynomials::Legendre(n-1,x[i]),2));
    //     w[n-i-1] = w[i];
    // }

    // if (n%2 != 0) {
    //     x[n_2] = 0;
    //     w[n_2] = 2/(n*(n-1)*std::pow(Polynomials::Legendre(n-1,x[n_2]),2));
    // } 

    // nodes.push_back(x);
    // weights.push_back(w);

    // // DEBUG SUMMARY
    // TRACE_MSG("MasterElement.setLGLOrder : Var %i - x, w array pushed", Var)
    // #if RELEASE==0
    //     std::stringstream printArr;
    //     f64 SUM = 0.;
    //     // print x array
    //     printArr << std::fixed << std::setprecision( 4 );
    //     printArr << "MasterElement.setLGLOrder : x = [ ";
    //     for (u8 i=0; i<x.size(); i++) printArr << x[i] << " ";   
    //     printArr << "]";
    //     DEBUG_MSG("%s", printArr.str().c_str())
    //     
    //     printArr.str(std::string());
    //     // print w array + sum
    //     printArr << std::fixed << std::setprecision( 4 );
    //     printArr << "MasterElement.setLGLOrder : w = [ ";
    //     for (u8 i=0; i<w.size(); i++) {printArr << w[i] << " "; SUM += w[i];}  
    //     printArr << "], SUM = " << SUM;
    //     DEBUG_MSG("%s", printArr.str().c_str())
    // #endif

    // // ----------------------------- //
    // // LGL-Lagranges and Derivatives //
    // // ----------------------------- // 

    // EigenDefs::Array1D<f64> y = EigenDefs::Array1D<f64>::Zero(x.rows());
    // std::vector<Polynomials::PolyInterp1D> lagrange_;
    // std::vector<Polynomials::PolyInterp1D> d1lagrange_;
    // for (u8 i=0; i<x.rows(); i++){
    //     DEBUG_MSG("MasterElement.setLGLOrder : -----------")
    //     DEBUG_MSG("MasterElement.setLGLOrder : n = %i", i)
    //     DEBUG_MSG("MasterElement.setLGLOrder : -----------")
    //     y.setZero();
    //     y[i] = 1.;

    //     Polynomials::PolyInterp1D   lagrange__(x,y);
    //     Polynomials::PolyInterp1D d1lagrange__ = lagrange__.derivative();

    //     // Push to subvector
    //     lagrange_.push_back(lagrange__);
    //     d1lagrange_.push_back(d1lagrange__);

    //     // DEBUG SUMMARY
    //     TRACE_MSG("MasterElement.setLGLOrder : P_%i, dP_%i pushed to subvector, n=%i", polyOrder, polyOrder, i)
    //     #if RELEASE==0
    //         std::stringstream printArr;
    //         EigenDefs::Array1D<f64> out; 
    //         // print p(x)
    //         out = lagrange_[i](x);
    //         printArr << std::fixed << std::setprecision( 4 );
    //         printArr << "MasterElement.setLGLOrder : P(x)   = [ ";
    //         for (u64 j=0; j<out.size(); j++) printArr << out[j] << " ";
    //         printArr << "] ?= [ ";
    //         for (u64 j=0; j<out.size(); j++) printArr << y[j] << " ";
    //         printArr << "]";
    //         DEBUG_MSG("%s", printArr.str().c_str())

    //         printArr.str(std::string());
    //         // print d1p(x)
    //         out = d1lagrange_[i](x);
    //         printArr << std::fixed << std::setprecision( 4 );
    //         printArr << "MasterElement.setLGLOrder : d1P(x) = [ ";
    //         for (u64 j=0; j<out.size(); j++) printArr << out[j] << " ";
    //         printArr << "]";
    //         DEBUG_MSG("%s", printArr.str().c_str())
    //     #endif
    // }
    // // Push to vector
    // lagrange.push_back(lagrange_);
    // d1lagrange.push_back(d1lagrange_);
    // TRACE_MSG("MasterElement.setLGLOrder : { P_%i, dP_%i } subvector pushed to main std::vector", polyOrder, polyOrder)

    // INFO_MSG("Variable %i - FEM space set to piecewise LGL-Lagrange polynomials of order %i", Var, polyOrder)

}

Geometry::Geometry(EigenDefs::Array1D<f64> x1) : nDims(1) {

    MasterElement = Mesh::MasterElement(nDims);
    
    INFO_MSG("%iD cartesian grid established", nDims)

}

Geometry::Geometry(EigenDefs::Array1D<f64> x1, EigenDefs::Array1D<f64> x2) : nDims(2) {

    MasterElement = Mesh::MasterElement(nDims);
    
    INFO_MSG("%iD cartesian grid established", nDims)
}

Geometry::Geometry(EigenDefs::Array1D<f64> x1, EigenDefs::Array1D<f64> x2, EigenDefs::Array1D<f64> x3) : nDims(3) {

    MasterElement = Mesh::MasterElement(nDims);
    
    INFO_MSG("%iD cartesian grid established", nDims)
}

} // end Mesh
