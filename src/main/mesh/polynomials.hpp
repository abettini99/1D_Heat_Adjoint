#pragma once

#include "CoreIncludes.hpp"

/************************************************************************************************************************ 
 *  @brief Any polynomial-related functions/classes are represented in this namespace.
 * 
 *  @details
 *  This namespace serves to identify any-and-all operations related to polynomials.
 ************************************************************************************************************************/
namespace Polynomials{

/************************************************************************************************************************ 
 *  @brief Returns the n-th Legendre polynomial between (-1,1) at specified location.
 * 
 *  @details
 *  Uses Bonnet's formula to recursively construct the Legendre polynomial. Adopted from <a href="https://colab.research.google.com/github/caiociardelli/sphglltools/blob/main/doc/L3_Gauss_Lobatto_Legendre_quadrature.ipynb#scrollTo=Yi60qASPO7tg">here</a>.
 * 
 *  @param n       legendre polynomial order.
 *  @param xi      position to evaluate the legendre polynomial.
 * 
 *  @return float
 ************************************************************************************************************************/ 
f64    Legendre(u8 n, f64 xi);

/************************************************************************************************************************ 
 *  @brief Returns the first derivative of the n-th Legendre polynomial between (-1,1) at specified location.
 * 
 *  @details
 *  Adopted from <a href="https://colab.research.google.com/github/caiociardelli/sphglltools/blob/main/doc/L3_Gauss_Lobatto_Legendre_quadrature.ipynb#scrollTo=Yi60qASPO7tg">here</a>.
 * 
 *  @param n       legendre polynomial order.
 *  @param xi      position to evaluate the first derivative of the legendre polynomial.
 * 
 *  @return float
 ************************************************************************************************************************/ 
 f64  d1Legendre(u8 n, f64 xi);

/************************************************************************************************************************ 
 *  @brief Returns the second derivative of the n-th Legendre polynomial between (-1,1) at specified location.
 * 
 *  @details
 *  Adopted from <a href="https://colab.research.google.com/github/caiociardelli/sphglltools/blob/main/doc/L3_Gauss_Lobatto_Legendre_quadrature.ipynb#scrollTo=Yi60qASPO7tg">here</a>.
 * 
 *  @param n       legendre polynomial order.
 *  @param xi      position to evaluate the second derivative of the legendre polynomial.
 * 
 *  @return float
 ************************************************************************************************************************/ 
f64  d2Legendre(u8 n, f64 xi);

/************************************************************************************************************************ 
 *  @brief Returns the third derivative of the n-th Legendre polynomial between (-1,1) at specified location.
 * 
 *  @details
 *  Adopted from <a href="https://colab.research.google.com/github/caiociardelli/sphglltools/blob/main/doc/L3_Gauss_Lobatto_Legendre_quadrature.ipynb#scrollTo=Yi60qASPO7tg">here</a>.
 * 
 *  @param n       legendre polynomial order.
 *  @param xi      position to evaluate the third derivative of the legendre polynomial.
 * 
 *  @return float
 ************************************************************************************************************************/ 
f64  d3Legendre(u8 n, f64 xi);

/************************************************************************************************************************ 
 *  @brief An interpolating polynomial that go through a set of points.
 * 
 *  @details
 *  Given an array of X and Y, returns an interpolating polynomial in the form of an array of its coefficients. Calling
 *  the object with some X or an array of X will return the interpolated values.
 ************************************************************************************************************************/ 
class PolyInterp1D{

    public:
        // ---------------- //
        // member functions //
        // ---------------- // 

	    /**< Default construction that takes in an X and Y array and fits a polynomial through it */
        PolyInterp1D(EigenDefs::Array1D<f64> X, EigenDefs::Array1D<f64> Y);
	    
        /**< Default construction that takes in an array of coefficients C */
        PolyInterp1D(EigenDefs::Array1D<f64> coeffs_);

	    /**< Disabled construction using another PolyInterp1D */
        PolyInterp1D(const PolyInterp1D&) = delete;

        /**< Disabled construction by equating to another PolyInterp1D */
        PolyInterp1D& operator =(const PolyInterp1D&) = delete;

        /**< Overloading call operator -> Vector of X position should return interpolated polynomial values */
        EigenDefs::Array1D<f64> operator()(EigenDefs::Array1D<f64> X);

        /**< Overloading call operator -> X position should return interpolated polynomial value */
        f64 operator()(f64 X);

        /************************************************************************************************************************ 
		 *  @brief Returns the derivative of the polynomial as another polynomial.
		 * 
		 *  @details
		 *  Basically just multiplies each of the coefficients by the respective exponent of the attached independent variable.
		 * 
		 *  @return PolyInterp1D using the derivative's coefficients.
		 ************************************************************************************************************************/ 
        PolyInterp1D derivative();

        // ---------------- //
        // member variables //
        // ---------------- // 
        EigenDefs::Array1D<f64> coeffs; /**< Coefficients of interpolating polynomial.*/

    private:

};


} // end Polynomials
