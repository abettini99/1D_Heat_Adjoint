#pragma once

#include "CoreIncludes.hpp"
#include "polynomials.hpp"

/************************************************************************************************************************ 
 *  @brief Any mesh-related functions/classes are represented in this namespace.
 * 
 *  @details
 *  This namespace serves to identify any-and-all operations related to the mesh. In the FEM-sense, this can include
 *  element spaces, elements, assembly, etc.
 ************************************************************************************************************************/
namespace Mesh{

class Geometry{

    public:

	    /**< Reads in geometry file, setup the mesh automatically here. TODO replace with mesh file if necessary */
        Geometry(f64 boundaries[2], u32 nElements);
        
	    /**< Disabled construction using another Geometry */
        Geometry(const Geometry&) = delete;

        /**< Disabled construction by equating to another Geometry */
        Geometry& operator =(const Geometry&) = delete;

    private:

};

class MasterElement{

    public:

        // ---------------- //
        // member functions //
        // ---------------- // 
	
	    /**< SEM master element setup, takes in the number of variables per element. */
        MasterElement(u8 nVars_, u8 nDims_); 
        
	    /**< Disabled construction using another MasterElement */
        MasterElement(const MasterElement&) = delete;

        /**< Disabled construction by equating to another MasterElement */
        MasterElement& operator =(const MasterElement&) = delete;

	    /************************************************************************************************************************ 
		 *  @brief Sets the corresponding variable's FEM space to LGL-Lagrange polynomials
		 * 
		 *  @details
		 *  Uses Halley's method to calculate LGL nodes and weights.
		 *  Adopted from <a href="https://colab.research.google.com/github/caiociardelli/sphglltools/blob/main/doc/L3_Gauss_Lobatto_Legendre_quadrature.ipynb#scrollTo=Yi60qASPO7tg">here</a>.
		 * 
		 *  @param Var        Variable who's space is to be set.
		 *  @param polyOrder  Lagrange polynomial order.
		 * 
		 *  @return None
		 ************************************************************************************************************************/ 
		void setGLLOrder(u8 Var, u8 polyOrder, f64 epsilon=1e-15);
	    
	    /**< Sets the Gauss-Lagrange polynomial order */
	    void setGaussOrder(u8 Var, u8 polyOrder);

		// TMP //

		EigenDefs::Array1D<f64> TMP1(u8 Var) {
			return weights[Var];
		}

		EigenDefs::Array1D<f64> TMP2(u8 Var) {
			return nodes[Var];
		}

        // ---------------- //
        // member variables //
        // ---------------- // 

    private:

        // ---------------- //
        // member variables //
        // ---------------- // 
	    u8 nVars, nDims;
	    std::vector<EigenDefs::Array1D<f64>> weights;                    /**< Master element weights for reference quadrature  */
	    std::vector<EigenDefs::Array1D<f64>> nodes;                      /**< Master element nodes for reference quadrature */
		std::vector<std::vector<Polynomials::PolyInterp1D>>   lagrange;  /**< Lagrange functions that fit through master element nodes, access is lagrange[Var][nPoly] */
		std::vector<std::vector<Polynomials::PolyInterp1D>> d1lagrange;  /**< Lagrange functions that fit through master element nodes, access is d1lagrange[Var][nPoly] */

};

} // end Mesh
