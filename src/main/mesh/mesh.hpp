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

class MasterElement{

    public:

        // ---------------- //
        // member functions //
        // ---------------- // 
	
	    /**< Empty constructor */
        MasterElement() : nDims(0), nVars(0) {}; 

	    /**< SEM master element setup, takes in the d-cube argument. The master element is dependent on what elements exist in the geometry. */
        MasterElement(u8 nDims_); 

	    /************************************************************************************************************************ 
		 *  @brief Sets the number of variables that this master element handles.
		 * 
		 *  @details
		 *  Necessary to develop the polynomials defined on the master element. Each variable can use different polynomials.
		 *  We keep this separate from the constructor as the dimensions are set by the number of dimensions available in the
		 *  Geometry, but the number of variables is problem dependent (and not geometry dependent).
		 * 
		 *  @param nVars      Number of variables that this master element handles.
		 * 
		 *  @return None
		 ************************************************************************************************************************/ 
		void setnVars(u8 nVars_);

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
		void setLGLOrder(u8 Var, u8 polyOrder, f64 epsilon=1e-15);
	    
	    /**< Sets the Gauss-Lagrange polynomial order TODO: implement, used for pressure in a ((P_n^u)^3 U (P_{n-2}^p)) space scheme for NS*/
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
	    std::vector<EigenDefs::Array1D<f64>> weights;                    /**< Master element weights for reference quadrature */
	    std::vector<EigenDefs::Array1D<f64>> nodes;                      /**< Master element nodes for reference quadrature */
		std::vector<std::vector<Polynomials::PolyInterp1D>>   lagrange;  /**< Lagrange functions that fit through master element nodes, access is lagrange[Var][nPoly] */
		std::vector<std::vector<Polynomials::PolyInterp1D>> d1lagrange;  /**< Lagrange functions that fit through master element nodes, access is d1lagrange[Var][nPoly] */

};


class Geometry{

    public:

        // ---------------- //
        // member functions //
        // ---------------- // 

	    /**< 1D grid */
        Geometry(EigenDefs::Array1D<f64> x1);
		
	    /**< 2D tensor-grid TODO*/
        Geometry(EigenDefs::Array1D<f64> x1, EigenDefs::Array1D<f64> x2);

	    /**< 3D tensor-grid TODO*/
        Geometry(EigenDefs::Array1D<f64> x1, EigenDefs::Array1D<f64> x2, EigenDefs::Array1D<f64> x3);
        
		/**< TODO */
		// Geometry(FILE)

	    /**< Disabled construction using another Geometry */
        Geometry(const Geometry&) = delete;

        /**< Disabled construction by equating to another Geometry */
        Geometry& operator =(const Geometry&) = delete;

        // ---------------- //
        // member variables //
        // ---------------- // 

		std::vector<EigenDefs::Matrix<f64>> dx_dxi;
		std::vector<u64> nElems;
		u8  nVars;
		u8  nDims;
		Mesh::MasterElement MasterElement;

    private:

};

} // end Mesh
