#pragma once

#include "CoreIncludes.hpp"

/************************************************************************************************************************ 
 *  @brief Any mesh-related functions/classes are represented in this namespace.
 * 
 *  @details
 *  This namespace serves to identify any-and-all operations related to the mesh. In the FEM-sense, this can include
 *  element spaces, elements, assembly, etc.
 ************************************************************************************************************************/
namespace Mesh{

// struct Simplex{
//     u8 d; /**< the d in d-simplex: 2-simplex is a triangle, 3-simplex is a tetrahedron, etc. */
//     Simplex(u8 d_) : d(d_) {};
// };
// 
// struct Cube{
//     u8          d;    /**< the d in d-cube: 2-cube is a parallelogram, 3-cube is a parallelepiped, etc. */
//     const char* type; /**< Cube type: SEM, Standard */
//     Cube(u8 d_, const char* type_) : d(d_), type(type_) {};
// };

class Geometry{

    public:

	/** Reads in geometry file, setup the mesh automatically here. TODO replace with mesh file if necessary */
        Geometry(f64 boundaries[2], u32 nElements);
        
	/**< Disabled construction using another Geometry */
        Geometry(const Geometry&) = delete;

        /**< Disabled construction by equating to another Geometry */
        Geometry& operator =(const Geometry&) = delete;

    private:

};

class MasterElement{

    public:
	
	/** Master element setup, takes in the number of variables per element, and the general shape. */
        MasterElement(u8 nVars_, u8 nDims_); 

	/** Sets the GLL-Lagrange polynomial order */
	void setGLLOrder(u8 Var, u8 polyOrder);
	
	/** Sets the Gauss-Lagrange polynomial order */
	void setGaussOrder(u8 Var, u8 polyOrder);

    private:
	u8 nVars, nDims;
	std::vector<EigenDefs::Array1D<f64>> weights;
	std::vector<EigenDefs::Array1D<f64>> nodes;
	

};

} // namespace Mesh
