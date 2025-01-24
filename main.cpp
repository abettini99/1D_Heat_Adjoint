/************************************************************************************************************************
 * In these two header files, we include:
 *    - Typedefs:   used to create an additional name for another data type, but does not create a new type.
 *    - Type alias: same as typedefs, but allows for templates, e.g. Vector<f64>, Vector<i32>.
 ************************************************************************************************************************/
#include "CoreIncludes.hpp"
#include "solver/CG.hpp"
#include "solver/BiCGstab_l_.hpp"
#include "mesh/mesh.hpp"
#include "mesh/polynomials.hpp"
#include "mesh/valueSource.hpp"

#include <vector>
#include <fstream>
#include <iostream>

/************************************************************************************************************************
 * Solve -div(grad(u)) = f, using FDM
 ************************************************************************************************************************/
int main(){

    //## ================== ##//
    //## Provide parameters ##//
    //## ================== ##//
    const u32 imax  = 101;              /**< #gridpoints in x */
    const u32 jmax  = 101;              /**< #gridpoints in y */
    const f64 Lx[2] = {0., 1.*EIGEN_PI}; /**< domain endpoints in x */
    const f64 Ly[2] = {0., 1.*EIGEN_PI}; /**< domain endpoints in y */

    //## ==================== ##//
    //## Calculate parameters ##//
    //## ==================== ##//
    const u32 n = (imax-2)*(jmax-2);     /**< sparse matrix size component (n,n), boundaries excluded */

    Mesh::MasterElement mastElem = Mesh::MasterElement(4,1);
    mastElem.setGLLOrder(0,7);
    mastElem.setGLLOrder(1,7);
    mastElem.setGLLOrder(2,7);
    mastElem.setGLLOrder(3,7);

    f64 x    = -0.9999;
    f64 d0 = Polynomials::  Legendre(3, x);
    f64 d1 = Polynomials::d1Legendre(3, x);
    f64 d2 = Polynomials::d2Legendre(3, x);
    f64 d3 = Polynomials::d3Legendre(3, x);

    INFO_MSG("%-4f %-4f %-4f %-4f %-4f", x, d0, d1, d2, d3)


    //## ============= ##//
    //## Problem Setup ##//
    //## ============= ##//
    // Declare and initalise 2D gridpoint list
    // Mesh::gridStruct grid;
    // grid.x.setLinSpaced(imax, Lx[0], Lx[1]);
    // grid.y.setLinSpaced(jmax, Ly[0], Ly[1]);

    //Mesh::Simplex simplex = Mesh::Simplex(1);
    //Mesh::Cube    cube    = Mesh::Cube(3,"SEM");

    // //## =============== ##//
    // //## Export solution ##//
    // //## =============== ##//
    // std::ofstream dataFile("data.bin", std::ios::out | std::ios::binary | std::ios::trunc); /**< Data output file, not using std::ios::app */ 
    // dataFile.write(reinterpret_cast<const char*>(&jmax), sizeof(jmax));
    // dataFile.write(reinterpret_cast<const char*>(&imax), sizeof(imax));

    // u32 idx = 0;
    // for (u32 j=0; j<jmax; j++){
    //     for (u32 i=0; i<imax; i++){
    //         // Plotting / postprocessing currently does not need to be done in such high precision
    //         f32 xPoint = (float) grid.x[i]; /**< f32 x-position */
    //         f32 yPoint = (float) grid.y[j]; /**< f32 y-position */
    //         f32 uPoint;                     /**< f32 u-solution */
    //         if      (i==0     ) {uPoint = (float) boundaries.West[j]; }
    //         else if (i==imax-1) {uPoint = (float) boundaries.East[j]; }
    //         else if (j==0     ) {uPoint = (float) boundaries.South[i];}
    //         else if (j==jmax-1) {uPoint = (float) boundaries.North[i];}
    //         else                {uPoint = (float) u[idx]; idx++;      }

    //         // INFO_MSG("idx=%d, uPoint=%f", idx, uPoint)
    //         dataFile.write(reinterpret_cast<const char*>(&xPoint), sizeof(xPoint));
    //         dataFile.write(reinterpret_cast<const char*>(&yPoint), sizeof(yPoint));
    //         dataFile.write(reinterpret_cast<const char*>(&uPoint), sizeof(uPoint));

    //     }
    // }
    // dataFile.close();

    // INFO_MSG("Solution saved.");

    return EXIT_SUCCESS;
}
