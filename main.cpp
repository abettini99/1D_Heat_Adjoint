/************************************************************************************************************************
 * In these two header files, we include:
 *    - Typedefs:   used to create an additional name for another data type, but does not create a new type.
 *    - Type alias: same as typedefs, but allows for templates, e.g. Vector<f64>, Vector<i32>.
 ************************************************************************************************************************/
#include "CoreIncludes.hpp"
#include "mesh/mesh.hpp"
#include "mesh/polynomials.hpp"

#include <mpi.h>
#include <vector>
#include <fstream>
#include <iostream>

/************************************************************************************************************************
 * Solve -div(grad(u)) = f, using FDM
 ************************************************************************************************************************/
int main(int argc, char *argv[]){
    i32 rankid, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (rankid == 0) {
        //## ================== ##//
        //## Provide parameters ##//
        //## ================== ##//
        const u32 nElemx  = 10;                /**< #Elements in x */
        const f64 Lx[2]   = {0., 1.*EIGEN_PI}; /**< domain endpoints in x */

        EigenDefs::Array1D<f64> x1e = EigenDefs::Array1D<f64>::LinSpaced(nElemx+1, Lx[0], Lx[1]); /**< x1 Endpoints */
        Mesh::Geometry Domain = Mesh::Geometry(x1e, x1e, x1e);
        Domain.MasterElement.setnVars(1);
        Domain.MasterElement.setLGLOrder(0,7);

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
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
