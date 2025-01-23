#include "CoreIncludes.hpp"
#include "mesh.hpp"

namespace Mesh{

MasterElement::MasterElement(u8 nVars_, u8 nDims_) : nVars(nVars_), nDims(nDims_) {

    CHECK_FATAL_ASSERT(nVars > 0, "Number of variables must be bigger than 0")
    CHECK_FATAL_ASSERT(nDims > 0, "Number of dimensions must be bigger than 0")

}

} // end Mesh
