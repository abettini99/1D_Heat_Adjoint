#pragma once

#include "definesStandard.hpp"
#include "Eigen/Sparse"

/************************************************************************************************************************ 
 *  @brief (Re)defined Eigen's typedefs are placed into this namespace.
 * 
 *  @details
 *  All personalised typedefs (for Cpp, "typedef" is "using") related to the Eigen library are placed in this namespace.
 ************************************************************************************************************************/
namespace EigenDefs{

// ------------ //
// matrix types //
// ------------ // 
/************************************************************************************************************************ 
 *  @brief Dynamically-sized matrix of type Type, e.g. i32, f64, i8, etc... use when size of matrix (n > 4, m > 4)
 * 
 *  @details
 *  References about matrix types and optimization (why there is Matrix22, Matrix33, Matrix44):
 *  @link https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
 *  @link https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html#title10
 *  @link https://eigen.tuxfamily.org/dox/group__matrixtypedefs.html#gaf20e523ca57ee8ef0a945cd4703d2bfd
 ************************************************************************************************************************/
template<typename Type> using Matrix   = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>;
/** Compile-time optimized matrix of size (2,2) of type Type, e.g. i32, f64, i8, etc... */
template<typename Type> using Matrix22 = Eigen::Matrix<Type, 2, 2>;
/** Compile-time optimized matrix of size (3,3) of type Type, e.g. i32, f64, i8, etc... */
template<typename Type> using Matrix33 = Eigen::Matrix<Type, 3, 3>;
/** Compile-time optimized matrix of size (4,4) of type Type, e.g. i32, f64, i8, etc... */
template<typename Type> using Matrix44 = Eigen::Matrix<Type, 4, 4>;



// ------------ //
// vector types //
// ------------ // 
/************************************************************************************************************************ 
 *  @brief Dynamically-sized vector of type Type, e.g. i32, f64, i8, etc... use when size of vector > 4
 * 
 *  @details
 *  References about vector types and optimization (why there is Vector2, Vector3, Vector4):
 *  @link https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
 *  @link https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html#title10
 *  @link https://eigen.tuxfamily.org/dox/group__matrixtypedefs.html#gaf20e523ca57ee8ef0a945cd4703d2bfd
 ************************************************************************************************************************/
template<typename Type> using Vector  = Eigen::Matrix<Type, Eigen::Dynamic, 1>;
/** Compile-time optimized vector of size 2 of type Type, e.g. i32, f64, i8, etc... */
template<typename Type> using Vector2 = Eigen::Matrix<Type, 2, 1>;
/** Compile-time optimized vector of size 3 of type Type, e.g. i32, f64, i8, etc... */
template<typename Type> using Vector3 = Eigen::Matrix<Type, 3, 1>;
/** Compile-time optimized vector of size 4 of type Type, e.g. i32, f64, i8, etc... */
template<typename Type> using Vector4 = Eigen::Matrix<Type, 4, 1>;



// ----------- //
// array types //
// ----------- // 
/** Dynamically-sized 1D array of type Type, e.g. i32, f64, i8, etc... */
template<typename Type> using Array1D = Eigen::Array<Type, Eigen::Dynamic, 1>;
/** Dynamically-sized 2D array of type Type, e.g. i32, f64, i8, etc... */
template<typename Type> using Array2D = Eigen::Array<Type, Eigen::Dynamic, Eigen::Dynamic>;


} // namespace EigenDefs