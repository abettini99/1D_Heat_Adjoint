#pragma once

// ------------------ //
// unsigned int types //
// ------------------ // 
/** 8-bit unsigned integer */
using u8  = unsigned char;
/** 16-bit unsigned integer */
using u16 = unsigned short;
/** 32-bit unsigned integer */
using u32 = unsigned int;
/** 64-bit unsigned integer */
using u64 = unsigned long long;



// ---------------- //
// signed int types //
// ---------------- // 
/** 8-bit signed integer */
using i8  = signed char;
/** 16-bit signed integer */
using i16 = signed short;
/** 32-bit signed integer */
using i32 = signed int;
/** 64-bit signed integer */
using i64 = signed long long;



// -------------------- //
// floating point types //
// -------------------- // 
/** 32-bit floating point */
using f32 = float;
/** 64-bit floating point */
using f64 = double;



// ------------- //
// boolean types //
// ------------- // 
/** 32-bit boolean */
using b32 = int;
/** 8-bit boolean */
using b8 =  char;



// ------- //
// boolean //
// ------- // 
/** true boolean (in text form) */
#define TRUE  1
/** false boolean (in text form) */
#define FALSE 0



// ---------- //
// exit codes //
// ---------- // 
/** exit code upon success */
#define EXIT_SUCCESS          0
/** exit code upon unknown failure */
#define EXIT_FAILURE_UNKNOWN -1
/** exit code upon assertion failure */
#define EXIT_FAILURE_ASSERTION 1
/** exit code upon iteration failure */
#define EXIT_FAILURE_ITERATION 2


// ----------------- //
// static assertions //
// ----------------- // 
// define static assertions
#if defined(__clang__) || defined(__gcc__)
/** compile time assertion (clang, gcc) */
#define STATIC_ASSERT _Static_assert
#else
/** compile time assertion */
#define STATIC_ASSERT static_assert
#endif

// check all types are correctly size
STATIC_ASSERT(sizeof(u8)  == 1,  "u8 should be 1 byte.");
STATIC_ASSERT(sizeof(u16) == 2, "u16 should be 2 byte.");
STATIC_ASSERT(sizeof(u32) == 4, "u32 should be 4 byte.");
STATIC_ASSERT(sizeof(u64) == 8, "u64 should be 8 byte.");

STATIC_ASSERT(sizeof(i8)  == 1,  "i8 should be 1 byte.");
STATIC_ASSERT(sizeof(i16) == 2, "i16 should be 2 byte.");
STATIC_ASSERT(sizeof(i32) == 4, "i32 should be 4 byte.");
STATIC_ASSERT(sizeof(i64) == 8, "i64 should be 8 byte.");

STATIC_ASSERT(sizeof(f32) == 4, "f32 should be 4 byte.");
STATIC_ASSERT(sizeof(f64) == 8, "f64 should be 8 byte.");
