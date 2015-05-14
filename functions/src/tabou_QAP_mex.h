/*****************************************************************/
// Implementation of the robust taboo search of: E. Taillard
// "Robust taboo search for the quadratic assignment problem", 
// Parallel Computing 17, 1991, 443-455.
//
// Data file format: 
//  n,
// (nxn) flow matrix,
// (nxn) distance matrix
//
// Copyright : E. Taillard, 1990-2004
// This code can be freely used for non-commercial purpose.
// Any use of this implementation or a modification of the code
// must acknowledge the work of E. Taillard
/****************************************************************/

#ifndef TABOUQAP
#define TABOUQAP
#include "mex.h"

/**
 * Convert a 1-D array to a type_matrix object
 */
int64_T** array_to_type_matrix(int64_T* arr, int64_T Q);

/**
 * Destroy (free) a type_matrix object
 */
void destroy_type_matrix(int64_T** matrix, int64_T Q);

/**
 * Generate random solution to the QAP problem
 */
void generate_random_solution(int64_T n, int64_T* sol);

/**
 * Swap 2 int variables
 */
void transpose(int64_T& a, int64_T& b);

/**
 * Return an integer between low and high 
 */
int64_T unif(int64_T low, int64_T high);

/**
 * L'Ecuyer random number generator
 */
const int64_T m = 2147483647;
const int64_T m2 = 2145483479;

const int64_T a12 = 63308;
const int64_T q12 = 33921;
const int64_T r12 = 12979;

const int64_T a13 = -183326;
const int64_T q13 = 11714;
const int64_T r13 = 2883;

const int64_T a21 = 86098;
const int64_T q21 = 24919;
const int64_T r21 = 7417;

const int64_T a23 = -539608;
const int64_T q23 = 3976;
const int64_T r23 = 2071;

const double invm = 4.656612873077393e-10;

int64_T x10 = 12345, x11 = 67890, x12 = 13579, 
        x20 = 24680, x21 = 98765, x22 = 43210; // initializae the values of seeds

double rando();

#endif