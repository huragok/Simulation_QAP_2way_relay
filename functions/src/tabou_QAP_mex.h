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
int64_T** array_to_type_matrix(int64_T* arr, int Q);

/**
 * Destroy (free) a type_matrix object
 */
void destroy_type_matrix(int64_T** matrix, int Q);

#endif