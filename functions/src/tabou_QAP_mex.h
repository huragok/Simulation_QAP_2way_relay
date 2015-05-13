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

#include <fstream>
using namespace std;

const long infinite = 999999999;

typedef int*   type_vector;
typedef long** type_matrix;

/*************** L'Ecuyer random number generator ***************/
const long m = 2147483647; const long m2 = 2145483479; 
const long a12= 63308; const long q12=33921; const long r12=12979; 
const long a13=-183326; const long q13=11714; const long r13=2883; 
const long a21= 86098; const long q21=24919; const long r21= 7417; 
const long a23=-539608; const long q23= 3976; const long r23=2071;
const double invm = 4.656612873077393e-10;

long x10 = 12345, x11 = 67890, x12 = 13579, // init. de la
     x20 = 24680, x21 = 98765, x22 = 43210; // valeur des germes

double rando();

/*********** return an integer between low and high *************/
long unif(long low, long high);

void transpose(int & a, int & b);

int min(long a, long b);

/*--------------------------------------------------------------*/
/*       compute the cost difference if elements i and j        */
/*         are transposed in permutation (solution) p           */
/*--------------------------------------------------------------*/
long compute_delta(int n, type_matrix & a, type_matrix & b, type_vector & p, int i, int j);

/*--------------------------------------------------------------*/
/*      Idem, but the value of delta[i][j] is supposed to       */
/*    be known before the transposition of elements r and s     */
/*--------------------------------------------------------------*/
long compute_delta_part(type_matrix & a, type_matrix & b, type_vector & p, type_matrix & delta, int i, int j, int r, int s);

/*
 * The worker function for tabou search
 */
void tabu_search(long n,                  // problem size
                 type_matrix & a,         // flows matrix
                 type_matrix & b,         // distance matrix
                 type_vector & best_sol,  // best solution found
                 long & best_cost,        // cost of best solution
                 long min_size,           // parameter 1 (< n^2/2)
                 long max_size,           // parameter 2 (< n^2/2)
                 long aspiration,         // parameter 3 (> n^2/2)
                 long nr_iterations)      // number of iterations 
                 
void generate_random_solution(long n, type_vector& p);

/*
 * Convert an array to a type_matrix
 */
type_matrix array_to_type_matrix(long* arr, long n);
#endif