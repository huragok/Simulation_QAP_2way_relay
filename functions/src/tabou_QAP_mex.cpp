#include "tabou_QAP_mex.h"
using namespace std;

/**
 * The gateway routine that for the tabu search QAP solver.
 *
 * @param nlhs Number of output arguments, at most 2
 * @param plhs Pointer array to the output arguments: plhs[0] is 1-by-Q mxArray (map) and plhs[0] is a 1-by-1 mxArray (cost). Both should be int64_T
 * @param nrhs Number of input arguments, must be 3
 * @param prhs Pointer array to the input arguments: prhs[0] is a Q-by-Q mxArray (flow matrix), prhs[1] is a Q-by-Q mxArray (distance matrix) and prhs[2] is a 1-by-1 mxArray (number of iterations). All should be of type int64_T
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check the number of arguments
    if (nrhs != 3)
    {
        mexErrMsgTxt("Wrong number of input arguments.");
    }
    else if (nlhs > 2)
    {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    // Check the type of input arguments
    if (!mxIsInt64(prhs[0]) || !mxIsInt64(prhs[1]) || !mxIsInt64(prhs[2]))
    {
        mexErrMsgTxt("All 3 input arguments must be of type int64.");
    }
    
    // Read and check the size of the input matrices 
    int64_T Q = mxGetM(prhs[0]);
    if (Q != mxGetN(prhs[0]) || Q != mxGetM(prhs[1]) || Q != mxGetN(prhs[1]))
    {
        mexErrMsgTxt("Wrong input matrix size: both the flow and distance matrix must be of size Q-by-Q.");
    }
    
    // Convert the flow and distance matrix from type mxArray to type_matrix
    int64_T** flow = array_to_type_matrix((int64_T*)mxGetData(prhs[0]), Q);
    int64_T** dist = array_to_type_matrix((int64_T*)mxGetData(prhs[1]), Q);
    
    // Read the number of iterations
    int64_T nitr = *(int64_T*)mxGetData(prhs[2]);
    
    // Create the ouput arguments
    plhs[0] = mxCreateNumericMatrix(1, Q, mxINT64_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);    
    int64_T* map = (int64_T*)mxGetData(plhs[0]);
    int64_T* cost = (int64_T*)mxGetData(plhs[1]);
    int64_T cost_tmp;
    
    // Generate the random solution
    int64_T* solution = (int64_T*)mxMalloc((Q + 1) * sizeof(int64_T));
    
    generate_random_solution(Q, solution);
    
// Call the worker functions to perform the tabu search
    //for (int64_T k = 1; k <= Q; k++)
    //{
    //    map[k] = solution[k];
    //}
    //*cost = x10;
    
    tabu_search(Q, 
                flow, dist,
                solution, cost_tmp,
                9 * Q / 10, 11 * Q / 10, 
                Q * Q * 2,
                nitr);
    
    *cost = cost_tmp;             
    // Free memories
    mxFree(solution);
    destroy_type_matrix(flow, Q);
    destroy_type_matrix(dist, Q);
    return;
}

int64_T** array_to_type_matrix(int64_T* arr, int64_T Q)
{
    int64_T** matrix = (int64_T**)mxMalloc(Q * sizeof(int64_T*));
    for (int64_T i = 0; i < Q; i++)
    {
        matrix[i] = (int64_T*)mxMalloc(Q * sizeof(int64_T));
    }
    
    for (int64_T r = 0; r < Q; r++)
    {
        for (int64_T c = 0; c < Q; c++)
        {
            matrix[r][c] = arr[c * Q + r];
        }
    }
    return(matrix);
}

void destroy_type_matrix(int64_T** matrix, int64_T Q)
{
    for (int64_T i = 0; i < Q; i++)
    {
        mxFree(matrix[i]);
    }
    mxFree(matrix);
    return;
}

void generate_random_solution(int64_T Q, int64_T* sol)
{
	int64_T i;
	for (i = 0; i <= Q; i++)
    {
        sol[i] = i;
    }
	for (i = 1; i < Q; i++)
    {
        transpose(sol[i], sol[unif(i, Q)]);
    }
    return;
}

void transpose(int64_T& a, int64_T& b)
{
    int64_T temp = a;
    a = b;
    b = temp;
    return;
}

int64_T unif(int64_T low, int64_T high)
{
    return(low + int64_T(double(high - low + 1) * rando()));
}

double rando()
{
    int64_T h, p12, p13, p21, p23;
    
	h = x10 / q13;
    p13 = -a13 * (x10 - h * q13) - h * r13;
	h = x11 / q12;
    p12 = a12 * (x11 - h * q12) - h * r12;  
	if (p13 < 0) p13 = p13 + m;
    if (p12 < 0) p12 = p12 + m;
	
    x10 = x11;
    x11 = x12;
    x12 = p12 - p13;
    if (x12 < 0) x12 = x12 + m;
    
	h = x20 / q23;
    p23 = -a23 * (x20 - h * q23) - h * r23;
	h = x22 / q21;
    p21 = a21 * (x22 - h * q21) - h * r21;
	if (p23 < 0) p23 = p23 + m2;
    if (p21 < 0) p21 = p21 + m2;
    
	x20 = x21;
    x21 = x22;
    x22 = p21-p23;
    if(x22 < 0) x22 = x22 + m2;
  
    if (x12 < x22)
    {
        h = x12 - x22 + m;
    }
    else
    {
        h = x12 - x22;
    }
    
	if (h == 0)
    {
        return(1.0);
    }
    else
    {
        return(h * invm);
    }
 }

void tabu_search(int64_T n, 
                 int64_T** a, int64_T** b,
                 int64_T* best_sol, int64_T& best_cost,
                 int64_T min_size, int64_T max_size,
                 int64_T aspiration,
                 int64_T nr_iterations)
{
    int64_T* p;                        // current solution
	int64_T** delta;                    // store move costs
	int64_T** tabu_list;                // tabu status
	int64_T current_iteration;               // current iteration
	int64_T current_cost;                    // current sol. value
	int64_T i, j, k, i_retained, j_retained;  // indices
    
    /***************** dynamic memory allocation *******************/
	p = (int64_T*)mxMalloc((n + 1) * sizeof(int64_T));
	delta = (int64_T**)mxMalloc((n + 1) * sizeof(int64_T*));
	for (i = 1; i <= n; i++)
    {
        delta[i] = (int64_T*)mxMalloc((n + 1) * sizeof(int64_T));
    }
	tabu_list = (int64_T**)mxMalloc((n + 1) * sizeof(int64_T*));
	for (i = 1; i <= n; i++)
    {
        tabu_list[i] = (int64_T*)mxMalloc((n + 1) * sizeof(int64_T));
    }
    
     /***************** dynamic memory release *******************/
    mxFree(p);
    for (i=1; i <= n; i++)
    {
        mxFree(tabu_list[i]);
    }
    mxFree(tabu_list);
    for (i=1; i <= n; i++)
    {
        mxFree(delta[i]);
    }
    mxFree(delta);
    
    return;
}