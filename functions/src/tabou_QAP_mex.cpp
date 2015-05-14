#include "tabou_QAP_mex.h"
using namespace std;

/**
 * The gateway routine that for the tabu search QAP solver.
 *
 * @param nlhs Number of output arguments, at most 2
 * @param plhs Pointer array to the output arguments: plhs[0] is 1-by-Q mxArray (map) and plhs[0] is a 1-by-1 mxArray (cost). plhs[0] should be int16_T and plhs[1] should be int64_T
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
    int Q = mxGetM(prhs[0]);
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
    plhs[0] = mxCreateNumericMatrix(1, Q, mxINT16_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);    
    size_t* map = (size_t*)mxGetData(plhs[0]);
    int64_T* cost = (int64_T*)mxGetData(plhs[1]);
    
    // Free memories
    destroy_type_matrix(flow, Q);
    destroy_type_matrix(dist, Q);
    return;
}

int64_T** array_to_type_matrix(int64_T* arr, int Q)
{
    int64_T** matrix = (int64_T**)mxMalloc(Q * sizeof(int64_T*));
    for (int i = 0; i < Q; i++)
    {
        matrix[i] = (int64_T*)mxMalloc(Q * sizeof(int64_T));
    }
    
    for (int r = 0; r < Q; r++)
    {
        for (int c = 0; c < Q; c++)
        {
            matrix[r][c] = arr[c * Q + r];
        }
    }
    return(matrix);
}

void destroy_type_matrix(int64_T** matrix, int Q)
{
    for (size_t i = 0; i < Q; i++)
    {
        mxFree(matrix[i]);
    }
    mxFree(matrix);
    return;
}