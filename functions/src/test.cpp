#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mxClassID category;
    category = mxGetClassID(prhs[0]);
    switch (category)
    {
        case mxINT8_CLASS:   mexPrintf("a");   break; 
        case mxUINT8_CLASS:  mexPrintf("b");  break;
        case mxINT16_CLASS:  mexPrintf("c");  break;
        case mxUINT16_CLASS: mexPrintf("d"); break;
        case mxINT32_CLASS:  mexPrintf("e");  break;
        case mxUINT32_CLASS: mexPrintf("f"); break;
        case mxINT64_CLASS:  mexPrintf("g");  break;
        case mxUINT64_CLASS: mexPrintf("h"); break;
        case mxSINGLE_CLASS: mexPrintf("i"); break; 
        case mxDOUBLE_CLASS: mexPrintf("j"); break;
        default: break;
	}
    int64_T* input = (int64_T *)mxGetData(prhs[0]); // Size of the problem
    //mexPrintf("%l", *input);
 
    int64_T M = mxGetM(prhs[0]);
    int64_T N = mxGetN(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
    double* output = mxGetPr(plhs[0]);
    output[0] = (double)(M);
    output[1] = (double)(N);
    return;
}