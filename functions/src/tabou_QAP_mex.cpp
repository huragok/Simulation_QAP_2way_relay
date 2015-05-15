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
    
    // Return the results to matlab
    *cost = cost_tmp;
    for (int64_T i = 0; i < Q; i++)
    {
        map[i] = solution[i + 1];
    }
    
    // Free memories
    mxFree(solution);
    destroy_type_matrix(flow, Q);
    destroy_type_matrix(dist, Q);
    return;
}

int64_T** array_to_type_matrix(int64_T* arr, int64_T Q)
{
    int64_T** matrix = (int64_T**)mxMalloc((Q + 1) * sizeof(int64_T*));
    for (int64_T i = 1; i <= Q; i++)
    {
        matrix[i] = (int64_T*)mxMalloc((Q + 1) * sizeof(int64_T));
    }
    
    for (int64_T r = 1; r <= Q; r++)
    {
        for (int64_T c = 1; c <= Q; c++)
        {
            matrix[r][c] = arr[(c - 1) * Q + (r - 1)];
        }
    }
    return(matrix);
}

void destroy_type_matrix(int64_T** matrix, int64_T Q)
{
    for (int64_T i = 1; i <= Q; i++)
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

int64_T min(int64_T a, int64_T b)
{
    if (a < b)
    {
        return(a);
    }
    else
    {
        return(b);
    }
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

int64_T compute_delta(int64_T n, int64_T** a, int64_T** b, int64_T* p, int64_T i, int64_T j)
{
    int64_T d, k;
	d = (a[i][i] - a[j][j]) * (b[p[j]][p[j]] - b[p[i]][p[i]]) +
        (a[i][j] - a[j][i]) * (b[p[j]][p[i]] - b[p[i]][p[j]]);
	for (k = 1; k <= n; k++)
    {
        if (k!=i && k!=j)
        {
            d += (a[k][i] - a[k][j]) * (b[p[k]][p[j]] - b[p[k]][p[i]]) +
                 (a[i][k] - a[j][k]) * (b[p[j]][p[k]] - b[p[i]][p[k]]);
        }
    }
	return(d);
}

int64_T compute_delta_part(int64_T** a, int64_T** b, int64_T* p, int64_T** delta, int64_T i, int64_T j, int64_T r, int64_T s)
{
    return(delta[i][j] + 
           (a[r][i]-a[r][j]+a[s][j]-a[s][i]) * (b[p[s]][p[i]]-b[p[s]][p[j]]+b[p[r]][p[j]]-b[p[r]][p[i]]) +
           (a[i][r]-a[j][r]+a[j][s]-a[i][s]) * (b[p[i]][p[s]]-b[p[j]][p[s]]+b[p[j]][p[r]]-b[p[i]][p[r]]));
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
    
    /************** current solution initialization ****************/
	for (i = 1; i <= n; i++)
    {
        p[i] = best_sol[i];
    }
    
    /********** initialization of current solution value ***********/
	/**************** and matrix of cost of moves  *****************/
	current_cost = 0;
	for (i = 1; i <= n; i = i + 1)
    {
        for (j = 1; j <= n; j = j + 1)
        {
            current_cost = current_cost + a[i][j] * b[p[i]][p[j]];
            if (i < j)
            {
                delta[i][j] = compute_delta(n, a, b, p, i, j);
            }
        }
	}
	best_cost = current_cost;
    
	/****************** tabu list initialization *******************/
	for (i = 1; i <= n; i++) 
    {
        for (j = 1; j <= n; j++)
        {
            tabu_list[i][j] = -(n * i + j);
        }
    }
    
    /******************** main tabu search loop ********************/
	for (current_iteration = 1; current_iteration <= nr_iterations; current_iteration++)
	{
        /** find best move (i_retained, j_retained) **/
        i_retained = infinite;       // in case all moves are tabu
        int64_T min_delta = infinite;   // retained move cost
        bool autorized;               // move not tabu?
        bool aspired;                 // move forced?
        bool already_aspired = false; // in case many moves forced

        for (i = 1; i < n; i++)
        {
            for (j = i+1; j <= n; j++)
            {
                autorized = (tabu_list[i][p[j]] < current_iteration) || 
                            (tabu_list[j][p[i]] < current_iteration);

                aspired = (tabu_list[i][p[j]] < current_iteration-aspiration)||
                          (tabu_list[j][p[i]] < current_iteration-aspiration)||
                          (current_cost + delta[i][j] < best_cost);                

                if ((aspired && !already_aspired) || // first move aspired
                    (aspired && already_aspired && (delta[i][j] < min_delta)) || // many move aspired and take best one
                    (!aspired && !already_aspired && (delta[i][j] < min_delta) && autorized)) // no move aspired yet
                {
                    i_retained = i;
                    j_retained = j;
                    min_delta = delta[i][j];
                    if (aspired) 
                    {
                        already_aspired = true;
                    };
                }
            }
        }

        if (i_retained == infinite)
        {
            mexWarnMsgTxt("All moves are tabu!");
        }
        else 
        {
            /** transpose elements in pos. i_retained and j_retained **/
            transpose(p[i_retained], p[j_retained]);
            // update solution value
            current_cost += delta[i_retained][j_retained];
            // forbid reverse move for a random number of iterations
            tabu_list[i_retained][p[j_retained]] = current_iteration + unif(min_size,max_size);
            tabu_list[j_retained][p[i_retained]] = current_iteration + unif(min_size,max_size);

            // best solution improved ?
            if (current_cost < best_cost)
            {
                best_cost = current_cost;
                for (k = 1; k <= n; k++)
                {
                    best_sol[k] = p[k];
                }
                // Uncomment the next line to prompt to matlab command window when ever a better solution is achieved.
                // mexPrintf("Solution of the value %f found at iter. %f\n", (double)(best_cost), (double)(current_iteration));
            }

            // update matrix of the move costs
            for (i = 1; i < n; i++)
            {
                for (j = i+1; j <= n; j++)
                {
                    if (i != i_retained && i != j_retained && j != i_retained && j != j_retained)
                    {
                        delta[i][j] = compute_delta_part(a, b, p, delta, i, j, i_retained, j_retained);
                    }
                    else
                    {
                        delta[i][j] = compute_delta(n, a, b, p, i, j);
                    }
                }
            }
        }
	}
   
    /***************** dynamic memory release *******************/
    mxFree(p);
    destroy_type_matrix(tabu_list, n);
    destroy_type_matrix(delta, n);
    
    return;
}