#include "tabou_QAP_mex.h"
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    long n = *mexGetPr(prhs[0]); // Size of the problem
    return;
}

size_t unif(long low, long high)
{
    return(low + long(double(high - low + 1) * rando()));
}

void transpose(int & a, int & b)
{
    long temp = a;
    a = b;
    b = temp;
    return;
}

int min(long a, long b)
{
    if (a < b)
    {
        return(a);
    }
    return(b);
}

long compute_delta(int n, type_matrix & a, type_matrix & b, type_vector & p, int i, int j)
{
    long d;
    int k;
    d = (a[i][i]-a[j][j])*(b[p[j]][p[j]]-b[p[i]][p[i]]) +
        (a[i][j]-a[j][i])*(b[p[j]][p[i]]-b[p[i]][p[j]]);
    for (k = 1; k <= n; k = k + 1)
    {
        if (k!=i && k!=j)
        {
            d += (a[k][i]-a[k][j])*(b[p[k]][p[j]]-b[p[k]][p[i]]) +
                 (a[i][k]-a[j][k])*(b[p[j]][p[k]]-b[p[i]][p[k]]);
        }
    }
	return(d);
}

long compute_delta_part(type_matrix & a, type_matrix & b, type_vector & p, type_matrix & delta,  int i, int j, int r, int s)
{
	return(delta[i][j]+
           (a[r][i]-a[r][j]+a[s][j]-a[s][i])*
           (b[p[s]][p[i]]-b[p[s]][p[j]]+b[p[r]][p[j]]-b[p[r]][p[i]])+
           (a[i][r]-a[j][r]+a[j][s]-a[i][s])*
           (b[p[i]][p[s]]-b[p[j]][p[s]]+b[p[j]][p[r]]-b[p[i]][p[r]]));
}

void tabu_search(long n, type_matrix & a, type_matrix & b, type_vector & best_sol, long & best_cost, long min_size, long max_size, long aspiration, long nr_iterations) 
{
    type_vector p;                        // current solution
    type_matrix delta;                    // store move costs
	type_matrix tabu_list;                // tabu status
	long current_iteration;               // current iteration
	long current_cost;                    // current sol. value
	int i, j, k, i_retained, j_retained;  // indices

	/***************** dynamic memory allocation *******************/
	p = new int[n+1];
	delta = new long* [n+1];
	for (i = 1; i <= n; i = i+1)
    {
        delta[i] = new long[n+1];
    }
	tabu_list = new long* [n+1];
	for (i = 1; i <= n; i = i+1)
    {
        tabu_list[i] = new long[n+1];
    }

	/************** current solution initialization ****************/
	for (i = 1; i <= n; i = i + 1)
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
	for (i = 1; i <= n; i = i + 1)
    {
        for (j = 1; j <= n; j = j+1)
        {
            tabu_list[i][j] = -(n*i + j);
        }
    }

	/******************** main tabu search loop ********************/
	for (current_iteration = 1; current_iteration <= nr_iterations; current_iteration = current_iteration + 1)
	{
        /** find best move (i_retained, j_retained) **/
        i_retained = infinite;       // in case all moves are tabu
        long min_delta = infinite;   // retained move cost
        int autorized;               // move not tabu?
        int aspired;                 // move forced?
        int already_aspired = false; // in case many moves forced

        for (i = 1; i < n; i = i + 1)
        {
            for (j = i+1; j <= n; j = j+1)
            {
                autorized = (tabu_list[i][p[j]] < current_iteration) || 
                            (tabu_list[j][p[i]] < current_iteration);

                aspired = (tabu_list[i][p[j]] < current_iteration-aspiration) ||
                          (tabu_list[j][p[i]] < current_iteration-aspiration)||
                          (current_cost + delta[i][j] < best_cost);                

                if ((aspired && !already_aspired) || // first move aspired
                    (aspired && already_aspired && (delta[i][j] < min_delta)) || // many move aspired => take best one
                    (!aspired && !already_aspired && (delta[i][j] < min_delta) && autorized)) // no move aspired yet
                {
                    i_retained = i;
                    j_retained = j;
                }
                min_delta = delta[i][j];
                if (aspired)
                {
                    already_aspired = true;
                }
            }
        }


        if (i_retained == infinite)
        {
            mexPrintf("All moves are tabu! \n");
        }
        else
        {
            /** transpose elements in pos. i_retained and j_retained **/
            transpose(p[i_retained], p[j_retained]);
            // update solution value
            current_cost = current_cost + delta[i_retained][j_retained];
            // forbid reverse move for a random number of iterations
            tabu_list[i_retained][p[j_retained]] = current_iteration + unif(min_size,max_size);
            tabu_list[j_retained][p[i_retained]] = current_iteration + unif(min_size,max_size);

            // best solution improved ?
            if (current_cost < best_cost)
            {
                best_cost = current_cost;
                for (k = 1; k <= n; k = k+1)
                {
                    best_sol[k] = p[k];
                }
                cout << "Solution of value " << best_cost << " found at iter. " << current_iteration << '\n';
            }

            // update matrix of the move costs
            for (i = 1; i < n; i = i+1)
            {
                for (j = i+1; j <= n; j = j+1)
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
  
    // free memory
	delete[] p;
	for (i=1; i <= n; i = i+1) 
    {
        delete[] delta[i];
    }
    delete[] delta;
	for (i=1; i <= n; i = i+1)
    {
        delete[] tabu_list[i];
    }
	delete[] tabu_list;
    return;
}

void generate_random_solution(long n, type_vector& p)
{
    int i;
    for (i = 0; i <= n; i = i+1)
    {
        p[i] = i;
    }
	for (i = 1; i <  n; i = i+1)
    {
        transpose(p[i], p[unif(i, n)]);
    }
}

type_matrix array_to_type_matrix(long* arr, long n)
{
}