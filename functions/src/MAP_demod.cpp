/** 
 * file: Map_demod.cpp
 * description: function of Chase combining HARQ MAP detector, and its mex function 
 * edited by: Wenhao Wu
 * created date: 05/26/15
 *
 * replaced file: Map_demod.cpp
 * description: function of MIMO MAP detector, and its mex function 
 * author: Weiliang Zeng
 * created date: 03/01/10
 * modified: 
 *	on 03/09/10: add accurate and piece-wise approaximate method to calculate log(e^x+e^y)
 */

#include "mex.h"
#include <complex>
#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>

using namespace std;

#define INF  10000  // infinity

/**
 * Inner product of two double vecters.
 * @param vec_a 1-by-L double array, the first vector
 * @param vec_b 1-by-L double array, the second vector
 * @param L int scalar, the length of the vectors
 * @return a double scalar of the inner product of vec_a and vec_b
 */
inline double vec_sum_mul(double* vec_a,  double* vec_b, int L)
{
	double tmp = 0;
	for(int i = 0; i < L; i++)
	{
		tmp = tmp + vec_a[i] * vec_b[i]; 
	}
	return tmp;
}

/**
 * Jacbian logarithm. 
 * correction function r(x) = log(1+e^-x), use piease wise linear function to approximate r(x) as table II in [1]
 * @param x double scalar
 * @param y double scalar
 * @return log(e^x + e^y) = max(x,y) + log(1+e^(-|x-y|))
 * @see [1] Xiao-Yu Hu,Evangelos Eleftheriou,Dieter-Michael Arnold,and Ajay Dholakia. "Efficient Implementations of the Sum-Product Algorithm for Decoding LDPC Codes". GlobalCom , 2001:ppV2 1036-1036E
 */
inline double jacln(double x, double y)
{
	// find larger one of x and y
	double max_input, min_input;
	if (x >= y)
	{  
		max_input = x;
		min_input = y;
	}
	else
	{
		max_input = y;
		min_input = x;
	}

	// find correction r(max_input - min_input)
   double corr_input, corr_output, result;
   corr_input = max_input - min_input;
	if ((corr_input >= 0) && (corr_input < 0.5))
	{
		corr_output = -corr_input * 0.5 + 0.7;
	}
	else if ((corr_input >= 0.5) && (corr_input < 1.6))
	{
		corr_output = -corr_input * 0.25 + 0.575;
	}
	else if ((corr_input >= 1.6) && (corr_input < 2.2))
	{
		corr_output = -corr_input * 0.125 + 0.375;
	}
	else if ((corr_input >= 2.2) && (corr_input < 3.2))
	{
		corr_output = -corr_input * 0.0625 + 0.2375;
	}
	else if ((corr_input >= 3.2) && (corr_input < 4.4))
	{
		corr_output = -corr_input * 0.03125 + 0.1375;
	}
	else
	{
		corr_output = 0;
	}	

	// log(e^x + e^y) = max(x,y) + log(1+e^(-|x-y|)).
	result = max_input + corr_output;
    return (result);
}

/**
 * The worker function of Chase combining HARQ MAP detector. 
 * @param LextDemodulation 1-by-(Nbps * n_symbol) double matrix, the extrinsic LLR 
 * @param rx_signal M-by-n_symbol complex<double> matrix, the received signal corresponding to 1 LDPC frame and M (re)transmisstions
 * @param chnl_eq M-by-n_symbol complex<double> matrix, equivalent channel when the noises or the M (re)transmissions are normalized to have equal power
 * @param bit_mat_anti Q-by-Nbps double matrix, take value in {-1, 1}, the antipodal matrix of all possible bit vectors corresponding to 1 constellation points, logic 1 is mapped to -1, and logic 0 is mapped to 1
 * @param prio_LLR_vec prior 1-by-(Nbps * n_symbol) double matrix, the a priori knowledge on each inner coded bit  
 * @param sym_mod_mat Q-by-M complex<double> matrix, each row corresponding to the M symbols across all (re)transmissions corresponding to the same row in bit_mat_anti, 
 * @param noise power double scalar, the noise power for each (re)transmission
 * @param M int scalar, number of transmissions
 * @param Nbps int scalar number of bits per symbol
 * @param n_symbol number of symbols in one LDPC frame
 * @see 
 *  [1] B. Hochwald, "Achieving near-capacity on a multiple-antenna channel", IEEE Transactions on communications, Mar, 2003
 *  [2] John Thompson,"Extending a fixed-complexity sphere decoder to obtain likelihood information for turbo-mimo systems", IEEE Transactions on vehicular technology, Sep,2008
 *  [3] C.Xiao, and Y.R. Zheng, "on the mutual information and power allocation for vector gaussian channels with finite discrete inputs"
 *  [4] matlab function written by Mingxi Wang. LextDemodulation = MAP_demodulate(y, LextC, chnl_eq, 2*variance, bit_mat_anti, sym_mod_mat, Nr, Ns, Ntime, M)
 */
void MAP_demod_c(double *LextDemodulation, complex <double> *rx_signal, complex <double> *chnl_eq, double *bit_mat_anti, double *prio_LLR_vec, complex <double> *sym_mod_mat, double noise_power, int M, int Nbps, int n_symbol)
{
	double Le_p1_tmp1, Le_p1_tmp2, Log_sum_p1; 
	double Le_n1_tmp1, Le_n1_tmp2, Log_sum_n1; 
	double double_tmp1, double_tmp2, real_tmp, imag_tmp; 
	complex <double> cmplx_tmp, cmplx_tmp1;
    complex <double> *chnl_sym_mod;

    int Q = int(pow(2.0, Nbps));
    chnl_sym_mod = new complex <double> [M * Q]; // Channel multiplied by each of the Q symbols, M-by-Q matrix

	//calculate extrinsic LLR of each bit
	for (int i_symbol = 0; i_symbol < n_symbol; i_symbol++) // the "i_symbol" th symbol in rx_signal
	{
        //calculate h * s, chnl_eq * sym_mod_mat 
        for(int m = 0; m < M; m++)
        {
            for(int q = 0; q < Q; q++)
            {
                chnl_sym_mod[m * Q + q] = chnl_eq[m * n_symbol + i_symbol] * sym_mod_mat[q * M + m];
            }
        }
	
		for (int i_bit = 0; i_bit < Nbps; i_bit++) // the "i_bit" th bit corresponding to this symbol
		{
			Log_sum_p1 = double(-INF);
			Log_sum_n1 = double(-INF);
			
			for(int q = 0; q < Q; q++) //q-th row in bit_mat_anti
			{ 
				//calculate P(bk=+1|y),when bit value = +1
				if (int(bit_mat_anti[q * Nbps + i_bit]) == 1) { // find k th bit bit_mat_anti(q, i_bit) = +1 in bit_mat_anti 
					//calculate 0.5*x[k]*La[k]
					Le_p1_tmp1 = vec_sum_mul(bit_mat_anti + Nbps * q, prio_LLR_vec + i_symbol * Nbps, Nbps); 
					Le_p1_tmp1 = Le_p1_tmp1 - bit_mat_anti[q * Nbps + i_bit] * prio_LLR_vec[i_symbol * Nbps + i_bit];	
					Le_p1_tmp1 = Le_p1_tmp1 / 2;
					////calculate -||y- h.*s||^2/sigma^2
					Le_p1_tmp2 = 0;
					for(int m = 0; m < M; m++)
					{
						cmplx_tmp = rx_signal[m * n_symbol + i_symbol] - chnl_sym_mod[m * Q + q]; //y(m)- h(m) * s(m) 
						real_tmp = real(cmplx_tmp);
						imag_tmp = imag(cmplx_tmp);
						double_tmp1 = real_tmp * real_tmp + imag_tmp * imag_tmp ; // |y(m)- h(m) * s(m)|^2
						Le_p1_tmp2 = Le_p1_tmp2 + double_tmp1; 
					}
					Le_p1_tmp2 = -Le_p1_tmp2 / noise_power; 
					Le_p1_tmp1 = Le_p1_tmp1 + Le_p1_tmp2;
					Log_sum_p1 = jacln(Log_sum_p1, Le_p1_tmp1);  // jacobian logarithm, log(e^Log_sum_p1, e^Le_p1_tmp1);
				}
				else if (int(bit_mat_anti[q * Nbps + i_bit]) == -1) { // find k th bit bit_mat_anti(q, i_bit) = -1 in bit_mat_anti 
					////calculate 0.5*x[k]*La[k]
					Le_n1_tmp1 = vec_sum_mul(bit_mat_anti + Nbps * q, prio_LLR_vec + i_symbol * Nbps, Nbps); 
					Le_n1_tmp1 = Le_n1_tmp1 - bit_mat_anti[q * Nbps + i_bit] * prio_LLR_vec[i_symbol * Nbps + i_bit];	
					Le_n1_tmp1 = Le_n1_tmp1 / 2;
					////calculate -||y- h.*s||^2/sigma^2
					Le_n1_tmp2 = 0;
					for (int m = 0; m < M; m++)
					{
						cmplx_tmp = rx_signal[m * n_symbol + i_symbol] - chnl_sym_mod[m * Q + q]; //y(m)- h(m) * s(m)
						real_tmp = real(cmplx_tmp);
						imag_tmp = imag(cmplx_tmp);
						double_tmp2 = real_tmp * real_tmp + imag_tmp * imag_tmp ; // |y(m)- h(m) * s(m)|^2
						Le_n1_tmp2 = Le_n1_tmp2 + double_tmp2; 
					}
					Le_n1_tmp2 = -Le_n1_tmp2 / noise_power; 
					Le_n1_tmp1 = Le_n1_tmp1 + Le_n1_tmp2;
					Log_sum_n1 = jacln(Log_sum_n1, Le_n1_tmp1);  // jacobian logarithm, log(e^Log_sum_p1, e^Le_p1_tmp1);
				}	
			}

			LextDemodulation[i_symbol * Nbps + i_bit] = Log_sum_p1 - Log_sum_n1; //  			
		}	
	}
	delete [] chnl_sym_mod;
	return;
}


/**
 * LextDemodulation = MAP_demod(rx_signal, chnl_eq, bit_mat_anti, LextC, sym_mod_mat, noise_power)
 * The cmex gateway routine to the function of Chase combining HARQ MAP detector 
 * @param rx_signal M-by-n_symbol complex matrix, the received signal corresponding to 1 LDPC frame and M (re)transmisstions
 * @param chnl_eq M-by-n_symbol complex matrix, equivalent channel when the noises or the M (re)transmissions are normalized to have equal power
 * @param bit_mat_anti Q-by-Nbps real matrix, take value in {-1, 1}, the antipodal matrix of all possible bit vectors corresponding to 1 constellation points, logic 1 is mapped to -1, and logic 0 is mapped to 1
 * @param prio_LLR_vec prior 1-by-(Nbps * n_symbol) real matrix, the a priori knowledge on each inner coded bit  
 * @param sym_mod_mat Q-by-M complex matrix, each row corresponding to the M symbols across all (re)transmissions corresponding to the same row in bit_mat_anti, 
 * @param noise power real scalar, the noise power for each (re)transmission
 * @return 1-by-(Nbps * n_symbol) matrix, the extrinsic LLR 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// check input and output
	// check number of inputs and outputs
	if (nrhs != 6){  
		mexErrMsgTxt("number of inputs errors)");
	}
	if (nlhs > 1){
		mexErrMsgTxt("number of outputs errors");
	}
	//check inputs and outputs are real or complex
	if ((!mxIsComplex(prhs[0])) || (!mxIsComplex(prhs[1])) || (!mxIsComplex(prhs[4]))){
		mexErrMsgTxt("y, chnl_eq, and sym_mod_mat must be complex.");
	}
	else if((mxIsComplex(prhs[2])) || (mxIsComplex(prhs[3])) || (mxIsComplex(prhs[5]))) {
		mexErrMsgTxt("bit_mat_anti, LextC and variance must be real.");
	}

	// Retrieve the dimensions of the index data and check for consistency
	int M = mxGetM(prhs[0]);
	int n_symbol = mxGetN(prhs[0]);
	int Nbps = mxGetN(prhs[2]);
	int Q = int(pow(2.0, Nbps));
	if (mxGetM(prhs[1]) != M || mxGetN(prhs[1]) != n_symbol ||
		mxGetM(prhs[2]) != Q || mxGetN(prhs[2]) != Nbps ||
		mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != Nbps * n_symbol ||
		mxGetM(prhs[4]) != Q || mxGetN(prhs[4]) != M ||
		mxGetM(prhs[5]) != 1 || mxGetN(prhs[5]) != 1)
	{
		mexErrMsgTxt("Input arguments sizes are not consistent!");
	}

    //Retrieve the input data 
	//get rx_signal and convert to matrix order in C language(cross first, then down)
	double *rx_signal_real, *rx_signal_imag;
	complex <double> *rx_signal_cmplx; 

	rx_signal_real = mxGetPr(prhs[0]);     
	rx_signal_imag = mxGetPi(prhs[0]);
	rx_signal_cmplx = new complex <double> [M * n_symbol];
	for(int i_symbol = 0; i_symbol < n_symbol; i_symbol++)
	{
		for(int m = 0; m < M; m++)
		{
			rx_signal_cmplx[m * n_symbol + i_symbol] = complex <double>(rx_signal_real[i_symbol * M + m], rx_signal_imag[i_symbol * M + m]); 
		}	
	}

	//get chnl_eq and convert to C matrix(cross first, then down)
	double *chnl_eq_real, *chnl_eq_imag;
	complex <double> *chnl_eq_cmplx; 
	
	chnl_eq_real = mxGetPr(prhs[1]);     
	chnl_eq_imag = mxGetPi(prhs[1]);

	chnl_eq_cmplx = new complex <double> [M * n_symbol];
	for(int m = 0; m < M; m++) // convert to complex class with order in C
	{
        for(int i_symbol = 0; i_symbol < n_symbol; i_symbol++)
        {
            chnl_eq_cmplx[m * n_symbol + i_symbol] = complex <double>(chnl_eq_real[i_symbol * M + m], chnl_eq_imag[i_symbol * M + m]);
        }
	}	

    //get bit_mat_anti
	double *bit_mat_anti, *bit_mat_anti_c;
	bit_mat_anti = mxGetPr(prhs[2]);     

	bit_mat_anti_c = new double [Q * Nbps];
    for(int i_bit = 0; i_bit < Nbps; i_bit++)   
	{
		for(int q = 0; q < Q; q++)
		{
			bit_mat_anti_c[q * Nbps + i_bit] = bit_mat_anti[i_bit * Q + q];
		}	
	}

	//get prio_LLR_mat and convert to a vector form
	double *prio_LLR_mat, *prio_LLR_vec;

	prio_LLR_mat = mxGetPr(prhs[3]);

	prio_LLR_vec = new double [n_symbol * Nbps];
	for(int i_symbol = 0; i_symbol < n_symbol; i_symbol++)
	{
		for (int i_bit = 0; i_bit < Nbps; i_bit++)
		{
			prio_LLR_vec[i_symbol * Nbps + i_bit] = prio_LLR_mat[i_symbol * Nbps + i_bit];
		}
	}

	//get sym_mod_mat, and convert to C matrix(cross first, then down)
	double *sym_mod_mat_real, *sym_mod_mat_imag;
	complex <double> *sym_mod_mat_cmplx;  
	
	sym_mod_mat_real = mxGetPr(prhs[4]);     
	sym_mod_mat_imag = mxGetPi(prhs[4]);	

	sym_mod_mat_cmplx = new complex <double> [Q * M];
    for(int m = 0; m < M; m++)   // convert to complex class
	{
		for(int q = 0; q < Q; q++)
		{
			sym_mod_mat_cmplx[q * M + m]= complex <double>(sym_mod_mat_real[m * Q + q], sym_mod_mat_imag[m * Q + q]); 
		}	
	}

    double noise_power = *mxGetPr(prhs[5]);

    plhs[0] = mxCreateDoubleMatrix(1, n_symbol * Nbps, mxREAL); //Create an mxArray for the output data 
	double *LextDemodulation = mxGetPr(plhs[0]); // Create a pointer to the output data
	
	MAP_demod_c(LextDemodulation, rx_signal_cmplx, chnl_eq_cmplx, bit_mat_anti_c, prio_LLR_vec, sym_mod_mat_cmplx, noise_power, M, Nbps, n_symbol);

    delete [] rx_signal_cmplx;
    delete [] chnl_eq_cmplx;
	delete [] bit_mat_anti_c;
	delete [] prio_LLR_vec;
	delete [] sym_mod_mat_cmplx;
}
