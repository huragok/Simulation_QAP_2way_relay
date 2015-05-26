//file: Map_demod.cpp
//description: function of MIMO MAP detector, and its mex function 
//author: Weiliang Zeng
//created date: 03/01/10
//modified: 
//on 03/09/10: add accurate and piece-wise approaximate method to calculate log(e^x+e^y) 

#include "mex.h"
#include <complex>
#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>

using namespace std;

#define INF  10000  // infinity

// description: return summation of multiplication of two vecters. sum(a(i)*b(i)), i = 1,2...L.  L is the length of vecter a and b
inline complex <double> vec_sum_mul(complex<double> *vec_a,  complex<double> *vec_b, int L)
{
    int i; 
	complex <double> tmp = complex <double> (0,0);
	for(i = 0; i < L; i++)
	{
          tmp = tmp + vec_a[i]*vec_b[i]; 
	}
	return tmp;
}

// description: return summation of multiplication of two vecters. sum(a(i)*b(i)), i = 1,2...L.  L is the length of vecter a and b
inline double vec_sum_mul(double*vec_a,  double*vec_b, int L)
{
    int i; 
	double tmp = 0;
	for(i = 0; i < L; i++)
	{
          tmp = tmp + vec_a[i]*vec_b[i]; 
	}
	return tmp;
}

//description: return jacbian logarithm. log(e^x + e^y) = max(x,y) + log(1+e^(-|x-y|)).
//correction function r(x) = log(1+e^-x), use piease wise linear function to approximate r(x) as table II in [1]
//reference:
//[1] Xiao-Yu Hu,Evangelos Eleftheriou,Dieter-Michael Arnold,and Ajay Dholakia. "Efficient Implementations of the Sum-Product Algorithm for Decoding LDPC Codes". GlobalCom , 2001:ppV2 1036-1036E
inline double jacln(double x, double y)
{
	// find larger one of x and y
	double max_input, min_input;
	if (x >= y){  
	     max_input = x;
         min_input = y;
	}
	else{
	     max_input = y;
         min_input = x;
	}

	// find correction r(max_input - min_input)
   double corr_input, corr_output, result;
   corr_input = max_input - min_input;
	if ((corr_input >= 0) &&(corr_input < 0.5)) {
		corr_output = -corr_input*0.5+0.7;
	}
	else if ((corr_input >= 0.5) &&(corr_input < 1.6)){
		corr_output = -corr_input*0.25+0.575;
	}
	else if ((corr_input >= 1.6) &&(corr_input < 2.2)){
		corr_output = -corr_input*0.125+0.375;
	}
	else if ((corr_input >= 2.2) &&(corr_input < 3.2)){
		corr_output = -corr_input*0.0625+0.2375;
	}
	else if ((corr_input >= 3.2) &&(corr_input < 4.4)){
		corr_output = -corr_input*0.03125+0.1375;
	}
	else {
		corr_output = 0;
	}	

	// log(e^x + e^y) = max(x,y) + log(1+e^(-|x-y|)).
	result = max_input + corr_output;
    return (result);

}


// description: function of MIMO MAP detector 
//LextDemodulation: real extrinsinc LLR in matrix, Ns x Mc*K*L, down first then cross as in matlab mex function
//rx_signal: received signal y = Hs+n, Nr*K x L, L is number of blocks
//chnl_eq: equivalent channel matrix, He = H*G: Nr*K x Ns*K, 
//bit_mat_anti:antipodal matrix of bit vectors,2^Ns*Mc*K x Ns*Mc*K, logic 1 is mapped to -1, and logic 0 is mapped to 1
//prio_LLR_vec: prio information La(b) of multiple bit streams in matrix, Ns x Mc*K*L,   
//sym_mod_mat: matrix of modulated symbols, corresponding to bit_mat_anti, 2^Ns*Mc*K x Ns*k 
//noise power: total noise power
//Nr: the number of receive antenna
//Ns: the number of bit streams
//K: number of temporal vectors in linear precoding 
//Mc: bit per symbol of modulation
//block_num: number of blocks of received signal y, L
//reference:
//[1]B. Hochwald, "Achieving near-capacity on a multiple-antenna channel",
//IEEE Transactions on communications, Mar, 2003
//[2]John Thompson,"Extending a fixed-complexity sphere decoder to obtain 
//likelihood information for turbo-mimo systems", IEEE Transactions on vehicular technology, Sep,2008
//[3]C.Xiao, and Y.R. Zheng, "on the mutual information and power allocation for vector gaussian channels with finite discrete inputs"
//[4] matlab function written by Mingxi Wang. LextDemodulation = MAP_demodulate(y, LextC, chnl_eq, 2*variance, bit_mat_anti, sym_mod_mat, Nr, Ns, Ntime, M); 
                
void MAP_demod_c(double *LextDemodulation, complex <double> *rx_signal, complex <double> *chnl_eq, double *bit_mat_anti,double *prio_LLR_vec,  complex <double> *sym_mod_mat, double noise_power,   int Nr, int Ns, int K, int Mc, int block_num)
{
	int i,j,k, row_index, col_index;
	int block_index, time_index, bit_index_c, bit_index_mex;
	int Ns_K_Mc, bit_mat_row_num;
	double Le_p1_tmp1, Le_p1_tmp2, Log_sum_p1, sum_exp_p1; 
	double Le_n1_tmp1, Le_n1_tmp2, Log_sum_n1, sum_exp_n1; 
	double double_tmp1, real_tmp, imag_tmp; 
	complex <double> cmplx_tmp, cmplx_tmp1;
    complex <double> *chnl_sym_mod;

	Ns_K_Mc = Ns*K*Mc;
    bit_mat_row_num = int(pow(2.0,Ns_K_Mc));
    chnl_sym_mod = new complex <double> [Nr*K*bit_mat_row_num];

    //calculate H*s, chnl_eq* sym_mod_mat 
	for(i = 0; i< Nr*K; i++)
	{
          for(j = 0; j < bit_mat_row_num; j++)
		  {
                 chnl_sym_mod[i*bit_mat_row_num+j] = vec_sum_mul(chnl_eq+Ns*K*i, sym_mod_mat+Ns*K*j, Ns*K);;
		  }
	}
	
	//calculate extrinsic LLR of each bit
	for (block_index = 0; block_index < block_num; block_index++) // the "block_index" th block in rx_signal
	{
	    for (time_index = 0; time_index < K; time_index++) // the "time_index"th vecter in each block
		{
            for (i = 0; i < Ns*Mc; i++) 
			{
				Log_sum_p1 = double(-INF);
				Log_sum_n1 = double(-INF);
				/*sum_exp_p1 = 0; 
				sum_exp_n1 = 0; */
				row_index = int (i % Ns) ;
				col_index = int ((i - row_index) /Ns);
                bit_index_c = time_index*Ns*Mc + row_index*Mc + col_index; 
				
				for(j= 0; j< bit_mat_row_num; j++) //j th row in bit_mat_anti
				{ 
					//calculate P(bk=+1|y),when bit value = +1
					if (int(bit_mat_anti[j*Ns_K_Mc+bit_index_c]) == 1) { // find k th bit bit_mat_anti[j][bit_index_c] = +1 in bit_mat_anti 
					    //calculate 0.5*b[k]*La[k]
						Le_p1_tmp1 = vec_sum_mul(bit_mat_anti+Ns_K_Mc*j, prio_LLR_vec +block_index*Ns_K_Mc, Ns_K_Mc); 
					    Le_p1_tmp1 = Le_p1_tmp1 - bit_mat_anti[j*Ns_K_Mc+bit_index_c]*prio_LLR_vec[block_index*Ns_K_Mc+bit_index_c];	
						Le_p1_tmp1 = Le_p1_tmp1 /2;
						////calculate -||y- H*s||^2/sigma^2
                        Le_p1_tmp2 = 0;
						for(k = 0; k< Nr*K; k++)
						{
						    cmplx_tmp = rx_signal[k*block_num + block_index] - chnl_sym_mod[k*bit_mat_row_num+j]; //y(k)- H*s 
                            real_tmp = real(cmplx_tmp);
							imag_tmp = imag(cmplx_tmp);
							double_tmp1 = real_tmp*real_tmp + imag_tmp *imag_tmp ; //|y(k)- H*s|^2
							Le_p1_tmp2 = Le_p1_tmp2 + double_tmp1; 
						}
						Le_p1_tmp2 = -Le_p1_tmp2/noise_power; 
                        Le_p1_tmp1 = Le_p1_tmp1+ Le_p1_tmp2;
						//sum_exp_p1 = sum_exp_p1 + exp(Le_p1_tmp1); // sum(exp(-|| y - H*s||^2 + 0.5*b*La))
						Log_sum_p1 = jacln(Log_sum_p1, Le_p1_tmp1);  // jacobian logarithm, log(e^Log_sum_p1, e^Le_p1_tmp1);
					}
					else if (int(bit_mat_anti[j*Ns_K_Mc+bit_index_c]) == -1) { //find k th bit bit_mat_anti[j][bit_index_c] = -1 in bit_mat_anti 
						////calculate 0.5*b[k]*La[k]
						Le_n1_tmp1 = vec_sum_mul(bit_mat_anti+Ns_K_Mc*j, prio_LLR_vec +block_index*Ns_K_Mc, Ns_K_Mc); 
					    Le_n1_tmp1 = Le_n1_tmp1 - bit_mat_anti[j*Ns_K_Mc+bit_index_c]*prio_LLR_vec[block_index*Ns_K_Mc+bit_index_c];	
						Le_n1_tmp1 = Le_n1_tmp1 /2;
						//////calculate -||y- H*s||^2/sigma^2
                        Le_n1_tmp2 = 0;
						for(k = 0; k< Nr*K; k++)
						{
						    cmplx_tmp = rx_signal[k*block_num + block_index] - chnl_sym_mod[k*bit_mat_row_num+j]; //y(k)- H*s 
							real_tmp = real(cmplx_tmp);
							imag_tmp = imag(cmplx_tmp);
							double_tmp1 = real_tmp*real_tmp + imag_tmp *imag_tmp ; //|y(k)- H*s|^2
							Le_n1_tmp2 = Le_n1_tmp2 + double_tmp1; 
						}
						Le_n1_tmp2 = -Le_n1_tmp2/noise_power; 
						Le_n1_tmp1 = Le_n1_tmp1+ Le_n1_tmp2;
						//sum_exp_n1 = sum_exp_n1 + exp(Le_n1_tmp1);  // sum(exp(-|| y - H*s||^2 + 0.5*b*La))
						Log_sum_n1 = jacln(Log_sum_n1, Le_n1_tmp1);  // jacobian logarithm, log(e^Log_sum_p1, e^Le_p1_tmp1);
					}	
	            }
				//Log_sum_p1 = log(sum_exp_p1); // sum(exp(-|| y - H*s||^2 + 0.5*b*La)) for bit = +1
	   //         Log_sum_n1 = log(sum_exp_n1); // sum(exp(-|| y - H*s||^2 + 0.5*b*La)) for bit = -1

				bit_index_mex = block_index*Ns_K_Mc + time_index*Ns*Mc+ col_index*Ns+row_index; // down first, then cross, 2D array as stored in matlab mex 		
				LextDemodulation[bit_index_mex] = Log_sum_p1 - Log_sum_n1 ; //  			
			}	
		}
	}

}


//description: mex function of void MAP_demod_c(double *LextDemodulation, complex <double> *rx_signal, complex <double> *chnl_eq, double *bit_mat_anti,double *prio_LLR_vec,  complex <double> *sym_mod_mat, double noise_power,   int Nr, int Ns, int K, int Mc, int block_num)
// interface between matlab and C program
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, j,k;
	int row_num,col_num;
	double *LextDemodulation;

	// check input and output
	// check number of inputs and outputs
	if (nrhs  != 6){  
		mexErrMsgTxt("number of inputs errors)");
	}
	if (nlhs  != 1){
		mexErrMsgTxt("number of outputs errors");
	}
	////check inputs and outputs are real or complex
	//if ( (!mxIsComplex(prhs[0]))||(!mxIsComplex(prhs[1]))||(!mxIsComplex(prhs[4])) ){
	//	mexErrMsgTxt("y, chnl_eq, and sym_mod_mat must be complex");
	//}
	//else if( (mxIsComplex(prhs[0]))||(mxIsComplex(prhs[1]))||(mxIsComplex(prhs[4]))) {
	//	mexErrMsgTxt("bit_mat_anti, LextC and 2*variance must be complex");
	//}

    //Retrieve the input data 
	//get rx_signal and convert to matrix order in C language(cross first, then down)
	int Nr_K, block_num;
	double *rx_signal_real,*rx_signal_imag;
	complex <double> *rx_signal_cmplx; 

	rx_signal_real = mxGetPr(prhs[0]);     
	rx_signal_imag = mxGetPi(prhs[0]);
    Nr_K = int(mxGetM(prhs[0])); //get number of rows
	block_num  = int(mxGetN(prhs[0]));
	rx_signal_cmplx = new complex <double> [Nr_K*block_num];
	for(i= 0; i< block_num; i++)
	{
		for(k= 0; k< Nr_K; k++)
		{
			rx_signal_cmplx[k*block_num+i]= complex <double>(rx_signal_real[i*Nr_K+ k], rx_signal_imag[i*Nr_K+k]); 
		}	
	}

	//get chnl_eq and convert to C matrix(cross first, then down)
	int Ns_K;
	double *chnl_eq_real, *chnl_eq_imag;
	complex <double> *chnl_eq_cmplx; 
	
	chnl_eq_real = mxGetPr(prhs[1]);     
	chnl_eq_imag = mxGetPi(prhs[1]);
	Ns_K = int(mxGetN(prhs[1]));

	chnl_eq_cmplx = new complex <double> [Nr_K*Ns_K];
    for(i= 0; i< Ns_K; i++)   // convert to complex class with order in C
	{
		for(k= 0; k< Nr_K; k++)
		{
			chnl_eq_cmplx[k*Ns_K+i]= complex <double>(chnl_eq_real[i*Nr_K+ k], chnl_eq_imag[i*Nr_K+k]); 
		}	
	}

    //get bit_mat_anti
	int Ns_K_Mc, Mc;
	double *bit_mat_anti, *bit_mat_anti_c;
	
	bit_mat_anti = mxGetPr(prhs[2]);     
	Ns_K_Mc = int(mxGetN(prhs[2]));
	Mc = int(Ns_K_Mc/Ns_K);
	row_num = int(pow(2.0,Ns_K_Mc));

	bit_mat_anti_c = new double  [row_num*Ns_K_Mc];
    for(i= 0; i< Ns_K_Mc; i++)   
	{
		for(k= 0; k< row_num; k++)
		{
			bit_mat_anti_c[k*Ns_K_Mc+i] = bit_mat_anti[i*row_num+k]; //double(bit_mat_anti[i*row_num+k]);  
		}	
	}


	//get prio_LLR_mat and convert to a vector form
    int Ns, K, Nr; 
	double *prio_LLR_mat, *prio_LLR_vec;

	prio_LLR_mat = mxGetPr(prhs[3]);     
    Ns = int(mxGetM(prhs[3])); 
	col_num  = int(mxGetN(prhs[3]));
    K = int(Ns_K/Ns);
    Nr = int(Nr_K/K);

	prio_LLR_vec = new double [Ns*block_num*K*Mc];
	for(i= 0; i< block_num*K; i++)
	{
		for(k= 0; k< Ns; k++)
		{
			for (j = 0; j< Mc; j++)
			{
			      prio_LLR_vec[i*Ns*Mc+k*Mc+j] = prio_LLR_mat[i*Ns*Mc+j*Ns+k];
			}
		}	
	}

	//get sym_mod_mat, and convert to C matrix(cross first, then down)
	double *sym_mod_mat_real, *sym_mod_mat_imag;
	complex <double> *sym_mod_mat_cmplx;  
	
	sym_mod_mat_real = mxGetPr(prhs[4]);     
	sym_mod_mat_imag = mxGetPi(prhs[4]);	
	row_num = int(mxGetM(prhs[4]));
	col_num = int(mxGetN(prhs[4]));

	sym_mod_mat_cmplx = new complex <double> [row_num*col_num];
    for(i= 0; i< col_num; i++)   // convert to complex class
	{
		for(k= 0; k< row_num; k++)
		{
			sym_mod_mat_cmplx[k*col_num+i]= complex <double>(sym_mod_mat_real[i*row_num+ k], sym_mod_mat_imag[i*row_num+k]); 
		}	
	}

    double *noise_power; // get noise_power
	noise_power = mxGetPr(prhs[5]);

    plhs[0] = mxCreateDoubleMatrix(Ns,block_num*K*Mc,  mxREAL); //Create an mxArray for the output data 
	//plhs[0] = mxCreateDoubleMatrix(1,3,  mxREAL); //for test
	LextDemodulation = mxGetPr(plhs[0]); // Create a pointer to the output data
	
	MAP_demod_c(LextDemodulation, rx_signal_cmplx, chnl_eq_cmplx, bit_mat_anti_c,prio_LLR_vec,  sym_mod_mat_cmplx, *noise_power,  Nr, Ns, K, Mc, block_num);

    delete [] rx_signal_cmplx;
    delete [] chnl_eq_cmplx;
	delete [] bit_mat_anti_c;
	delete [] prio_LLR_vec;
	delete [] sym_mod_mat_cmplx;
}
