/*===================================================
*
* DAScompound32x.cpp is used to generate corresponding .MEX file to
* 		conduct beamforming much quicker in MATLAB.
* This function is for 32 elements (hommade probe) and 32T32R only and the method is compund method.
* 		[rf] = DAScompound32(pt,N_active,fs,cs,xele,x,y,alpha)
* Finished by Yuxiang Ma on 13 Aug, 2019.
* 
*====================================================*/

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
* Notes!!!
* 1. The data type of all inputs should be double (int16, int32 or other type will lead to memory overread or memory disorder). 
* 2. If the program crash, please first check the input form pt. (Wrong format of them will surely cause a memory overread).
* 3. The program has been checked for many times, there is surely no bugs in it. So if you have any trouble with the program, please check the input.
* 4. All mexPrintf lines which are commented out can help you debug the program, so please don't delete them.
* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

/*     include librarys     */
#include <math.h>
#include <stdlib.h>
#include "mex.h"
//#include "mexAdapter.hpp"

/*     define namespace     */
// using namespace matlab::data;
// using matlab::mex::ArgumentList;

/*     values in mxArrays is stored in a column-wise order.    */

void DAScomNtNr128fh(double *pt, double* alpha, double *xele, double *x,
	double *y, unsigned int fs, unsigned int cs, double *rf,
	unsigned int Nf, unsigned int Nx, unsigned int M, unsigned int N, int Nt)
{
	unsigned int pos[128];
//	unsigned int Nc
//	Nc = 64;
	for (unsigned int nc = 0; nc < Nx; nc++){
		pos[nc] = nc; 
}
	/*   compute the delay   */
	
	unsigned int i;
	double valsum1;
	double dist;
	
	double *sa, *ca;
	sa = (double *)malloc(sizeof(double)*Nf);
	ca = (double *)malloc(sizeof(double)*Nf);
	for (unsigned int nf = 0; nf < Nf; nf++) {
		sa[nf] = sin(alpha[nf]);
		ca[nf] = cos(alpha[nf]);
	}
	
	double ***Delay;
	Delay = (double ***)malloc(sizeof(double **)*Nx);
	/*	Delay time of backward propagation	*/
	double dt;
	
	for (unsigned int i = 0; i < Nx; i++) {
		Delay[i] = (double **)malloc(sizeof(double *)*N);
		for (unsigned int j = 0; j < N; j++) {
			Delay[i][j] = (double *)malloc(sizeof(double)*M);
			for (unsigned int k = 0; k < M; k++) {
				dist = sqrt((x[j] - xele[i])*(x[j] - xele[i]) + y[k] * y[k]);
				dt = dist * fs / cs;
				Delay[i][j][k] = (unsigned int)dt;
			}
		}
	}
	
	int delay;
	for (unsigned int nf = 0; nf < Nf; nf++) {
		for (unsigned int nx = 0; nx < N ; nx++) {
			for (unsigned int ny = 0; ny < M; ny++) {
				valsum1 = 0;
				for (unsigned int n = 0; n < Nx; n++) {
					if ( (xele[n] - x[nx])/y[ny] > 1 || (xele[n] - x[nx])/y[ny] < -1  )
						continue;
					if (alpha[nf] >= 0 ){
						delay = (int)(Delay[n][nx][ny] + (y[ny]*ca[nf] + x[nx]*sa[nf])*fs / cs + 0.5);
					}
					else{
						delay = (int)(Delay[n][nx][ny] + (y[ny]*ca[nf] + (x[nx] - xele[Nx-1])*sa[nf])*fs / cs + 0.5);
					}
					if (delay >= 0){
						if (delay <= Nt-1){
							valsum1 = valsum1 + pt[ pos[n]*Nt*Nf + delay + nf*Nt];
						}
					}
				}
				rf[nx*M + ny] = rf[nx*M + ny] + valsum1;
			}
		}
	}
	
	for (unsigned int i = 0; i < Nx; i++) {
		for (unsigned int j = 0; j < N; j++) {
			free(Delay[i][j]);
		}
		free(Delay[i]);
	}
	free(Delay);
	
	return;
}

/* [rf] = DAScompound32(pt,fs,cs,xele,x,y,alpha)     			*/
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])
{
/* define the hyper-parameters */
// unsigned int Nx = 32;
unsigned int Nx;
Nx = (unsigned int)*mxGetPr(prhs[7]);
unsigned int Nf;
double *xele;

/* input setting */
double *pt, *rf;
double *x, *y;
double *alpha;
unsigned int fs, cs, Ntim;

unsigned int rn, cn;
unsigned int M, N;
/* check for proper number of arguments */
if (nrhs != 8) {
	mexErrMsgIdAndTxt("MATLAB:DASBF:invalidNumInputs",
		"Seven input arguments required.");
}
else if (nlhs > 1) {
	mexErrMsgIdAndTxt("MATLAB:DASBF:maxlhs",
		"Too many output arguments.");
}

rn = (unsigned int)mxGetM(prhs[0]);
cn = (unsigned int)mxGetN(prhs[0]);
M  = (unsigned int)mxGetN(prhs[5]);
N  = (unsigned int)mxGetN(prhs[4]);

/* create matrices for the return argument */
plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);

/* assign pointrs to the input and output parameters */
pt = (double *)mxGetPr(prhs[0]);
/* Na = (unsigned int)*mxGetPr(prhs[1]); */
xele = mxGetPr(prhs[3]);
//Ntim = (unsigned int)*mxGetPr(prhs[2]);
x = (double *)mxGetPr(prhs[4]);
y = (double *)mxGetPr(prhs[5]);
fs = (unsigned int)*mxGetPr(prhs[1]);
cs = (unsigned int)*mxGetPr(prhs[2]);
alpha = (double *)mxGetPr(prhs[6]);
// mexPrintf("%d\n fs",fs);
// mexPrintf("%d\n cs",cs);

Nf = (unsigned int)mxGetN(prhs[6]);
Ntim = rn/Nf;
// mexPrintf("%d\n Nf",Nf);

/*
if (rn != Nf*Ntim){
	mexPrintf("row number is %d, which is not consistent to Ntim and Nfra", rn);
	return;
}
*/


rf = mxGetPr(plhs[0]);

/* Do the actual computations in a subroutine */
DAScomNtNr128fh(pt,alpha, xele, x, y, fs, cs, rf, Nf, Nx, M, N, (int)Ntim);

//mexPrintf("20e%dt%drf Beamforming is completed !\n",Na, Na);

return;
}
