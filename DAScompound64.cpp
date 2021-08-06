/*     include librarys     */
//#include <math.h>
#include <stdlib.h>
#include "mex.h"
// #include "ap_int.h"
// #include "ap_fixed.h"

#define Nt 2560
#define N  253
#define M  390
#define Nf 21
#define cs 1540
#define fs 31250000
#define H  390
#define W  23

//	const float alpha[21] = { -0.17453292,
//							  -0.15707964,
//							  -0.13962634,
//							  -0.12217305,
//							  -0.10471976,
//							  -0.087266460,
//							  -0.069813170,
//							  -0.052359879,
//							  -0.034906585,
//							  -0.017453292,
//							   0,
//							   0.017453292,
//							   0.034906585,
//							   0.052359879,
//							   0.069813170,
//							   0.087266460,
//							   0.10471976,
//							   0.12217305,
//							   0.13962634,
//							   0.15707964,
//							   0.17453292 };
//
//	const float sa[21]    = { -0.17364818,
//							  -0.15643446,
//							  -0.1391731,
//							  -0.1218693,
//							  -0.1045285,
//							  -0.08715574,
//							  -0.06975647,
//							  -0.05233596,
//							  -0.03489950,
//							  -0.01745241,
//							   0,
//							   0.01745241,
//							   0.03489950,
//							   0.05233596,
//							   0.06975647,
//							   0.08715574,
//							   0.10452846,
//							   0.12186934,
//							   0.13917311,
//							   0.15643446,
//							   0.17364818 };
//
//	const float ca[21]    = {  0.98480773,
//							   0.98768836,
//							   0.99026805,
//							   0.99254614,
//							   0.99452192,
//							   0.99619472,
//							   0.99756408,
//							   0.99862951,
//							   0.99939084,
//							   0.99984771,
//							   1,
//							   0.99984771,
//							   0.99939084,
//							   0.99862951,
//							   0.99756408,
//							   0.99619472,
//							   0.99452192,
//							   0.99254614,
//							   0.99026805,
//							   0.98768836,
//							   0.98480773 };

short sqrt(int in) {

	int op, one, res;
	// int res;
	int i, j;

	op = in;
	res = 0;
	one = 4194304;

	sqrt1: for (i = 0; i < 12; ++i) {
		if (one > op)
			one = one >> 2;
		else continue;
	}

	sqrt2: for (j = 0; j < 12; ++j) {
		if (one != 0) {
			if (op >= res + one) {
				op = op - (res + one);
				res = res + 2*one;
				res = res >> 1;
				one = one >> 2;
			}
			else {
				res = res >> 1;
				one = one >> 2;
			}
		}
		else continue;
	}
	return (short)res;
}

void computeDelay(short xele, short *x, short *y, int R, short *out)
{
	// static short dist;
	int sq;
	int i;

	delay_col: for (i = 0; i < N; ++i) {
// #pragma HLS pipeline II=1
		sq = (x[i] - xele)*(x[i] - xele) + y[R] * y[R];
		out[i] = sqrt(sq);
	}
}

void computeAcc(short *pt, short *y, short *Delay, int R, int *acc_out)
{

	short delay;
	int i;

	rf_row: for (i = 0; i < N; ++i) {
// #pragma HLS pipeline II=2
		delay = (Delay[i] + y[R]);
		if (delay >= 0 && delay <= Nt-1) acc_out[i] = pt[delay];
		else acc_out[i] = 0;
	}
}

void DAScompound64(short *pt, short *xele, short *x, short *y, int *rf)
{
// #pragma HLS interface s_axilite port=pt bundle=DIn
// #pragma HLS interface s_axilite port=xele bundle=Xe
// #pragma HLS interface s_axilite port=x bundle=X
// #pragma HLS interface s_axilite port=y bundle=Y
// #pragma HLS interface axis port=rf name=M_AXIS_DOut
// #pragma HLS interface s_axilite register port=return
	
	short *Delay = (short *)malloc(sizeof(short)*N);
	int *acc = (int *)malloc(sizeof(int)*N);

	int k, R;

	row_fill: for (R = 0; R < M; ++R) {
// #pragma HLS loop_flatten off

		computeDelay(xele[0], x, y, R, Delay);

		computeAcc(pt, y, Delay, R, acc);

		rf_out_col: for (k = 0; k < N; ++k) {
			rf[R + k*M]  = acc[k];
		}
	}

	free(Delay);
	free(acc);
}

/* [rf] = DAScompound64(short pt, float xele, float x, float y, char Num) */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])
{
/* define the hyper-parameters */
// int N = 253;
// int M = 390;

/* input setting */
short *pt;
int   *rf;
short *x, *y, *xele;

/* check for proper number of arguments */
if (nrhs != 4) {
	mexErrMsgIdAndTxt("MATLAB:DASBF:invalidNumInputs",
		"Four input arguments required.");
}
else if (nlhs > 1) {
	mexErrMsgIdAndTxt("MATLAB:DASBF:maxlhs",
		"Too many output arguments.");
}

/* create matrices for the return argument */
plhs[0] = mxCreateNumericMatrix(M, N, mxINT32_CLASS, mxREAL);

/* assign pointrs to the input and output parameters */
pt   = (short *)mxGetPr(prhs[0]);
xele = (short *)mxGetPr(prhs[1]);
x    = (short *)mxGetPr(prhs[2]);
y    = (short *)mxGetPr(prhs[3]);
rf   = (int *)mxGetPr(plhs[0]);

/* Do the actual computations in a subroutine */
DAScompound64(pt, xele, x, y, rf);

return;
}
