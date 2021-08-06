#include <stdlib.h>
#include "ap_axi_sdata.h"
#include "hls_stream.h"

#define Nt 2560
#define N  253
#define M  390
#define Nf 21
#define cs 1540
#define fs 31250000
#define H  390
#define W  23

ap_int<12> sqrt(ap_int<24> in) {

	static ap_uint<24> op, one, res;
	static int i, j;

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
	return (ap_int<12>)res;
}

void computeDelay(hls::stream<ap_axis<16,1,1,1>> *xele, hls::stream<ap_axis<16,1,1,1>> *x, hls::stream<ap_axis<16,1,1,1>> *y, int R, ap_int<12> out[N])
{
	static ap_int<24> sq;
	static ap_int<16> xr, xer, yr;
	static int i;

	xer = xele[0].read().data;
	yr = y[R].read().data;

	delay_col: for (i = 0; i < N; ++i) {
#pragma HLS pipeline II=1
		xr = x[i].read().data;
		sq = (xr - xer)*(xr - xer) + yr*yr;
		out[i] = sqrt(sq);
	}
}

void computeAcc(hls::stream<ap_axis<16,1,1,1>> *pt, hls::stream<ap_axis<16,1,1,1>> *y, ap_int<12> Delay[N], int R, hls::stream<ap_axis<16,1,1,1>> *acc_out)
{

	static ap_int<12> delay;
	static ap_int<16> yr;
	static ap_axis<16,1,1,1> temp;
	static int i;

	yr = y[R].read().data;

	temp.data = 0;
	temp.keep = 1;
	temp.strb = 1;
	temp.user = 1;
	temp.last = 0;
	temp.id = 0;
	temp.dest = 1;

	rf_row: for (i = 0; i < N; ++i) {
#pragma HLS pipeline II=1
		delay = (Delay[i] + yr);
		if (delay >= 0 && delay <= Nt-1) acc_out[R*N + i].write(pt[delay].read());
		else acc_out[R*N + i].write(temp);
	}
}

void DAScompound64(hls::stream<ap_axis<16,1,1,1>> *DIn, hls::stream<ap_axis<16,1,1,1>> *Xele, hls::stream<ap_axis<16,1,1,1>> *X, hls::stream<ap_axis<16,1,1,1>> *Y, hls::stream<ap_axis<16,1,1,1>> *DOut)
{
#pragma HLS interface axis port=DIn
#pragma HLS interface axis port=Xele
#pragma HLS interface axis port=X
#pragma HLS interface axis port=Y
#pragma HLS interface axis port=DOut
#pragma HLS interface s_axilite register port=return
	
	ap_int<12> Delay[N];
	int k, R;

	/*   Delay time of backward propagation   */

//	config_load: {
//		load_loop: for (k = 0; k < Nt; ++k) {
//#pragma HLS pipeline
////#pragma HLS unroll
//			if (k < N) X[k].read_nb(x[k]);
//			if (k < M) Y[k].read_nb(y[k]);
//			DIn[k].read(pt[k]);
//		}
//		 Xele[0].read_nb(xele);
//	}

	row_fill: for (R = 0; R < M; ++R) {
#pragma HLS loop_flatten off

		computeDelay(Xele/*xele*/, X/*x*/, Y/*y*/, R, Delay);

		computeAcc(DIn/*pt*/, Y/*y*/, Delay, R, DOut);
	}
}
