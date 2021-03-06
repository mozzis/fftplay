/*
 * Copyright (c) 2003, 2006 Matteo Frigo
 * Copyright (c) 2003, 2006 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Fri Jan 27 20:42:32 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 12 -dif -name hb_12 -include hb.h */

/*
 * This function contains 118 FP additions, 68 FP multiplications,
 * (or, 72 additions, 22 multiplications, 46 fused multiply/add),
 * 66 stack variables, and 48 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hb.h"

static const R *hb_12(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 22, MAKE_VOLATILE_STRIDE(ios)) {
	  E T18, T20, T21, T1b, T2a, T1s, T29, T1p, Tz, T11, TD, Tb, Tg, T23, T1f;
	  E Tl, TN, TI, T1i, T24, T1z, T2d, T1w, T2c;
	  {
	       E T5, Ta, Ty, Tt;
	       {
		    E T1, Tp, T6, Tu, T7, T1o, T4, T17, Ts, T8, Tv, Tw;
		    T1 = rio[0];
		    Tp = iio[0];
		    T6 = iio[-WS(ios, 6)];
		    Tu = rio[WS(ios, 6)];
		    {
			 E T2, T3, Tq, Tr;
			 T2 = rio[WS(ios, 4)];
			 T3 = iio[-WS(ios, 8)];
			 Tq = iio[-WS(ios, 4)];
			 Tr = rio[WS(ios, 8)];
			 T7 = iio[-WS(ios, 10)];
			 T1o = T2 - T3;
			 T4 = T2 + T3;
			 T17 = Tr + Tq;
			 Ts = Tq - Tr;
			 T8 = rio[WS(ios, 2)];
			 Tv = rio[WS(ios, 10)];
			 Tw = iio[-WS(ios, 2)];
		    }
		    {
			 E T1r, T1a, T19, T1q, T9, Tx, T16, T1n;
			 T5 = T1 + T4;
			 T16 = FNMS(KP500000000, T4, T1);
			 T1r = T7 - T8;
			 T9 = T7 + T8;
			 T1a = Tv + Tw;
			 Tx = Tv - Tw;
			 T18 = FMA(KP866025403, T17, T16);
			 T20 = FNMS(KP866025403, T17, T16);
			 T19 = FNMS(KP500000000, T9, T6);
			 Ta = T6 + T9;
			 Ty = Tu + Tx;
			 T1q = FNMS(KP500000000, Tx, Tu);
			 T1n = FNMS(KP500000000, Ts, Tp);
			 Tt = Tp + Ts;
			 T21 = FMA(KP866025403, T1a, T19);
			 T1b = FNMS(KP866025403, T1a, T19);
			 T2a = FNMS(KP866025403, T1r, T1q);
			 T1s = FMA(KP866025403, T1r, T1q);
			 T29 = FMA(KP866025403, T1o, T1n);
			 T1p = FNMS(KP866025403, T1o, T1n);
		    }
	       }
	       {
		    E Tc, TE, Th, TM, Ti, Tf, T1v, TH, T1e, Tj, TJ, TK;
		    Tc = rio[WS(ios, 3)];
		    Tz = Tt + Ty;
		    T11 = Tt - Ty;
		    TE = iio[-WS(ios, 3)];
		    TD = T5 - Ta;
		    Tb = T5 + Ta;
		    Th = iio[-WS(ios, 9)];
		    TM = rio[WS(ios, 9)];
		    {
			 E Td, Te, TF, TG;
			 Td = iio[-WS(ios, 7)];
			 Te = iio[-WS(ios, 11)];
			 TF = rio[WS(ios, 7)];
			 TG = rio[WS(ios, 11)];
			 Ti = rio[WS(ios, 1)];
			 Tf = Td + Te;
			 T1v = Td - Te;
			 TH = TF + TG;
			 T1e = TF - TG;
			 Tj = rio[WS(ios, 5)];
			 TJ = iio[-WS(ios, 5)];
			 TK = iio[-WS(ios, 1)];
		    }
		    {
			 E T1y, T1h, T1g, T1x, Tk, TL, T1d, T1u;
			 T1d = FNMS(KP500000000, Tf, Tc);
			 Tg = Tc + Tf;
			 Tk = Ti + Tj;
			 T1y = Ti - Tj;
			 TL = TJ + TK;
			 T1h = TJ - TK;
			 T23 = FMA(KP866025403, T1e, T1d);
			 T1f = FNMS(KP866025403, T1e, T1d);
			 Tl = Th + Tk;
			 T1g = FNMS(KP500000000, Tk, Th);
			 T1x = FMA(KP500000000, TL, TM);
			 TN = TL - TM;
			 TI = TE - TH;
			 T1u = FMA(KP500000000, TH, TE);
			 T1i = FNMS(KP866025403, T1h, T1g);
			 T24 = FMA(KP866025403, T1h, T1g);
			 T1z = FMA(KP866025403, T1y, T1x);
			 T2d = FNMS(KP866025403, T1y, T1x);
			 T1w = FNMS(KP866025403, T1v, T1u);
			 T2c = FMA(KP866025403, T1v, T1u);
		    }
	       }
	  }
	  {
	       E TY, T13, TX, T10;
	       {
		    E Tn, T12, TC, Tm, To, TS, TP, TO;
		    Tn = W[16];
		    T12 = TI + TN;
		    TO = TI - TN;
		    TC = W[17];
		    Tm = Tg + Tl;
		    To = Tg - Tl;
		    TS = TD + TO;
		    TP = TD - TO;
		    {
			 E TV, TU, TW, TT;
			 {
			      E TB, TR, TA, TQ;
			      TV = Tz - To;
			      TA = To + Tz;
			      rio[0] = Tb + Tm;
			      TQ = Tn * TP;
			      TB = Tn * TA;
			      TR = W[4];
			      rio[WS(ios, 9)] = FNMS(TC, TA, TQ);
			      TU = W[5];
			      iio[-WS(ios, 2)] = FMA(TC, TP, TB);
			      TW = TR * TV;
			      TT = TR * TS;
			 }
			 iio[-WS(ios, 8)] = FMA(TU, TS, TW);
			 rio[WS(ios, 3)] = FNMS(TU, TV, TT);
			 TY = Tb - Tm;
			 T13 = T11 - T12;
			 TX = W[10];
			 T10 = W[11];
			 iio[-WS(ios, 11)] = T11 + T12;
		    }
	       }
	       {
		    E T1c, T1A, T1t, T1j, T22, T2e, T2b, T2B, T2q, T25, T2s, T2y, T2C, T2z, T2w;
		    E T2A;
		    {
			 E T1X, T1M, T1O, T1U, T1Y, T1V, T1S, T1W, T1P, T1Q;
			 {
			      E T1K, TZ, T14, T1L;
			      T1c = T18 + T1b;
			      T1K = T18 - T1b;
			      TZ = TX * TY;
			      T14 = T10 * TY;
			      T1L = T1w + T1z;
			      T1A = T1w - T1z;
			      T1t = T1p - T1s;
			      T1P = T1p + T1s;
			      rio[WS(ios, 6)] = FNMS(T10, T13, TZ);
			      iio[-WS(ios, 5)] = FMA(TX, T13, T14);
			      T1X = T1K + T1L;
			      T1M = T1K - T1L;
			      T1Q = T1f - T1i;
			      T1j = T1f + T1i;
			 }
			 {
			      E T1J, T1T, T1R, T1N;
			      T1J = W[8];
			      T1O = W[9];
			      T1T = W[20];
			      T1U = T1P - T1Q;
			      T1R = T1P + T1Q;
			      T1N = T1J * T1M;
			      T1Y = T1T * T1X;
			      T1V = T1T * T1U;
			      T1S = T1J * T1R;
			      rio[WS(ios, 5)] = FNMS(T1O, T1R, T1N);
			      T1W = W[21];
			 }
			 {
			      E T2t, T2u, T2o, T2p;
			      T2o = T20 - T21;
			      T22 = T20 + T21;
			      iio[-WS(ios, 6)] = FMA(T1O, T1M, T1S);
			      T2p = T2c + T2d;
			      T2e = T2c - T2d;
			      rio[WS(ios, 11)] = FNMS(T1W, T1U, T1Y);
			      iio[0] = FMA(T1W, T1X, T1V);
			      T2b = T29 - T2a;
			      T2t = T29 + T2a;
			      T2B = T2o + T2p;
			      T2q = T2o - T2p;
			      T2u = T23 - T24;
			      T25 = T23 + T24;
			      {
				   E T2n, T2x, T2v, T2r;
				   T2n = W[0];
				   T2s = W[1];
				   T2x = W[12];
				   T2y = T2t - T2u;
				   T2v = T2t + T2u;
				   T2r = T2n * T2q;
				   T2C = T2x * T2B;
				   T2z = T2x * T2y;
				   T2w = T2n * T2v;
				   rio[WS(ios, 1)] = FNMS(T2s, T2v, T2r);
				   T2A = W[13];
			      }
			 }
		    }
		    {
			 E T2i, T2h, T2l, T2j, T2k, T26;
			 iio[-WS(ios, 10)] = FMA(T2s, T2q, T2w);
			 rio[WS(ios, 7)] = FNMS(T2A, T2y, T2C);
			 iio[-WS(ios, 4)] = FMA(T2A, T2B, T2z);
			 T2i = T22 + T25;
			 T26 = T22 - T25;
			 {
			      E T1Z, T28, T2f, T27, T2g;
			      T1Z = W[18];
			      T28 = W[19];
			      T2h = W[6];
			      T2l = T2b + T2e;
			      T2f = T2b - T2e;
			      T27 = T1Z * T26;
			      T2g = T28 * T26;
			      T2j = T2h * T2i;
			      T2k = W[7];
			      rio[WS(ios, 10)] = FNMS(T28, T2f, T27);
			      iio[-WS(ios, 1)] = FMA(T1Z, T2f, T2g);
			 }
			 {
			      E T1k, T1E, T1H, T1B, T2m, T15, T1m;
			      rio[WS(ios, 4)] = FNMS(T2k, T2l, T2j);
			      T2m = T2k * T2i;
			      iio[-WS(ios, 7)] = FMA(T2h, T2l, T2m);
			      T1k = T1c - T1j;
			      T1E = T1c + T1j;
			      T1H = T1t + T1A;
			      T1B = T1t - T1A;
			      T15 = W[2];
			      T1m = W[3];
			      {
				   E T1D, T1G, T1l, T1C, T1F, T1I;
				   T1D = W[14];
				   T1G = W[15];
				   T1l = T15 * T1k;
				   T1C = T1m * T1k;
				   T1F = T1D * T1E;
				   T1I = T1G * T1E;
				   rio[WS(ios, 2)] = FNMS(T1m, T1B, T1l);
				   iio[-WS(ios, 9)] = FMA(T15, T1B, T1C);
				   rio[WS(ios, 8)] = FNMS(T1G, T1H, T1F);
				   iio[-WS(ios, 3)] = FMA(T1D, T1H, T1I);
			      }
			 }
		    }
	       }
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 12},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 12, "hb_12", twinstr, &GENUS, {72, 22, 46, 0}, 0, 0, 0 };

void X(codelet_hb_12) (planner *p) {
     X(khc2hc_register) (p, hb_12, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc -compact -variables 4 -pipeline-latency 4 -sign 1 -n 12 -dif -name hb_12 -include hb.h */

/*
 * This function contains 118 FP additions, 60 FP multiplications,
 * (or, 88 additions, 30 multiplications, 30 fused multiply/add),
 * 39 stack variables, and 48 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hb.h"

static const R *hb_12(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 22, MAKE_VOLATILE_STRIDE(ios)) {
	  E T5, Tt, T12, T1M, T1i, T1U, Tl, TM, T1c, T1Y, T1s, T1Q, Ta, Ty, T15;
	  E T1N, T1l, T1V, Tg, TH, T19, T1X, T1p, T1P;
	  {
	       E T1, Tp, T4, T1g, Ts, T11, T10, T1h;
	       T1 = rio[0];
	       Tp = iio[0];
	       {
		    E T2, T3, Tq, Tr;
		    T2 = rio[WS(ios, 4)];
		    T3 = iio[-WS(ios, 8)];
		    T4 = T2 + T3;
		    T1g = KP866025403 * (T2 - T3);
		    Tq = rio[WS(ios, 8)];
		    Tr = iio[-WS(ios, 4)];
		    Ts = Tq - Tr;
		    T11 = KP866025403 * (Tq + Tr);
	       }
	       T5 = T1 + T4;
	       Tt = Tp - Ts;
	       T10 = FNMS(KP500000000, T4, T1);
	       T12 = T10 - T11;
	       T1M = T10 + T11;
	       T1h = FMA(KP500000000, Ts, Tp);
	       T1i = T1g + T1h;
	       T1U = T1h - T1g;
	  }
	  {
	       E Th, TL, Tk, T1a, TK, T1r, T1b, T1q;
	       Th = iio[-WS(ios, 9)];
	       TL = rio[WS(ios, 9)];
	       {
		    E Ti, Tj, TI, TJ;
		    Ti = rio[WS(ios, 1)];
		    Tj = rio[WS(ios, 5)];
		    Tk = Ti + Tj;
		    T1a = KP866025403 * (Ti - Tj);
		    TI = iio[-WS(ios, 5)];
		    TJ = iio[-WS(ios, 1)];
		    TK = TI + TJ;
		    T1r = KP866025403 * (TI - TJ);
	       }
	       Tl = Th + Tk;
	       TM = TK - TL;
	       T1b = FMA(KP500000000, TK, TL);
	       T1c = T1a - T1b;
	       T1Y = T1a + T1b;
	       T1q = FNMS(KP500000000, Tk, Th);
	       T1s = T1q + T1r;
	       T1Q = T1q - T1r;
	  }
	  {
	       E T6, Tx, T9, T1j, Tw, T14, T13, T1k;
	       T6 = iio[-WS(ios, 6)];
	       Tx = rio[WS(ios, 6)];
	       {
		    E T7, T8, Tu, Tv;
		    T7 = iio[-WS(ios, 10)];
		    T8 = rio[WS(ios, 2)];
		    T9 = T7 + T8;
		    T1j = KP866025403 * (T7 - T8);
		    Tu = rio[WS(ios, 10)];
		    Tv = iio[-WS(ios, 2)];
		    Tw = Tu - Tv;
		    T14 = KP866025403 * (Tu + Tv);
	       }
	       Ta = T6 + T9;
	       Ty = Tw + Tx;
	       T13 = FNMS(KP500000000, T9, T6);
	       T15 = T13 + T14;
	       T1N = T13 - T14;
	       T1k = FMS(KP500000000, Tw, Tx);
	       T1l = T1j + T1k;
	       T1V = T1k - T1j;
	  }
	  {
	       E Tc, TD, Tf, T17, TG, T1o, T18, T1n;
	       Tc = rio[WS(ios, 3)];
	       TD = iio[-WS(ios, 3)];
	       {
		    E Td, Te, TE, TF;
		    Td = iio[-WS(ios, 7)];
		    Te = iio[-WS(ios, 11)];
		    Tf = Td + Te;
		    T17 = KP866025403 * (Td - Te);
		    TE = rio[WS(ios, 7)];
		    TF = rio[WS(ios, 11)];
		    TG = TE + TF;
		    T1o = KP866025403 * (TE - TF);
	       }
	       Tg = Tc + Tf;
	       TH = TD - TG;
	       T18 = FMA(KP500000000, TG, TD);
	       T19 = T17 + T18;
	       T1X = T18 - T17;
	       T1n = FNMS(KP500000000, Tf, Tc);
	       T1p = T1n + T1o;
	       T1P = T1n - T1o;
	  }
	  {
	       E Tb, Tm, TU, TW, TX, TY, TT, TV;
	       Tb = T5 + Ta;
	       Tm = Tg + Tl;
	       TU = Tb - Tm;
	       TW = Tt - Ty;
	       TX = TH + TM;
	       TY = TW - TX;
	       rio[0] = Tb + Tm;
	       iio[-WS(ios, 11)] = TW + TX;
	       TT = W[10];
	       TV = W[11];
	       rio[WS(ios, 6)] = FNMS(TV, TY, TT * TU);
	       iio[-WS(ios, 5)] = FMA(TV, TU, TT * TY);
	  }
	  {
	       E T28, T2g, T2c, T2e;
	       {
		    E T26, T27, T2a, T2b;
		    T26 = T1M - T1N;
		    T27 = T1X + T1Y;
		    T28 = T26 - T27;
		    T2g = T26 + T27;
		    T2a = T1U - T1V;
		    T2b = T1P - T1Q;
		    T2c = T2a + T2b;
		    T2e = T2a - T2b;
	       }
	       {
		    E T25, T29, T2d, T2f;
		    T25 = W[8];
		    T29 = W[9];
		    rio[WS(ios, 5)] = FNMS(T29, T2c, T25 * T28);
		    iio[-WS(ios, 6)] = FMA(T25, T2c, T29 * T28);
		    T2d = W[20];
		    T2f = W[21];
		    iio[0] = FMA(T2d, T2e, T2f * T2g);
		    rio[WS(ios, 11)] = FNMS(T2f, T2e, T2d * T2g);
	       }
	  }
	  {
	       E TA, TS, TO, TQ;
	       {
		    E To, Tz, TC, TN;
		    To = Tg - Tl;
		    Tz = Tt + Ty;
		    TA = To + Tz;
		    TS = Tz - To;
		    TC = T5 - Ta;
		    TN = TH - TM;
		    TO = TC - TN;
		    TQ = TC + TN;
	       }
	       {
		    E Tn, TB, TP, TR;
		    Tn = W[16];
		    TB = W[17];
		    iio[-WS(ios, 2)] = FMA(Tn, TA, TB * TO);
		    rio[WS(ios, 9)] = FNMS(TB, TA, Tn * TO);
		    TP = W[4];
		    TR = W[5];
		    rio[WS(ios, 3)] = FNMS(TR, TS, TP * TQ);
		    iio[-WS(ios, 8)] = FMA(TP, TS, TR * TQ);
	       }
	  }
	  {
	       E T1S, T22, T20, T24;
	       {
		    E T1O, T1R, T1W, T1Z;
		    T1O = T1M + T1N;
		    T1R = T1P + T1Q;
		    T1S = T1O - T1R;
		    T22 = T1O + T1R;
		    T1W = T1U + T1V;
		    T1Z = T1X - T1Y;
		    T20 = T1W - T1Z;
		    T24 = T1W + T1Z;
	       }
	       {
		    E T1L, T1T, T21, T23;
		    T1L = W[2];
		    T1T = W[3];
		    rio[WS(ios, 2)] = FNMS(T1T, T20, T1L * T1S);
		    iio[-WS(ios, 9)] = FMA(T1T, T1S, T1L * T20);
		    T21 = W[14];
		    T23 = W[15];
		    rio[WS(ios, 8)] = FNMS(T23, T24, T21 * T22);
		    iio[-WS(ios, 3)] = FMA(T23, T22, T21 * T24);
	       }
	  }
	  {
	       E T1C, T1I, T1G, T1K;
	       {
		    E T1A, T1B, T1E, T1F;
		    T1A = T12 + T15;
		    T1B = T1p + T1s;
		    T1C = T1A - T1B;
		    T1I = T1A + T1B;
		    T1E = T1i + T1l;
		    T1F = T19 + T1c;
		    T1G = T1E - T1F;
		    T1K = T1E + T1F;
	       }
	       {
		    E T1z, T1D, T1H, T1J;
		    T1z = W[18];
		    T1D = W[19];
		    rio[WS(ios, 10)] = FNMS(T1D, T1G, T1z * T1C);
		    iio[-WS(ios, 1)] = FMA(T1D, T1C, T1z * T1G);
		    T1H = W[6];
		    T1J = W[7];
		    rio[WS(ios, 4)] = FNMS(T1J, T1K, T1H * T1I);
		    iio[-WS(ios, 7)] = FMA(T1J, T1I, T1H * T1K);
	       }
	  }
	  {
	       E T1e, T1y, T1u, T1w;
	       {
		    E T16, T1d, T1m, T1t;
		    T16 = T12 - T15;
		    T1d = T19 - T1c;
		    T1e = T16 - T1d;
		    T1y = T16 + T1d;
		    T1m = T1i - T1l;
		    T1t = T1p - T1s;
		    T1u = T1m + T1t;
		    T1w = T1m - T1t;
	       }
	       {
		    E TZ, T1f, T1v, T1x;
		    TZ = W[0];
		    T1f = W[1];
		    rio[WS(ios, 1)] = FNMS(T1f, T1u, TZ * T1e);
		    iio[-WS(ios, 10)] = FMA(TZ, T1u, T1f * T1e);
		    T1v = W[12];
		    T1x = W[13];
		    iio[-WS(ios, 4)] = FMA(T1v, T1w, T1x * T1y);
		    rio[WS(ios, 7)] = FNMS(T1x, T1w, T1v * T1y);
	       }
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 12},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 12, "hb_12", twinstr, &GENUS, {88, 30, 30, 0}, 0, 0, 0 };

void X(codelet_hb_12) (planner *p) {
     X(khc2hc_register) (p, hb_12, &desc);
}
#endif				/* HAVE_FMA */
