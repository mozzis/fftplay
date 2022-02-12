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
/* Generated on Fri Jan 27 19:31:23 EST 2006 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -twiddle-log3 -precompute-twiddles -n 16 -name t2_16 -include t.h */

/*
 * This function contains 196 FP additions, 134 FP multiplications,
 * (or, 104 additions, 42 multiplications, 92 fused multiply/add),
 * 116 stack variables, and 64 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle.ml,v 1.23 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t.h"

static const R *t2_16(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     for (i = m; i > 0; i = i - 1, ri = ri + dist, ii = ii + dist, W = W + 8, MAKE_VOLATILE_STRIDE(ios)) {
	  E T3S, T3R;
	  {
	       E T3, T6, Tf, TM, T2, Th, T4, Tp, Tt, Ta, Tz, Ti, Tm, TC, TO;
	       E T1e, T1C, T1i, TT, TP, T1y, T5;
	       {
		    E TS, TN, Tl, Tg;
		    T3 = W[0];
		    T6 = W[1];
		    Tf = W[2];
		    TM = W[6];
		    T2 = W[4];
		    Th = W[3];
		    Tl = Tf * T6;
		    Tg = Tf * T3;
		    TS = TM * T6;
		    TN = TM * T3;
		    T4 = T2 * T3;
		    Tp = T2 * Tf;
		    Tt = T2 * Th;
		    Ta = T2 * T6;
		    Tz = FMA(Th, T6, Tg);
		    Ti = FNMS(Th, T6, Tg);
		    Tm = FMA(Th, T3, Tl);
		    TC = FNMS(Th, T3, Tl);
		    TO = W[7];
		    T1e = T2 * Ti;
		    T1C = T2 * TC;
		    T1i = T2 * Tm;
		    TT = FNMS(TO, T3, TS);
		    TP = FMA(TO, T6, TN);
		    T1y = T2 * Tz;
		    T5 = W[5];
	       }
	       {
		    E TW, TZ, T3L, T3A, T1U, Te, T2D, T1G, T3h, T2A, T3i, T2I, T2B, T1R, T3M;
		    E Tx, T3w, T1Z, TL, T26, T25, T37, T1d, T2o, T2l, T3c, T10, T2m, T1s, T3d;
		    E T2t, TX, TV, T2a;
		    {
			 E T1m, T1f, T1p, T1j, TI, Tu, TF, Tr, Tv, To, T1W, Ts, T1X;
			 {
			      E Tq, T1L, T1O, T1x, T1F, T2x, T2z;
			      {
				   E T1, T3z, Tc, Tb, T7, T1z, T1D, T8;
				   T1 = ri[0];
				   T3z = ii[0];
				   Tc = ii[WS(ios, 8)];
				   T1m = FNMS(T5, Tm, T1e);
				   T1f = FMA(T5, Tm, T1e);
				   T1p = FMA(T5, Ti, T1i);
				   T1j = FNMS(T5, Ti, T1i);
				   TW = FMA(T5, Th, Tp);
				   Tq = FNMS(T5, Th, Tp);
				   Tb = FNMS(T5, T3, Ta);
				   TI = FMA(T5, T3, Ta);
				   TZ = FNMS(T5, Tf, Tt);
				   Tu = FMA(T5, Tf, Tt);
				   T7 = FMA(T5, T6, T4);
				   TF = FNMS(T5, T6, T4);
				   T1L = FMA(T5, TC, T1y);
				   T1z = FNMS(T5, TC, T1y);
				   T1O = FNMS(T5, Tz, T1C);
				   T1D = FMA(T5, Tz, T1C);
				   T8 = ri[WS(ios, 8)];
				   {
					E T1u, T1w, T1A, T1v, T2w, T1B, T1E, T2y;
					T1u = ri[WS(ios, 15)];
					T1w = ii[WS(ios, 15)];
					{
					     E T3y, Td, T3x, T9;
					     T1A = ri[WS(ios, 7)];
					     T3x = Tb * T8;
					     T9 = T7 * T8;
					     T1v = TM * T1u;
					     T2w = TM * T1w;
					     T3y = FMA(T7, Tc, T3x);
					     Td = FNMS(Tb, Tc, T9);
					     T1B = T1z * T1A;
					     T1E = ii[WS(ios, 7)];
					     T3L = T3z - T3y;
					     T3A = T3y + T3z;
					     T1U = T1 - Td;
					     Te = T1 + Td;
					     T2y = T1z * T1E;
					}
					T1x = FMA(TO, T1w, T1v);
					T1F = FMA(T1D, T1E, T1B);
					T2x = FNMS(TO, T1u, T2w);
					T2z = FNMS(T1D, T1A, T2y);
				   }
			      }
			      {
				   E T1P, T1J, T1H, T1M;
				   T1P = ii[WS(ios, 11)];
				   T2D = T1x - T1F;
				   T1G = T1x + T1F;
				   T3h = T2x + T2z;
				   T2A = T2x - T2z;
				   T1J = ii[WS(ios, 3)];
				   T1H = ri[WS(ios, 3)];
				   T1M = ri[WS(ios, 11)];
				   {
					E Tj, Tk, Tn, T1V;
					{
					     E T2F, T1K, T2H, T1Q;
					     Tj = ri[WS(ios, 4)];
					     {
						  E T2E, T1I, T2G, T1N;
						  T2E = Tf * T1J;
						  T1I = Tf * T1H;
						  T2G = T1O * T1M;
						  T1N = T1L * T1M;
						  T2F = FNMS(Th, T1H, T2E);
						  T1K = FMA(Th, T1J, T1I);
						  T2H = FMA(T1L, T1P, T2G);
						  T1Q = FNMS(T1O, T1P, T1N);
						  Tk = Ti * Tj;
					     }
					     Tn = ii[WS(ios, 4)];
					     T3i = T2F + T2H;
					     T2I = T2F - T2H;
					     T2B = T1K - T1Q;
					     T1R = T1K + T1Q;
					     T1V = Ti * Tn;
					}
					Tr = ri[WS(ios, 12)];
					Tv = ii[WS(ios, 12)];
					To = FMA(Tm, Tn, Tk);
					T1W = FNMS(Tm, Tj, T1V);
					Ts = Tq * Tr;
					T1X = Tq * Tv;
				   }
			      }
			 }
			 {
			      E T18, T2i, T1c, T2k;
			      {
				   E TG, TH, TJ, TE, T22;
				   {
					E TD, T21, TB, TA, Tw, T1Y;
					TD = ii[WS(ios, 2)];
					TA = ri[WS(ios, 2)];
					Tw = FMA(Tu, Tv, Ts);
					T1Y = FNMS(Tu, Tr, T1X);
					TG = ri[WS(ios, 10)];
					T21 = TC * TA;
					TB = Tz * TA;
					T3M = To - Tw;
					Tx = To + Tw;
					T3w = T1W + T1Y;
					T1Z = T1W - T1Y;
					TH = TF * TG;
					TJ = ii[WS(ios, 10)];
					TE = FNMS(TC, TD, TB);
					T22 = FMA(Tz, TD, T21);
				   }
				   {
					E T15, T17, T16, T2h, T19, T1b, T1a, T2j, T24, TK, T23;
					T15 = ri[WS(ios, 1)];
					TK = FMA(TI, TJ, TH);
					T23 = TF * TJ;
					T17 = ii[WS(ios, 1)];
					T16 = T3 * T15;
					TL = TE + TK;
					T26 = TE - TK;
					T24 = FNMS(TI, TG, T23);
					T2h = T3 * T17;
					T19 = ri[WS(ios, 9)];
					T1b = ii[WS(ios, 9)];
					T25 = T22 - T24;
					T37 = T22 + T24;
					T1a = T2 * T19;
					T2j = T2 * T1b;
					T18 = FMA(T6, T17, T16);
					T2i = FNMS(T6, T15, T2h);
					T1c = FMA(T5, T1b, T1a);
					T2k = FNMS(T5, T19, T2j);
				   }
			      }
			      {
				   E T1n, T1q, T1l, T2q, T1o, T2r;
				   {
					E T1k, T1h, T2p, T1g;
					T1k = ii[WS(ios, 5)];
					T1g = ri[WS(ios, 5)];
					T1d = T18 + T1c;
					T2o = T18 - T1c;
					T2l = T2i - T2k;
					T3c = T2i + T2k;
					T1h = T1f * T1g;
					T2p = T1j * T1g;
					T1n = ri[WS(ios, 13)];
					T1q = ii[WS(ios, 13)];
					T1l = FNMS(T1j, T1k, T1h);
					T2q = FMA(T1f, T1k, T2p);
					T1o = T1m * T1n;
					T2r = T1m * T1q;
				   }
				   {
					E TU, T29, TR, TQ, T1r, T2s;
					TU = ii[WS(ios, 14)];
					TQ = ri[WS(ios, 14)];
					T1r = FMA(T1p, T1q, T1o);
					T2s = FNMS(T1p, T1n, T2r);
					T10 = ii[WS(ios, 6)];
					T29 = TT * TQ;
					TR = TP * TQ;
					T2m = T1l - T1r;
					T1s = T1l + T1r;
					T3d = T2q + T2s;
					T2t = T2q - T2s;
					TX = ri[WS(ios, 6)];
					TV = FNMS(TT, TU, TR);
					T2a = FMA(TP, TU, T29);
				   }
			      }
			 }
		    }
		    {
			 E T36, Ty, T3B, T3G, T3b, T3g, T3e, T3r, T2d, T28, T3D, T1T, T3v, T39, T3F;
			 E T13, T3s, T3j, T2b, TY;
			 T36 = Te - Tx;
			 Ty = Te + Tx;
			 T2b = TZ * TX;
			 TY = TW * TX;
			 T3B = T3w + T3A;
			 T3G = T3A - T3w;
			 {
			      E T1t, T2c, T11, T1S, T38, T12;
			      T3b = T1d - T1s;
			      T1t = T1d + T1s;
			      T2c = FMA(TW, T10, T2b);
			      T11 = FNMS(TZ, T10, TY);
			      T1S = T1G + T1R;
			      T3g = T1G - T1R;
			      T3e = T3c - T3d;
			      T3r = T3c + T3d;
			      T38 = T2a + T2c;
			      T2d = T2a - T2c;
			      T28 = TV - T11;
			      T12 = TV + T11;
			      T3D = T1S - T1t;
			      T1T = T1t + T1S;
			      T3v = T37 + T38;
			      T39 = T37 - T38;
			      T3F = T12 - TL;
			      T13 = TL + T12;
			      T3s = T3h + T3i;
			      T3j = T3h - T3i;
			 }
			 {
			      E T3m, T3a, T3J, T3H, T3E, T3C;
			      T3E = T3B - T3v;
			      T3C = T3v + T3B;
			      {
				   E T3q, T14, T3u, T3t;
				   T3q = Ty - T13;
				   T14 = Ty + T13;
				   T3u = T3r + T3s;
				   T3t = T3r - T3s;
				   ii[WS(ios, 4)] = T3D + T3E;
				   ii[WS(ios, 12)] = T3E - T3D;
				   ri[0] = T14 + T1T;
				   ri[WS(ios, 8)] = T14 - T1T;
				   ii[0] = T3u + T3C;
				   ii[WS(ios, 8)] = T3C - T3u;
				   ri[WS(ios, 4)] = T3q + T3t;
				   ri[WS(ios, 12)] = T3q - T3t;
			      }
			      T3m = T36 - T39;
			      T3a = T36 + T39;
			      T3J = T3G - T3F;
			      T3H = T3F + T3G;
			      {
				   E T2Q, T20, T3N, T3T, T2J, T2C, T3O, T2f, T34, T30, T2W, T2V, T3U, T2T, T2N;
				   E T2v;
				   {
					E T2R, T27, T2e, T2S;
					{
					     E T3n, T3f, T3o, T3k;
					     T2Q = T1U + T1Z;
					     T20 = T1U - T1Z;
					     T3n = T3e - T3b;
					     T3f = T3b + T3e;
					     T3o = T3g + T3j;
					     T3k = T3g - T3j;
					     T3N = T3L - T3M;
					     T3T = T3M + T3L;
					     {
						  E T3p, T3I, T3K, T3l;
						  T3p = T3n - T3o;
						  T3I = T3n + T3o;
						  T3K = T3k - T3f;
						  T3l = T3f + T3k;
						  ri[WS(ios, 6)] = FMA(KP707106781, T3p, T3m);
						  ri[WS(ios, 14)] = FNMS(KP707106781, T3p, T3m);
						  ii[WS(ios, 10)] = FNMS(KP707106781, T3I, T3H);
						  ii[WS(ios, 2)] = FMA(KP707106781, T3I, T3H);
						  ii[WS(ios, 14)] = FNMS(KP707106781, T3K, T3J);
						  ii[WS(ios, 6)] = FMA(KP707106781, T3K, T3J);
						  ri[WS(ios, 2)] = FMA(KP707106781, T3l, T3a);
						  ri[WS(ios, 10)] = FNMS(KP707106781, T3l, T3a);
						  T2R = T26 + T25;
						  T27 = T25 - T26;
						  T2e = T28 + T2d;
						  T2S = T28 - T2d;
					     }
					}
					{
					     E T2Y, T2Z, T2n, T2u;
					     T2J = T2D - T2I;
					     T2Y = T2D + T2I;
					     T2Z = T2A - T2B;
					     T2C = T2A + T2B;
					     T3O = T27 + T2e;
					     T2f = T27 - T2e;
					     T34 = FMA(KP414213562, T2Y, T2Z);
					     T30 = FNMS(KP414213562, T2Z, T2Y);
					     T2W = T2l - T2m;
					     T2n = T2l + T2m;
					     T2u = T2o - T2t;
					     T2V = T2o + T2t;
					     T3U = T2S - T2R;
					     T2T = T2R + T2S;
					     T2N = FNMS(KP414213562, T2n, T2u);
					     T2v = FMA(KP414213562, T2u, T2n);
					}
				   }
				   {
					E T33, T2X, T3X, T3Y;
					{
					     E T2M, T2g, T2O, T2K, T3V, T3W, T2P, T2L;
					     T2M = FNMS(KP707106781, T2f, T20);
					     T2g = FMA(KP707106781, T2f, T20);
					     T33 = FNMS(KP414213562, T2V, T2W);
					     T2X = FMA(KP414213562, T2W, T2V);
					     T2O = FMA(KP414213562, T2C, T2J);
					     T2K = FNMS(KP414213562, T2J, T2C);
					     T3V = FMA(KP707106781, T3U, T3T);
					     T3X = FNMS(KP707106781, T3U, T3T);
					     T3W = T2O - T2N;
					     T2P = T2N + T2O;
					     T3Y = T2v + T2K;
					     T2L = T2v - T2K;
					     ii[WS(ios, 11)] = FNMS(KP923879532, T3W, T3V);
					     ii[WS(ios, 3)] = FMA(KP923879532, T3W, T3V);
					     ri[WS(ios, 3)] = FMA(KP923879532, T2L, T2g);
					     ri[WS(ios, 11)] = FNMS(KP923879532, T2L, T2g);
					     ri[WS(ios, 15)] = FMA(KP923879532, T2P, T2M);
					     ri[WS(ios, 7)] = FNMS(KP923879532, T2P, T2M);
					}
					{
					     E T32, T3P, T3Q, T35, T2U, T31;
					     T32 = FNMS(KP707106781, T2T, T2Q);
					     T2U = FMA(KP707106781, T2T, T2Q);
					     T31 = T2X + T30;
					     T3S = T30 - T2X;
					     T3R = FNMS(KP707106781, T3O, T3N);
					     T3P = FMA(KP707106781, T3O, T3N);
					     ii[WS(ios, 15)] = FMA(KP923879532, T3Y, T3X);
					     ii[WS(ios, 7)] = FNMS(KP923879532, T3Y, T3X);
					     ri[WS(ios, 1)] = FMA(KP923879532, T31, T2U);
					     ri[WS(ios, 9)] = FNMS(KP923879532, T31, T2U);
					     T3Q = T33 + T34;
					     T35 = T33 - T34;
					     ii[WS(ios, 9)] = FNMS(KP923879532, T3Q, T3P);
					     ii[WS(ios, 1)] = FMA(KP923879532, T3Q, T3P);
					     ri[WS(ios, 5)] = FMA(KP923879532, T35, T32);
					     ri[WS(ios, 13)] = FNMS(KP923879532, T35, T32);
					}
				   }
			      }
			 }
		    }
	       }
	  }
	  ii[WS(ios, 13)] = FNMS(KP923879532, T3S, T3R);
	  ii[WS(ios, 5)] = FMA(KP923879532, T3S, T3R);
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 0, 1},
     {TW_CEXP, 0, 3},
     {TW_CEXP, 0, 9},
     {TW_CEXP, 0, 15},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 16, "t2_16", twinstr, &GENUS, {104, 42, 92, 0}, 0, 0, 0 };

void X(codelet_t2_16) (planner *p) {
     X(kdft_dit_register) (p, t2_16, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle -compact -variables 4 -pipeline-latency 4 -twiddle-log3 -precompute-twiddles -n 16 -name t2_16 -include t.h */

/*
 * This function contains 196 FP additions, 108 FP multiplications,
 * (or, 156 additions, 68 multiplications, 40 fused multiply/add),
 * 82 stack variables, and 64 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle.ml,v 1.23 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t.h"

static const R *t2_16(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     DK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     for (i = m; i > 0; i = i - 1, ri = ri + dist, ii = ii + dist, W = W + 8, MAKE_VOLATILE_STRIDE(ios)) {
	  E T3, T6, Tg, Ti, Tk, To, TE, TC, T5, T2, T8, TW, TJ, Tt, TU;
	  E Tc, Tx, TH, TN, TO, TP, TR, T1f, T1k, T1b, T1i, T1y, T1H, T1u, T1F;
	  {
	       E T7, Tw, Tb, Ts, T4, Tv, Ta, Tr;
	       {
		    E Th, Tn, Tj, Tm;
		    T3 = W[0];
		    T6 = W[1];
		    Tg = W[2];
		    Ti = W[3];
		    Th = Tg * T3;
		    Tn = Ti * T3;
		    Tj = Ti * T6;
		    Tm = Tg * T6;
		    Tk = Th - Tj;
		    To = Tm + Tn;
		    TE = Tm - Tn;
		    TC = Th + Tj;
		    T5 = W[5];
		    T7 = T5 * T6;
		    Tw = T5 * Tg;
		    Tb = T5 * T3;
		    Ts = T5 * Ti;
		    T2 = W[4];
		    T4 = T2 * T3;
		    Tv = T2 * Ti;
		    Ta = T2 * T6;
		    Tr = T2 * Tg;
	       }
	       T8 = T4 + T7;
	       TW = Tv - Tw;
	       TJ = Ta + Tb;
	       Tt = Tr - Ts;
	       TU = Tr + Ts;
	       Tc = Ta - Tb;
	       Tx = Tv + Tw;
	       TH = T4 - T7;
	       TN = W[6];
	       TO = W[7];
	       TP = FMA(TN, T3, TO * T6);
	       TR = FNMS(TO, T3, TN * T6);
	       {
		    E T1d, T1e, T19, T1a;
		    T1d = T2 * To;
		    T1e = T5 * Tk;
		    T1f = T1d - T1e;
		    T1k = T1d + T1e;
		    T19 = T2 * Tk;
		    T1a = T5 * To;
		    T1b = T19 + T1a;
		    T1i = T19 - T1a;
	       }
	       {
		    E T1w, T1x, T1s, T1t;
		    T1w = T2 * TE;
		    T1x = T5 * TC;
		    T1y = T1w + T1x;
		    T1H = T1w - T1x;
		    T1s = T2 * TC;
		    T1t = T5 * TE;
		    T1u = T1s - T1t;
		    T1F = T1s + T1t;
	       }
	  }
	  {
	       E Tf, T3r, T1N, T3e, TA, T3s, T1Q, T3b, TM, T2M, T1W, T2w, TZ, T2N, T21;
	       E T2x, T1B, T1K, T2V, T2W, T2X, T2Y, T2j, T2D, T2o, T2E, T18, T1n, T2Q, T2R;
	       E T2S, T2T, T28, T2A, T2d, T2B;
	       {
		    E T1, T3d, Te, T3c, T9, Td;
		    T1 = ri[0];
		    T3d = ii[0];
		    T9 = ri[WS(ios, 8)];
		    Td = ii[WS(ios, 8)];
		    Te = FNMS(Tc, Td, T8 * T9);
		    T3c = FMA(Tc, T9, T8 * Td);
		    Tf = T1 + Te;
		    T3r = T3d - T3c;
		    T1N = T1 - Te;
		    T3e = T3c + T3d;
	       }
	       {
		    E Tq, T1O, Tz, T1P;
		    {
			 E Tl, Tp, Tu, Ty;
			 Tl = ri[WS(ios, 4)];
			 Tp = ii[WS(ios, 4)];
			 Tq = FMA(Tk, Tl, To * Tp);
			 T1O = FNMS(To, Tl, Tk * Tp);
			 Tu = ri[WS(ios, 12)];
			 Ty = ii[WS(ios, 12)];
			 Tz = FMA(Tt, Tu, Tx * Ty);
			 T1P = FNMS(Tx, Tu, Tt * Ty);
		    }
		    TA = Tq + Tz;
		    T3s = Tq - Tz;
		    T1Q = T1O - T1P;
		    T3b = T1O + T1P;
	       }
	       {
		    E TG, T1S, TL, T1T, T1U, T1V;
		    {
			 E TD, TF, TI, TK;
			 TD = ri[WS(ios, 2)];
			 TF = ii[WS(ios, 2)];
			 TG = FNMS(TE, TF, TC * TD);
			 T1S = FMA(TE, TD, TC * TF);
			 TI = ri[WS(ios, 10)];
			 TK = ii[WS(ios, 10)];
			 TL = FMA(TH, TI, TJ * TK);
			 T1T = FNMS(TJ, TI, TH * TK);
		    }
		    TM = TG + TL;
		    T2M = T1S + T1T;
		    T1U = T1S - T1T;
		    T1V = TG - TL;
		    T1W = T1U - T1V;
		    T2w = T1V + T1U;
	       }
	       {
		    E TT, T1Y, TY, T1Z, T1X, T20;
		    {
			 E TQ, TS, TV, TX;
			 TQ = ri[WS(ios, 14)];
			 TS = ii[WS(ios, 14)];
			 TT = FNMS(TR, TS, TP * TQ);
			 T1Y = FMA(TR, TQ, TP * TS);
			 TV = ri[WS(ios, 6)];
			 TX = ii[WS(ios, 6)];
			 TY = FNMS(TW, TX, TU * TV);
			 T1Z = FMA(TW, TV, TU * TX);
		    }
		    TZ = TT + TY;
		    T2N = T1Y + T1Z;
		    T1X = TT - TY;
		    T20 = T1Y - T1Z;
		    T21 = T1X + T20;
		    T2x = T1X - T20;
	       }
	       {
		    E T1r, T2k, T1J, T2h, T1A, T2l, T1E, T2g;
		    {
			 E T1p, T1q, T1G, T1I;
			 T1p = ri[WS(ios, 15)];
			 T1q = ii[WS(ios, 15)];
			 T1r = FMA(TN, T1p, TO * T1q);
			 T2k = FNMS(TO, T1p, TN * T1q);
			 T1G = ri[WS(ios, 11)];
			 T1I = ii[WS(ios, 11)];
			 T1J = FNMS(T1H, T1I, T1F * T1G);
			 T2h = FMA(T1H, T1G, T1F * T1I);
		    }
		    {
			 E T1v, T1z, T1C, T1D;
			 T1v = ri[WS(ios, 7)];
			 T1z = ii[WS(ios, 7)];
			 T1A = FMA(T1u, T1v, T1y * T1z);
			 T2l = FNMS(T1y, T1v, T1u * T1z);
			 T1C = ri[WS(ios, 3)];
			 T1D = ii[WS(ios, 3)];
			 T1E = FMA(Tg, T1C, Ti * T1D);
			 T2g = FNMS(Ti, T1C, Tg * T1D);
		    }
		    T1B = T1r + T1A;
		    T1K = T1E + T1J;
		    T2V = T1B - T1K;
		    T2W = T2k + T2l;
		    T2X = T2g + T2h;
		    T2Y = T2W - T2X;
		    {
			 E T2f, T2i, T2m, T2n;
			 T2f = T1r - T1A;
			 T2i = T2g - T2h;
			 T2j = T2f - T2i;
			 T2D = T2f + T2i;
			 T2m = T2k - T2l;
			 T2n = T1E - T1J;
			 T2o = T2m + T2n;
			 T2E = T2m - T2n;
		    }
	       }
	       {
		    E T14, T24, T1m, T2b, T17, T25, T1h, T2a;
		    {
			 E T12, T13, T1j, T1l;
			 T12 = ri[WS(ios, 1)];
			 T13 = ii[WS(ios, 1)];
			 T14 = FMA(T3, T12, T6 * T13);
			 T24 = FNMS(T6, T12, T3 * T13);
			 T1j = ri[WS(ios, 13)];
			 T1l = ii[WS(ios, 13)];
			 T1m = FMA(T1i, T1j, T1k * T1l);
			 T2b = FNMS(T1k, T1j, T1i * T1l);
		    }
		    {
			 E T15, T16, T1c, T1g;
			 T15 = ri[WS(ios, 9)];
			 T16 = ii[WS(ios, 9)];
			 T17 = FMA(T2, T15, T5 * T16);
			 T25 = FNMS(T5, T15, T2 * T16);
			 T1c = ri[WS(ios, 5)];
			 T1g = ii[WS(ios, 5)];
			 T1h = FNMS(T1f, T1g, T1b * T1c);
			 T2a = FMA(T1f, T1c, T1b * T1g);
		    }
		    T18 = T14 + T17;
		    T1n = T1h + T1m;
		    T2Q = T18 - T1n;
		    T2R = T24 + T25;
		    T2S = T2a + T2b;
		    T2T = T2R - T2S;
		    {
			 E T26, T27, T29, T2c;
			 T26 = T24 - T25;
			 T27 = T1h - T1m;
			 T28 = T26 + T27;
			 T2A = T26 - T27;
			 T29 = T14 - T17;
			 T2c = T2a - T2b;
			 T2d = T29 - T2c;
			 T2B = T29 + T2c;
		    }
	       }
	       {
		    E T23, T2r, T3A, T3C, T2q, T3B, T2u, T3x;
		    {
			 E T1R, T22, T3y, T3z;
			 T1R = T1N - T1Q;
			 T22 = KP707106781 * (T1W - T21);
			 T23 = T1R + T22;
			 T2r = T1R - T22;
			 T3y = KP707106781 * (T2x - T2w);
			 T3z = T3s + T3r;
			 T3A = T3y + T3z;
			 T3C = T3z - T3y;
		    }
		    {
			 E T2e, T2p, T2s, T2t;
			 T2e = FMA(KP923879532, T28, KP382683432 * T2d);
			 T2p = FNMS(KP923879532, T2o, KP382683432 * T2j);
			 T2q = T2e + T2p;
			 T3B = T2p - T2e;
			 T2s = FNMS(KP923879532, T2d, KP382683432 * T28);
			 T2t = FMA(KP382683432, T2o, KP923879532 * T2j);
			 T2u = T2s - T2t;
			 T3x = T2s + T2t;
		    }
		    ri[WS(ios, 11)] = T23 - T2q;
		    ii[WS(ios, 11)] = T3A - T3x;
		    ri[WS(ios, 3)] = T23 + T2q;
		    ii[WS(ios, 3)] = T3x + T3A;
		    ri[WS(ios, 15)] = T2r - T2u;
		    ii[WS(ios, 15)] = T3C - T3B;
		    ri[WS(ios, 7)] = T2r + T2u;
		    ii[WS(ios, 7)] = T3B + T3C;
	       }
	       {
		    E T2P, T31, T3m, T3o, T30, T3n, T34, T3j;
		    {
			 E T2L, T2O, T3k, T3l;
			 T2L = Tf - TA;
			 T2O = T2M - T2N;
			 T2P = T2L + T2O;
			 T31 = T2L - T2O;
			 T3k = TZ - TM;
			 T3l = T3e - T3b;
			 T3m = T3k + T3l;
			 T3o = T3l - T3k;
		    }
		    {
			 E T2U, T2Z, T32, T33;
			 T2U = T2Q + T2T;
			 T2Z = T2V - T2Y;
			 T30 = KP707106781 * (T2U + T2Z);
			 T3n = KP707106781 * (T2Z - T2U);
			 T32 = T2T - T2Q;
			 T33 = T2V + T2Y;
			 T34 = KP707106781 * (T32 - T33);
			 T3j = KP707106781 * (T32 + T33);
		    }
		    ri[WS(ios, 10)] = T2P - T30;
		    ii[WS(ios, 10)] = T3m - T3j;
		    ri[WS(ios, 2)] = T2P + T30;
		    ii[WS(ios, 2)] = T3j + T3m;
		    ri[WS(ios, 14)] = T31 - T34;
		    ii[WS(ios, 14)] = T3o - T3n;
		    ri[WS(ios, 6)] = T31 + T34;
		    ii[WS(ios, 6)] = T3n + T3o;
	       }
	       {
		    E T2z, T2H, T3u, T3w, T2G, T3v, T2K, T3p;
		    {
			 E T2v, T2y, T3q, T3t;
			 T2v = T1N + T1Q;
			 T2y = KP707106781 * (T2w + T2x);
			 T2z = T2v + T2y;
			 T2H = T2v - T2y;
			 T3q = KP707106781 * (T1W + T21);
			 T3t = T3r - T3s;
			 T3u = T3q + T3t;
			 T3w = T3t - T3q;
		    }
		    {
			 E T2C, T2F, T2I, T2J;
			 T2C = FMA(KP382683432, T2A, KP923879532 * T2B);
			 T2F = FNMS(KP382683432, T2E, KP923879532 * T2D);
			 T2G = T2C + T2F;
			 T3v = T2F - T2C;
			 T2I = FNMS(KP382683432, T2B, KP923879532 * T2A);
			 T2J = FMA(KP923879532, T2E, KP382683432 * T2D);
			 T2K = T2I - T2J;
			 T3p = T2I + T2J;
		    }
		    ri[WS(ios, 9)] = T2z - T2G;
		    ii[WS(ios, 9)] = T3u - T3p;
		    ri[WS(ios, 1)] = T2z + T2G;
		    ii[WS(ios, 1)] = T3p + T3u;
		    ri[WS(ios, 13)] = T2H - T2K;
		    ii[WS(ios, 13)] = T3w - T3v;
		    ri[WS(ios, 5)] = T2H + T2K;
		    ii[WS(ios, 5)] = T3v + T3w;
	       }
	       {
		    E T11, T35, T3g, T3i, T1M, T3h, T38, T39;
		    {
			 E TB, T10, T3a, T3f;
			 TB = Tf + TA;
			 T10 = TM + TZ;
			 T11 = TB + T10;
			 T35 = TB - T10;
			 T3a = T2M + T2N;
			 T3f = T3b + T3e;
			 T3g = T3a + T3f;
			 T3i = T3f - T3a;
		    }
		    {
			 E T1o, T1L, T36, T37;
			 T1o = T18 + T1n;
			 T1L = T1B + T1K;
			 T1M = T1o + T1L;
			 T3h = T1L - T1o;
			 T36 = T2R + T2S;
			 T37 = T2W + T2X;
			 T38 = T36 - T37;
			 T39 = T36 + T37;
		    }
		    ri[WS(ios, 8)] = T11 - T1M;
		    ii[WS(ios, 8)] = T3g - T39;
		    ri[0] = T11 + T1M;
		    ii[0] = T39 + T3g;
		    ri[WS(ios, 12)] = T35 - T38;
		    ii[WS(ios, 12)] = T3i - T3h;
		    ri[WS(ios, 4)] = T35 + T38;
		    ii[WS(ios, 4)] = T3h + T3i;
	       }
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_CEXP, 0, 1},
     {TW_CEXP, 0, 3},
     {TW_CEXP, 0, 9},
     {TW_CEXP, 0, 15},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 16, "t2_16", twinstr, &GENUS, {156, 68, 40, 0}, 0, 0, 0 };

void X(codelet_t2_16) (planner *p) {
     X(kdft_dit_register) (p, t2_16, &desc);
}
#endif				/* HAVE_FMA */
