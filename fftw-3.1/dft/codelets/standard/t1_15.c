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
/* Generated on Fri Jan 27 19:25:30 EST 2006 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 15 -name t1_15 -include t.h */

/*
 * This function contains 184 FP additions, 140 FP multiplications,
 * (or, 72 additions, 28 multiplications, 112 fused multiply/add),
 * 89 stack variables, and 60 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle.ml,v 1.23 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t.h"

static const R *t1_15(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     DK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DK(KP618033988, +0.618033988749894848204586834365638117720309180);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = m; i > 0; i = i - 1, ri = ri + dist, ii = ii + dist, W = W + 28, MAKE_VOLATILE_STRIDE(ios)) {
	  E T2d, T2O, T2Q, T2m, T2k, T2l, T2P, T2n;
	  {
	       E T1G, T3u, T3k, T3t, T1B, Tf, T37, T1y, T2V, T2M, T2a, T2i, T39, Tz, T2X;
	       E T2t, T1O, T2e, T3a, TT, T10, T2Y, T2z, T1V, T2f, T2C, T12, T15, T14, T21;
	       E T1c, T1Y, T13;
	       {
		    E T2I, T1k, T1m, T1p, T1o, T28, T1w, T25, T1n;
		    {
			 E T1, T3j, T9, Tc, Tb, T1D, T7, T1E, Ta, T1j, T1i, T1h;
			 T1 = ri[0];
			 T3j = ii[0];
			 {
			      E T3, T6, T2, T5, T1C, T4, T8;
			      T3 = ri[WS(ios, 5)];
			      T6 = ii[WS(ios, 5)];
			      T2 = W[8];
			      T5 = W[9];
			      T9 = ri[WS(ios, 10)];
			      Tc = ii[WS(ios, 10)];
			      T1C = T2 * T6;
			      T4 = T2 * T3;
			      T8 = W[18];
			      Tb = W[19];
			      T1D = FNMS(T5, T3, T1C);
			      T7 = FMA(T5, T6, T4);
			      T1E = T8 * Tc;
			      Ta = T8 * T9;
			 }
			 {
			      E T1g, T1F, Td, T1f, T3i, Te, T2H;
			      T1g = ri[WS(ios, 9)];
			      T1j = ii[WS(ios, 9)];
			      T1F = FNMS(Tb, T9, T1E);
			      Td = FMA(Tb, Tc, Ta);
			      T1f = W[16];
			      T1i = W[17];
			      T1G = T1D - T1F;
			      T3i = T1D + T1F;
			      T3u = Td - T7;
			      Te = T7 + Td;
			      T2H = T1f * T1j;
			      T1h = T1f * T1g;
			      T3k = T3i + T3j;
			      T3t = FNMS(KP500000000, T3i, T3j);
			      T1B = FNMS(KP500000000, Te, T1);
			      Tf = T1 + Te;
			      T2I = FNMS(T1i, T1g, T2H);
			 }
			 T1k = FMA(T1i, T1j, T1h);
			 {
			      E T1s, T1v, T1r, T1u, T27, T1t, T1l;
			      T1s = ri[WS(ios, 4)];
			      T1v = ii[WS(ios, 4)];
			      T1r = W[6];
			      T1u = W[7];
			      T1m = ri[WS(ios, 14)];
			      T1p = ii[WS(ios, 14)];
			      T27 = T1r * T1v;
			      T1t = T1r * T1s;
			      T1l = W[26];
			      T1o = W[27];
			      T28 = FNMS(T1u, T1s, T27);
			      T1w = FMA(T1u, T1v, T1t);
			      T25 = T1l * T1p;
			      T1n = T1l * T1m;
			 }
		    }
		    {
			 E Tl, T2p, Tn, Tq, Tp, T1M, Tx, T1J, To;
			 {
			      E Th, Tk, T26, T1q, Tg, Tj;
			      Th = ri[WS(ios, 3)];
			      Tk = ii[WS(ios, 3)];
			      T26 = FNMS(T1o, T1m, T25);
			      T1q = FMA(T1o, T1p, T1n);
			      Tg = W[4];
			      Tj = W[5];
			      {
				   E T29, T2J, T1x, T2L;
				   T29 = T26 - T28;
				   T2J = T26 + T28;
				   T1x = T1q + T1w;
				   T2L = T1w - T1q;
				   {
					E T2o, Ti, T2K, T24;
					T2o = Tg * Tk;
					Ti = Tg * Th;
					T2K = FNMS(KP500000000, T2J, T2I);
					T37 = T2I + T2J;
					T24 = FNMS(KP500000000, T1x, T1k);
					T1y = T1k + T1x;
					Tl = FMA(Tj, Tk, Ti);
					T2V = FNMS(KP866025403, T2L, T2K);
					T2M = FMA(KP866025403, T2L, T2K);
					T2a = FNMS(KP866025403, T29, T24);
					T2i = FMA(KP866025403, T29, T24);
					T2p = FNMS(Tj, Th, T2o);
				   }
			      }
			 }
			 {
			      E Tt, Tw, Ts, Tv, T1L, Tu, Tm;
			      Tt = ri[WS(ios, 13)];
			      Tw = ii[WS(ios, 13)];
			      Ts = W[24];
			      Tv = W[25];
			      Tn = ri[WS(ios, 8)];
			      Tq = ii[WS(ios, 8)];
			      T1L = Ts * Tw;
			      Tu = Ts * Tt;
			      Tm = W[14];
			      Tp = W[15];
			      T1M = FNMS(Tv, Tt, T1L);
			      Tx = FMA(Tv, Tw, Tu);
			      T1J = Tm * Tq;
			      To = Tm * Tn;
			 }
			 {
			      E TF, T2v, TH, TK, TJ, T1T, TR, T1Q, TI;
			      {
				   E TB, TE, T1K, Tr, TA, TD;
				   TB = ri[WS(ios, 12)];
				   TE = ii[WS(ios, 12)];
				   T1K = FNMS(Tp, Tn, T1J);
				   Tr = FMA(Tp, Tq, To);
				   TA = W[22];
				   TD = W[23];
				   {
					E T1N, T2q, Ty, T2s;
					T1N = T1K - T1M;
					T2q = T1K + T1M;
					Ty = Tr + Tx;
					T2s = Tx - Tr;
					{
					     E T2u, TC, T2r, T1I;
					     T2u = TA * TE;
					     TC = TA * TB;
					     T2r = FNMS(KP500000000, T2q, T2p);
					     T39 = T2p + T2q;
					     T1I = FNMS(KP500000000, Ty, Tl);
					     Tz = Tl + Ty;
					     TF = FMA(TD, TE, TC);
					     T2X = FNMS(KP866025403, T2s, T2r);
					     T2t = FMA(KP866025403, T2s, T2r);
					     T1O = FNMS(KP866025403, T1N, T1I);
					     T2e = FMA(KP866025403, T1N, T1I);
					     T2v = FNMS(TD, TB, T2u);
					}
				   }
			      }
			      {
				   E TN, TQ, TM, TP, T1S, TO, TG;
				   TN = ri[WS(ios, 7)];
				   TQ = ii[WS(ios, 7)];
				   TM = W[12];
				   TP = W[13];
				   TH = ri[WS(ios, 2)];
				   TK = ii[WS(ios, 2)];
				   T1S = TM * TQ;
				   TO = TM * TN;
				   TG = W[2];
				   TJ = W[3];
				   T1T = FNMS(TP, TN, T1S);
				   TR = FMA(TP, TQ, TO);
				   T1Q = TG * TK;
				   TI = TG * TH;
			      }
			      {
				   E TW, TZ, T1R, TL, TV, TY;
				   TW = ri[WS(ios, 6)];
				   TZ = ii[WS(ios, 6)];
				   T1R = FNMS(TJ, TH, T1Q);
				   TL = FMA(TJ, TK, TI);
				   TV = W[10];
				   TY = W[11];
				   {
					E T1U, T2w, TS, T2y;
					T1U = T1R - T1T;
					T2w = T1R + T1T;
					TS = TL + TR;
					T2y = TR - TL;
					{
					     E T2B, TX, T2x, T1P;
					     T2B = TV * TZ;
					     TX = TV * TW;
					     T2x = FNMS(KP500000000, T2w, T2v);
					     T3a = T2v + T2w;
					     T1P = FNMS(KP500000000, TS, TF);
					     TT = TF + TS;
					     T10 = FMA(TY, TZ, TX);
					     T2Y = FNMS(KP866025403, T2y, T2x);
					     T2z = FMA(KP866025403, T2y, T2x);
					     T1V = FNMS(KP866025403, T1U, T1P);
					     T2f = FMA(KP866025403, T1U, T1P);
					     T2C = FNMS(TY, TW, T2B);
					}
				   }
			      }
			      {
				   E T18, T1b, T17, T1a, T20, T19, T11;
				   T18 = ri[WS(ios, 1)];
				   T1b = ii[WS(ios, 1)];
				   T17 = W[0];
				   T1a = W[1];
				   T12 = ri[WS(ios, 11)];
				   T15 = ii[WS(ios, 11)];
				   T20 = T17 * T1b;
				   T19 = T17 * T18;
				   T11 = W[20];
				   T14 = W[21];
				   T21 = FNMS(T1a, T18, T20);
				   T1c = FMA(T1a, T1b, T19);
				   T1Y = T11 * T15;
				   T13 = T11 * T12;
			      }
			 }
		    }
	       }
	       {
		    E T2G, T2h, T3J, T3I, T32, T30, T1H, T1W, T3P, T3O, T2b;
		    {
			 E T3f, T3b, T1Z, T16, T3p, TU;
			 T3f = T39 + T3a;
			 T3b = T39 - T3a;
			 T1Z = FNMS(T14, T12, T1Y);
			 T16 = FMA(T14, T15, T13);
			 T3p = Tz - TT;
			 TU = Tz + TT;
			 {
			      E T3g, T2U, T23, T3c, T3e, T3q, T3s, T1A, T34, T3r, T3n;
			      {
				   E T22, T1d, T2F, T2E, T36, T2D;
				   T22 = T1Z - T21;
				   T2D = T1Z + T21;
				   T1d = T16 + T1c;
				   T2F = T1c - T16;
				   T2E = FNMS(KP500000000, T2D, T2C);
				   T36 = T2C + T2D;
				   {
					E T1e, T1X, T38, T1z, T3o;
					T1e = T10 + T1d;
					T1X = FNMS(KP500000000, T1d, T10);
					T38 = T36 - T37;
					T3g = T36 + T37;
					T2G = FMA(KP866025403, T2F, T2E);
					T2U = FNMS(KP866025403, T2F, T2E);
					T1z = T1e + T1y;
					T3o = T1e - T1y;
					T2h = FMA(KP866025403, T22, T1X);
					T23 = FNMS(KP866025403, T22, T1X);
					T3c = FNMS(KP618033988, T3b, T38);
					T3e = FMA(KP618033988, T38, T3b);
					T3q = FNMS(KP618033988, T3p, T3o);
					T3s = FMA(KP618033988, T3o, T3p);
					T1A = TU + T1z;
					T34 = TU - T1z;
				   }
			      }
			      {
				   E T2W, T33, T3m, T3h, T2Z, T3d, T35, T3l;
				   T3J = T2U + T2V;
				   T2W = T2U - T2V;
				   ri[0] = Tf + T1A;
				   T33 = FNMS(KP250000000, T1A, Tf);
				   T3m = T3f - T3g;
				   T3h = T3f + T3g;
				   T2Z = T2X - T2Y;
				   T3I = T2X + T2Y;
				   T3d = FMA(KP559016994, T34, T33);
				   T35 = FNMS(KP559016994, T34, T33);
				   ii[0] = T3h + T3k;
				   T3l = FNMS(KP250000000, T3h, T3k);
				   ri[WS(ios, 3)] = FMA(KP951056516, T3c, T35);
				   ri[WS(ios, 12)] = FNMS(KP951056516, T3c, T35);
				   ri[WS(ios, 6)] = FMA(KP951056516, T3e, T3d);
				   ri[WS(ios, 9)] = FNMS(KP951056516, T3e, T3d);
				   T3r = FMA(KP559016994, T3m, T3l);
				   T3n = FNMS(KP559016994, T3m, T3l);
				   T32 = FMA(KP618033988, T2W, T2Z);
				   T30 = FNMS(KP618033988, T2Z, T2W);
			      }
			      ii[WS(ios, 12)] = FMA(KP951056516, T3q, T3n);
			      ii[WS(ios, 3)] = FNMS(KP951056516, T3q, T3n);
			      ii[WS(ios, 9)] = FMA(KP951056516, T3s, T3r);
			      ii[WS(ios, 6)] = FNMS(KP951056516, T3s, T3r);
			      T2d = FMA(KP866025403, T1G, T1B);
			      T1H = FNMS(KP866025403, T1G, T1B);
			      T1W = T1O + T1V;
			      T3P = T1O - T1V;
			      T3O = T23 - T2a;
			      T2b = T23 + T2a;
			 }
		    }
		    {
			 E T3H, T3v, T2S, T3Q, T3S, T2R, T2c;
			 T3H = FNMS(KP866025403, T3u, T3t);
			 T3v = FMA(KP866025403, T3u, T3t);
			 T2c = T1W + T2b;
			 T2S = T1W - T2b;
			 T3Q = FNMS(KP618033988, T3P, T3O);
			 T3S = FMA(KP618033988, T3O, T3P);
			 ri[WS(ios, 5)] = T1H + T2c;
			 T2R = FNMS(KP250000000, T2c, T1H);
			 {
			      E T2g, T2j, T3G, T3E, T2A, T2N, T3y, T3A, T3M, T3L, T3z, T3F, T3B;
			      {
				   E T3C, T3D, T31, T2T, T3K;
				   T2g = T2e + T2f;
				   T3C = T2e - T2f;
				   T3D = T2h - T2i;
				   T2j = T2h + T2i;
				   T31 = FMA(KP559016994, T2S, T2R);
				   T2T = FNMS(KP559016994, T2S, T2R);
				   T3K = T3I + T3J;
				   T3M = T3I - T3J;
				   ri[WS(ios, 8)] = FMA(KP951056516, T30, T2T);
				   ri[WS(ios, 2)] = FNMS(KP951056516, T30, T2T);
				   ri[WS(ios, 11)] = FMA(KP951056516, T32, T31);
				   ri[WS(ios, 14)] = FNMS(KP951056516, T32, T31);
				   ii[WS(ios, 5)] = T3K + T3H;
				   T3L = FNMS(KP250000000, T3K, T3H);
				   T3G = FNMS(KP618033988, T3C, T3D);
				   T3E = FMA(KP618033988, T3D, T3C);
			      }
			      {
				   E T3N, T3R, T3w, T3x;
				   T3N = FNMS(KP559016994, T3M, T3L);
				   T3R = FMA(KP559016994, T3M, T3L);
				   T3w = T2t + T2z;
				   T2A = T2t - T2z;
				   T2N = T2G - T2M;
				   T3x = T2G + T2M;
				   ii[WS(ios, 8)] = FNMS(KP951056516, T3Q, T3N);
				   ii[WS(ios, 2)] = FMA(KP951056516, T3Q, T3N);
				   ii[WS(ios, 14)] = FMA(KP951056516, T3S, T3R);
				   ii[WS(ios, 11)] = FNMS(KP951056516, T3S, T3R);
				   T3y = T3w + T3x;
				   T3A = T3w - T3x;
			      }
			      ii[WS(ios, 10)] = T3y + T3v;
			      T3z = FNMS(KP250000000, T3y, T3v);
			      T2O = FMA(KP618033988, T2N, T2A);
			      T2Q = FNMS(KP618033988, T2A, T2N);
			      T3F = FNMS(KP559016994, T3A, T3z);
			      T3B = FMA(KP559016994, T3A, T3z);
			      ii[WS(ios, 4)] = FMA(KP951056516, T3E, T3B);
			      ii[WS(ios, 1)] = FNMS(KP951056516, T3E, T3B);
			      ii[WS(ios, 13)] = FNMS(KP951056516, T3G, T3F);
			      ii[WS(ios, 7)] = FMA(KP951056516, T3G, T3F);
			      T2m = T2g - T2j;
			      T2k = T2g + T2j;
			 }
		    }
	       }
	  }
	  ri[WS(ios, 10)] = T2d + T2k;
	  T2l = FNMS(KP250000000, T2k, T2d);
	  T2P = FNMS(KP559016994, T2m, T2l);
	  T2n = FMA(KP559016994, T2m, T2l);
	  ri[WS(ios, 1)] = FMA(KP951056516, T2O, T2n);
	  ri[WS(ios, 4)] = FNMS(KP951056516, T2O, T2n);
	  ri[WS(ios, 13)] = FMA(KP951056516, T2Q, T2P);
	  ri[WS(ios, 7)] = FNMS(KP951056516, T2Q, T2P);
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 15},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 15, "t1_15", twinstr, &GENUS, {72, 28, 112, 0}, 0, 0, 0 };

void X(codelet_t1_15) (planner *p) {
     X(kdft_dit_register) (p, t1_15, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle -compact -variables 4 -pipeline-latency 4 -n 15 -name t1_15 -include t.h */

/*
 * This function contains 184 FP additions, 112 FP multiplications,
 * (or, 128 additions, 56 multiplications, 56 fused multiply/add),
 * 65 stack variables, and 60 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle.ml,v 1.23 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t.h"

static const R *t1_15(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     DK(KP587785252, +0.587785252292473129168705954639072768597652438);
     DK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     INT i;
     for (i = m; i > 0; i = i - 1, ri = ri + dist, ii = ii + dist, W = W + 28, MAKE_VOLATILE_STRIDE(ios)) {
	  E T1q, T34, Td, T1n, T2S, T35, T13, T1k, T1l, T2E, T2F, T2O, T1H, T1T, T2k;
	  E T2t, T2f, T2s, T1M, T1U, Tu, TL, TM, T2H, T2I, T2N, T1w, T1Q, T29, T2w;
	  E T24, T2v, T1B, T1R;
	  {
	       E T1, T2R, T6, T1o, Tb, T1p, Tc, T2Q;
	       T1 = ri[0];
	       T2R = ii[0];
	       {
		    E T3, T5, T2, T4;
		    T3 = ri[WS(ios, 5)];
		    T5 = ii[WS(ios, 5)];
		    T2 = W[8];
		    T4 = W[9];
		    T6 = FMA(T2, T3, T4 * T5);
		    T1o = FNMS(T4, T3, T2 * T5);
	       }
	       {
		    E T8, Ta, T7, T9;
		    T8 = ri[WS(ios, 10)];
		    Ta = ii[WS(ios, 10)];
		    T7 = W[18];
		    T9 = W[19];
		    Tb = FMA(T7, T8, T9 * Ta);
		    T1p = FNMS(T9, T8, T7 * Ta);
	       }
	       T1q = KP866025403 * (T1o - T1p);
	       T34 = KP866025403 * (Tb - T6);
	       Tc = T6 + Tb;
	       Td = T1 + Tc;
	       T1n = FNMS(KP500000000, Tc, T1);
	       T2Q = T1o + T1p;
	       T2S = T2Q + T2R;
	       T35 = FNMS(KP500000000, T2Q, T2R);
	  }
	  {
	       E TR, T2c, T18, T2h, TW, T1E, T11, T1F, T12, T2d, T1d, T1J, T1i, T1K, T1j;
	       E T2i;
	       {
		    E TO, TQ, TN, TP;
		    TO = ri[WS(ios, 6)];
		    TQ = ii[WS(ios, 6)];
		    TN = W[10];
		    TP = W[11];
		    TR = FMA(TN, TO, TP * TQ);
		    T2c = FNMS(TP, TO, TN * TQ);
	       }
	       {
		    E T15, T17, T14, T16;
		    T15 = ri[WS(ios, 9)];
		    T17 = ii[WS(ios, 9)];
		    T14 = W[16];
		    T16 = W[17];
		    T18 = FMA(T14, T15, T16 * T17);
		    T2h = FNMS(T16, T15, T14 * T17);
	       }
	       {
		    E TT, TV, TS, TU;
		    TT = ri[WS(ios, 11)];
		    TV = ii[WS(ios, 11)];
		    TS = W[20];
		    TU = W[21];
		    TW = FMA(TS, TT, TU * TV);
		    T1E = FNMS(TU, TT, TS * TV);
	       }
	       {
		    E TY, T10, TX, TZ;
		    TY = ri[WS(ios, 1)];
		    T10 = ii[WS(ios, 1)];
		    TX = W[0];
		    TZ = W[1];
		    T11 = FMA(TX, TY, TZ * T10);
		    T1F = FNMS(TZ, TY, TX * T10);
	       }
	       T12 = TW + T11;
	       T2d = T1E + T1F;
	       {
		    E T1a, T1c, T19, T1b;
		    T1a = ri[WS(ios, 14)];
		    T1c = ii[WS(ios, 14)];
		    T19 = W[26];
		    T1b = W[27];
		    T1d = FMA(T19, T1a, T1b * T1c);
		    T1J = FNMS(T1b, T1a, T19 * T1c);
	       }
	       {
		    E T1f, T1h, T1e, T1g;
		    T1f = ri[WS(ios, 4)];
		    T1h = ii[WS(ios, 4)];
		    T1e = W[6];
		    T1g = W[7];
		    T1i = FMA(T1e, T1f, T1g * T1h);
		    T1K = FNMS(T1g, T1f, T1e * T1h);
	       }
	       T1j = T1d + T1i;
	       T2i = T1J + T1K;
	       {
		    E T1D, T1G, T2g, T2j;
		    T13 = TR + T12;
		    T1k = T18 + T1j;
		    T1l = T13 + T1k;
		    T2E = T2c + T2d;
		    T2F = T2h + T2i;
		    T2O = T2E + T2F;
		    T1D = FNMS(KP500000000, T12, TR);
		    T1G = KP866025403 * (T1E - T1F);
		    T1H = T1D - T1G;
		    T1T = T1D + T1G;
		    T2g = KP866025403 * (T1i - T1d);
		    T2j = FNMS(KP500000000, T2i, T2h);
		    T2k = T2g + T2j;
		    T2t = T2j - T2g;
		    {
			 E T2b, T2e, T1I, T1L;
			 T2b = KP866025403 * (T11 - TW);
			 T2e = FNMS(KP500000000, T2d, T2c);
			 T2f = T2b + T2e;
			 T2s = T2e - T2b;
			 T1I = FNMS(KP500000000, T1j, T18);
			 T1L = KP866025403 * (T1J - T1K);
			 T1M = T1I - T1L;
			 T1U = T1I + T1L;
		    }
	       }
	  }
	  {
	       E Ti, T21, Tz, T26, Tn, T1t, Ts, T1u, Tt, T22, TE, T1y, TJ, T1z, TK;
	       E T27;
	       {
		    E Tf, Th, Te, Tg;
		    Tf = ri[WS(ios, 3)];
		    Th = ii[WS(ios, 3)];
		    Te = W[4];
		    Tg = W[5];
		    Ti = FMA(Te, Tf, Tg * Th);
		    T21 = FNMS(Tg, Tf, Te * Th);
	       }
	       {
		    E Tw, Ty, Tv, Tx;
		    Tw = ri[WS(ios, 12)];
		    Ty = ii[WS(ios, 12)];
		    Tv = W[22];
		    Tx = W[23];
		    Tz = FMA(Tv, Tw, Tx * Ty);
		    T26 = FNMS(Tx, Tw, Tv * Ty);
	       }
	       {
		    E Tk, Tm, Tj, Tl;
		    Tk = ri[WS(ios, 8)];
		    Tm = ii[WS(ios, 8)];
		    Tj = W[14];
		    Tl = W[15];
		    Tn = FMA(Tj, Tk, Tl * Tm);
		    T1t = FNMS(Tl, Tk, Tj * Tm);
	       }
	       {
		    E Tp, Tr, To, Tq;
		    Tp = ri[WS(ios, 13)];
		    Tr = ii[WS(ios, 13)];
		    To = W[24];
		    Tq = W[25];
		    Ts = FMA(To, Tp, Tq * Tr);
		    T1u = FNMS(Tq, Tp, To * Tr);
	       }
	       Tt = Tn + Ts;
	       T22 = T1t + T1u;
	       {
		    E TB, TD, TA, TC;
		    TB = ri[WS(ios, 2)];
		    TD = ii[WS(ios, 2)];
		    TA = W[2];
		    TC = W[3];
		    TE = FMA(TA, TB, TC * TD);
		    T1y = FNMS(TC, TB, TA * TD);
	       }
	       {
		    E TG, TI, TF, TH;
		    TG = ri[WS(ios, 7)];
		    TI = ii[WS(ios, 7)];
		    TF = W[12];
		    TH = W[13];
		    TJ = FMA(TF, TG, TH * TI);
		    T1z = FNMS(TH, TG, TF * TI);
	       }
	       TK = TE + TJ;
	       T27 = T1y + T1z;
	       {
		    E T1s, T1v, T25, T28;
		    Tu = Ti + Tt;
		    TL = Tz + TK;
		    TM = Tu + TL;
		    T2H = T21 + T22;
		    T2I = T26 + T27;
		    T2N = T2H + T2I;
		    T1s = FNMS(KP500000000, Tt, Ti);
		    T1v = KP866025403 * (T1t - T1u);
		    T1w = T1s - T1v;
		    T1Q = T1s + T1v;
		    T25 = KP866025403 * (TJ - TE);
		    T28 = FNMS(KP500000000, T27, T26);
		    T29 = T25 + T28;
		    T2w = T28 - T25;
		    {
			 E T20, T23, T1x, T1A;
			 T20 = KP866025403 * (Ts - Tn);
			 T23 = FNMS(KP500000000, T22, T21);
			 T24 = T20 + T23;
			 T2v = T23 - T20;
			 T1x = FNMS(KP500000000, TK, Tz);
			 T1A = KP866025403 * (T1y - T1z);
			 T1B = T1x - T1A;
			 T1R = T1x + T1A;
		    }
	       }
	  }
	  {
	       E T2C, T1m, T2B, T2K, T2M, T2G, T2J, T2L, T2D;
	       T2C = KP559016994 * (TM - T1l);
	       T1m = TM + T1l;
	       T2B = FNMS(KP250000000, T1m, Td);
	       T2G = T2E - T2F;
	       T2J = T2H - T2I;
	       T2K = FNMS(KP587785252, T2J, KP951056516 * T2G);
	       T2M = FMA(KP951056516, T2J, KP587785252 * T2G);
	       ri[0] = Td + T1m;
	       T2L = T2C + T2B;
	       ri[WS(ios, 9)] = T2L - T2M;
	       ri[WS(ios, 6)] = T2L + T2M;
	       T2D = T2B - T2C;
	       ri[WS(ios, 12)] = T2D - T2K;
	       ri[WS(ios, 3)] = T2D + T2K;
	  }
	  {
	       E T2U, T2P, T2T, T2Y, T30, T2W, T2X, T2Z, T2V;
	       T2U = KP559016994 * (T2N - T2O);
	       T2P = T2N + T2O;
	       T2T = FNMS(KP250000000, T2P, T2S);
	       T2W = T13 - T1k;
	       T2X = Tu - TL;
	       T2Y = FNMS(KP587785252, T2X, KP951056516 * T2W);
	       T30 = FMA(KP951056516, T2X, KP587785252 * T2W);
	       ii[0] = T2P + T2S;
	       T2Z = T2U + T2T;
	       ii[WS(ios, 6)] = T2Z - T30;
	       ii[WS(ios, 9)] = T30 + T2Z;
	       T2V = T2T - T2U;
	       ii[WS(ios, 3)] = T2V - T2Y;
	       ii[WS(ios, 12)] = T2Y + T2V;
	  }
	  {
	       E T2y, T2A, T1r, T1O, T2p, T2q, T2z, T2r;
	       {
		    E T2u, T2x, T1C, T1N;
		    T2u = T2s - T2t;
		    T2x = T2v - T2w;
		    T2y = FNMS(KP587785252, T2x, KP951056516 * T2u);
		    T2A = FMA(KP951056516, T2x, KP587785252 * T2u);
		    T1r = T1n - T1q;
		    T1C = T1w + T1B;
		    T1N = T1H + T1M;
		    T1O = T1C + T1N;
		    T2p = FNMS(KP250000000, T1O, T1r);
		    T2q = KP559016994 * (T1C - T1N);
	       }
	       ri[WS(ios, 5)] = T1r + T1O;
	       T2z = T2q + T2p;
	       ri[WS(ios, 14)] = T2z - T2A;
	       ri[WS(ios, 11)] = T2z + T2A;
	       T2r = T2p - T2q;
	       ri[WS(ios, 2)] = T2r - T2y;
	       ri[WS(ios, 8)] = T2r + T2y;
	  }
	  {
	       E T3h, T3q, T3i, T3l, T3m, T3n, T3p, T3o;
	       {
		    E T3f, T3g, T3j, T3k;
		    T3f = T1H - T1M;
		    T3g = T1w - T1B;
		    T3h = FNMS(KP587785252, T3g, KP951056516 * T3f);
		    T3q = FMA(KP951056516, T3g, KP587785252 * T3f);
		    T3i = T35 - T34;
		    T3j = T2v + T2w;
		    T3k = T2s + T2t;
		    T3l = T3j + T3k;
		    T3m = FNMS(KP250000000, T3l, T3i);
		    T3n = KP559016994 * (T3j - T3k);
	       }
	       ii[WS(ios, 5)] = T3l + T3i;
	       T3p = T3n + T3m;
	       ii[WS(ios, 11)] = T3p - T3q;
	       ii[WS(ios, 14)] = T3q + T3p;
	       T3o = T3m - T3n;
	       ii[WS(ios, 2)] = T3h + T3o;
	       ii[WS(ios, 8)] = T3o - T3h;
	  }
	  {
	       E T3c, T3d, T36, T37, T33, T38, T3e, T39;
	       {
		    E T3a, T3b, T31, T32;
		    T3a = T1Q - T1R;
		    T3b = T1T - T1U;
		    T3c = FMA(KP951056516, T3a, KP587785252 * T3b);
		    T3d = FNMS(KP587785252, T3a, KP951056516 * T3b);
		    T36 = T34 + T35;
		    T31 = T24 + T29;
		    T32 = T2f + T2k;
		    T37 = T31 + T32;
		    T33 = KP559016994 * (T31 - T32);
		    T38 = FNMS(KP250000000, T37, T36);
	       }
	       ii[WS(ios, 10)] = T37 + T36;
	       T3e = T38 - T33;
	       ii[WS(ios, 7)] = T3d + T3e;
	       ii[WS(ios, 13)] = T3e - T3d;
	       T39 = T33 + T38;
	       ii[WS(ios, 1)] = T39 - T3c;
	       ii[WS(ios, 4)] = T3c + T39;
	  }
	  {
	       E T2m, T2o, T1P, T1W, T1X, T1Y, T2n, T1Z;
	       {
		    E T2a, T2l, T1S, T1V;
		    T2a = T24 - T29;
		    T2l = T2f - T2k;
		    T2m = FMA(KP951056516, T2a, KP587785252 * T2l);
		    T2o = FNMS(KP587785252, T2a, KP951056516 * T2l);
		    T1P = T1n + T1q;
		    T1S = T1Q + T1R;
		    T1V = T1T + T1U;
		    T1W = T1S + T1V;
		    T1X = KP559016994 * (T1S - T1V);
		    T1Y = FNMS(KP250000000, T1W, T1P);
	       }
	       ri[WS(ios, 10)] = T1P + T1W;
	       T2n = T1Y - T1X;
	       ri[WS(ios, 7)] = T2n - T2o;
	       ri[WS(ios, 13)] = T2n + T2o;
	       T1Z = T1X + T1Y;
	       ri[WS(ios, 4)] = T1Z - T2m;
	       ri[WS(ios, 1)] = T1Z + T2m;
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 15},
     {TW_NEXT, 1, 0}
};

static const ct_desc desc = { 15, "t1_15", twinstr, &GENUS, {128, 56, 56, 0}, 0, 0, 0 };

void X(codelet_t1_15) (planner *p) {
     X(kdft_dit_register) (p, t1_15, &desc);
}
#endif				/* HAVE_FMA */
