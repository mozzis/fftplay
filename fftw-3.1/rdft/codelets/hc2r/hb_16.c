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
/* Generated on Fri Jan 27 20:42:44 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 16 -dif -name hb_16 -include hb.h */

/*
 * This function contains 174 FP additions, 100 FP multiplications,
 * (or, 104 additions, 30 multiplications, 70 fused multiply/add),
 * 83 stack variables, and 64 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hb.h"

static const R *hb_16(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 30, MAKE_VOLATILE_STRIDE(ios)) {
	  E T3v, T3s, T3u, T3w, T3t;
	  {
	       E T26, T3j, T2z, T36, T11, T1K, T18, T1L, T1C, Tf, T37, T2d, T1m, TE, T3k;
	       E T2C, T1J, Tu, T20, T1F, T1a, TN, T3n, T3e, T1b, TW, T2k, T2h, T2F, T2s;
	       E T3m, T3b;
	       {
		    E TD, Tw, T2B, T2A;
		    {
			 E T24, T3, T2y, T14, T2x, T6, T25, T17, Tb, T2b, Ta, T2a, Tz, Tc, TA;
			 E TB;
			 {
			      E T4, T5, T15, T16;
			      {
				   E T1, T2, T12, T13;
				   T1 = rio[0];
				   T2 = iio[-WS(ios, 8)];
				   T12 = iio[0];
				   T13 = rio[WS(ios, 8)];
				   T4 = rio[WS(ios, 4)];
				   T24 = T1 - T2;
				   T3 = T1 + T2;
				   T2y = T12 + T13;
				   T14 = T12 - T13;
				   T5 = iio[-WS(ios, 12)];
				   T15 = iio[-WS(ios, 4)];
				   T16 = rio[WS(ios, 12)];
			      }
			      {
				   E T8, T9, Tx, Ty;
				   T8 = rio[WS(ios, 2)];
				   T2x = T4 - T5;
				   T6 = T4 + T5;
				   T25 = T15 + T16;
				   T17 = T15 - T16;
				   T9 = iio[-WS(ios, 10)];
				   Tx = iio[-WS(ios, 2)];
				   Ty = rio[WS(ios, 10)];
				   Tb = iio[-WS(ios, 14)];
				   T2b = T8 - T9;
				   Ta = T8 + T9;
				   T2a = Tx + Ty;
				   Tz = Tx - Ty;
				   Tc = rio[WS(ios, 6)];
				   TA = iio[-WS(ios, 6)];
				   TB = rio[WS(ios, 14)];
			      }
			 }
			 {
			      E T27, T28, TC, Te, Td, T7, T29, T2c;
			      T26 = T24 - T25;
			      T3j = T24 + T25;
			      T27 = Tb - Tc;
			      Td = Tb + Tc;
			      T28 = TB + TA;
			      TC = TA - TB;
			      T2z = T2x + T2y;
			      T36 = T2y - T2x;
			      Te = Ta + Td;
			      T11 = Td - Ta;
			      T1K = T14 + T17;
			      T18 = T14 - T17;
			      TD = Tz - TC;
			      T1L = Tz + TC;
			      Tw = T3 - T6;
			      T7 = T3 + T6;
			      T2B = T27 + T28;
			      T29 = T27 - T28;
			      T2c = T2a - T2b;
			      T2A = T2b + T2a;
			      T1C = T7 - Te;
			      Tf = T7 + Te;
			      T37 = T2c + T29;
			      T2d = T29 - T2c;
			 }
		    }
		    {
			 E T2f, Ti, T2j, TI, T2i, Tl, T2g, TL, Tq, T2m, Tp, T2q, TR, Tr, TS;
			 E TT;
			 {
			      E Tj, Tk, TJ, TK;
			      {
				   E Tg, Th, TG, TH;
				   Tg = rio[WS(ios, 1)];
				   T1m = Tw - TD;
				   TE = Tw + TD;
				   T3k = T2A + T2B;
				   T2C = T2A - T2B;
				   Th = iio[-WS(ios, 9)];
				   TG = iio[-WS(ios, 1)];
				   TH = rio[WS(ios, 9)];
				   Tj = rio[WS(ios, 5)];
				   T2f = Tg - Th;
				   Ti = Tg + Th;
				   T2j = TG + TH;
				   TI = TG - TH;
				   Tk = iio[-WS(ios, 13)];
				   TJ = iio[-WS(ios, 5)];
				   TK = rio[WS(ios, 13)];
			      }
			      {
				   E Tn, To, TP, TQ;
				   Tn = iio[-WS(ios, 15)];
				   T2i = Tj - Tk;
				   Tl = Tj + Tk;
				   T2g = TJ + TK;
				   TL = TJ - TK;
				   To = rio[WS(ios, 7)];
				   TP = iio[-WS(ios, 7)];
				   TQ = rio[WS(ios, 15)];
				   Tq = rio[WS(ios, 3)];
				   T2m = Tn - To;
				   Tp = Tn + To;
				   T2q = TQ + TP;
				   TR = TP - TQ;
				   Tr = iio[-WS(ios, 11)];
				   TS = iio[-WS(ios, 3)];
				   TT = rio[WS(ios, 11)];
			      }
			 }
			 {
			      E TO, TV, T3c, T2r, T3d, T2o, T39, T3a;
			      {
				   E TF, Tm, T2p, T2n, TU, T1D, TM, Tt, Ts, T1E;
				   TF = Ti - Tl;
				   Tm = Ti + Tl;
				   T2p = Tq - Tr;
				   Ts = Tq + Tr;
				   T2n = TS + TT;
				   TU = TS - TT;
				   T1D = TI + TL;
				   TM = TI - TL;
				   TO = Tp - Ts;
				   Tt = Tp + Ts;
				   TV = TR - TU;
				   T1E = TR + TU;
				   T1J = Tt - Tm;
				   Tu = Tm + Tt;
				   T3c = T2p + T2q;
				   T2r = T2p - T2q;
				   T20 = T1D + T1E;
				   T1F = T1D - T1E;
				   T1a = TM - TF;
				   TN = TF + TM;
				   T3d = T2m + T2n;
				   T2o = T2m - T2n;
			      }
			      T3n = FMA(KP414213562, T3c, T3d);
			      T3e = FNMS(KP414213562, T3d, T3c);
			      T1b = TO + TV;
			      TW = TO - TV;
			      T2k = T2i + T2j;
			      T39 = T2j - T2i;
			      T3a = T2f + T2g;
			      T2h = T2f - T2g;
			      T2F = FNMS(KP414213562, T2o, T2r);
			      T2s = FMA(KP414213562, T2r, T2o);
			      T3m = FMA(KP414213562, T39, T3a);
			      T3b = FNMS(KP414213562, T3a, T39);
			 }
		    }
	       }
	       {
		    E T2E, T2l, T1c, T19, TX, T1z, T1v, T1y, T1x, T1A;
		    {
			 E T1M, T1W, T21, T1V, T1Y, T1Z;
			 rio[0] = Tf + Tu;
			 T1Z = T1L + T1K;
			 T1M = T1K - T1L;
			 T2E = FMA(KP414213562, T2h, T2k);
			 T2l = FNMS(KP414213562, T2k, T2h);
			 T1W = Tf - Tu;
			 T21 = T1Z - T20;
			 T1V = W[14];
			 T1Y = W[15];
			 iio[-WS(ios, 15)] = T20 + T1Z;
			 {
			      E T1G, T1T, T1N, T1P, T1B, T1U, T1I, T1H, T1O;
			      {
				   E T1S, T1R, T1X, T22, T1Q;
				   T1X = T1V * T1W;
				   T22 = T1Y * T1W;
				   T1G = T1C + T1F;
				   T1Q = T1C - T1F;
				   rio[WS(ios, 8)] = FNMS(T1Y, T21, T1X);
				   iio[-WS(ios, 7)] = FMA(T1V, T21, T22);
				   T1T = T1M - T1J;
				   T1N = T1J + T1M;
				   T1P = W[6];
				   T1S = W[7];
				   T1B = W[22];
				   T1R = T1P * T1Q;
				   T1U = T1S * T1Q;
				   T1I = W[23];
				   T1H = T1B * T1G;
				   rio[WS(ios, 4)] = FNMS(T1S, T1T, T1R);
			      }
			      iio[-WS(ios, 11)] = FMA(T1P, T1T, T1U);
			      T1O = T1I * T1G;
			      rio[WS(ios, 12)] = FNMS(T1I, T1N, T1H);
			      {
				   E T1r, T1s, T1w, T1o, T1n;
				   T1n = T1b - T1a;
				   T1c = T1a + T1b;
				   T19 = T11 + T18;
				   T1r = T18 - T11;
				   iio[-WS(ios, 3)] = FMA(T1B, T1N, T1O);
				   TX = TN + TW;
				   T1s = TN - TW;
				   T1w = FNMS(KP707106781, T1n, T1m);
				   T1o = FMA(KP707106781, T1n, T1m);
				   {
					E T1l, T1t, T1q, T1p, T1u;
					T1l = W[2];
					T1t = FMA(KP707106781, T1s, T1r);
					T1z = FNMS(KP707106781, T1s, T1r);
					T1q = W[3];
					T1p = T1l * T1o;
					T1v = W[18];
					T1y = W[19];
					T1u = T1q * T1o;
					rio[WS(ios, 2)] = FNMS(T1q, T1t, T1p);
					T1x = T1v * T1w;
					T1A = T1y * T1w;
					iio[-WS(ios, 13)] = FMA(T1l, T1t, T1u);
				   }
			      }
			 }
		    }
		    {
			 E T2V, T2R, T2Q, T2W, T2N, T2M, T2L;
			 {
			      E T1g, T1f, T1j, T1h, T1i, TY;
			      rio[WS(ios, 10)] = FNMS(T1y, T1z, T1x);
			      iio[-WS(ios, 5)] = FMA(T1v, T1z, T1A);
			      T1g = FNMS(KP707106781, TX, TE);
			      TY = FMA(KP707106781, TX, TE);
			      {
				   E Tv, T10, T1d, TZ, T1e;
				   Tv = W[26];
				   T10 = W[27];
				   T1f = W[10];
				   T1j = FNMS(KP707106781, T1c, T19);
				   T1d = FMA(KP707106781, T1c, T19);
				   TZ = Tv * TY;
				   T1e = T10 * TY;
				   T1h = T1f * T1g;
				   T1i = W[11];
				   rio[WS(ios, 14)] = FNMS(T10, T1d, TZ);
				   iio[-WS(ios, 1)] = FMA(Tv, T1d, T1e);
			      }
			      {
				   E T2u, T2K, T2H, T23, T2w;
				   {
					E T2e, T1k, T2t, T2D, T2G;
					T2e = FMA(KP707106781, T2d, T26);
					T2V = FNMS(KP707106781, T2d, T26);
					rio[WS(ios, 6)] = FNMS(T1i, T1j, T1h);
					T1k = T1i * T1g;
					T2R = T2s - T2l;
					T2t = T2l + T2s;
					T2D = FMA(KP707106781, T2C, T2z);
					T2Q = FNMS(KP707106781, T2C, T2z);
					T2W = T2E - T2F;
					T2G = T2E + T2F;
					iio[-WS(ios, 9)] = FMA(T1f, T1j, T1k);
					T2u = FMA(KP923879532, T2t, T2e);
					T2N = FNMS(KP923879532, T2t, T2e);
					T2K = FNMS(KP923879532, T2G, T2D);
					T2H = FMA(KP923879532, T2G, T2D);
				   }
				   T23 = W[0];
				   T2w = W[1];
				   {
					E T2J, T2I, T2v, T2O;
					T2J = W[16];
					T2M = W[17];
					T2I = T23 * T2H;
					T2v = T23 * T2u;
					T2O = T2J * T2N;
					T2L = T2J * T2K;
					iio[-WS(ios, 14)] = FMA(T2w, T2u, T2I);
					rio[WS(ios, 1)] = FNMS(T2w, T2H, T2v);
					rio[WS(ios, 9)] = FNMS(T2M, T2K, T2O);
				   }
			      }
			 }
			 iio[-WS(ios, 6)] = FMA(T2M, T2N, T2L);
			 {
			      E T33, T30, T32, T34, T31;
			      {
				   E T2P, T2S, T2X, T2U, T2T, T2Z, T2Y;
				   T2P = W[24];
				   T33 = FNMS(KP923879532, T2R, T2Q);
				   T2S = FMA(KP923879532, T2R, T2Q);
				   T30 = FNMS(KP923879532, T2W, T2V);
				   T2X = FMA(KP923879532, T2W, T2V);
				   T2U = W[25];
				   T2T = T2P * T2S;
				   T2Z = W[8];
				   T2Y = T2P * T2X;
				   T32 = W[9];
				   iio[-WS(ios, 2)] = FMA(T2U, T2X, T2T);
				   T34 = T2Z * T33;
				   T31 = T2Z * T30;
				   rio[WS(ios, 13)] = FNMS(T2U, T2S, T2Y);
			      }
			      {
				   E T3l, T3f, T38, T3o, T3L, T3I, T3K, T3M, T3J;
				   {
					E T3y, T3z, T3D, T3E;
					T3l = FMA(KP707106781, T3k, T3j);
					T3y = FNMS(KP707106781, T3k, T3j);
					iio[-WS(ios, 10)] = FMA(T32, T30, T34);
					rio[WS(ios, 5)] = FNMS(T32, T33, T31);
					T3z = T3b + T3e;
					T3f = T3b - T3e;
					T38 = FMA(KP707106781, T37, T36);
					T3D = FNMS(KP707106781, T37, T36);
					T3E = T3m - T3n;
					T3o = T3m + T3n;
					{
					     E T3x, T3A, T3F, T3C, T3B, T3H, T3G;
					     T3x = W[4];
					     T3L = FMA(KP923879532, T3z, T3y);
					     T3A = FNMS(KP923879532, T3z, T3y);
					     T3I = FNMS(KP923879532, T3E, T3D);
					     T3F = FMA(KP923879532, T3E, T3D);
					     T3C = W[5];
					     T3B = T3x * T3A;
					     T3H = W[20];
					     T3G = T3x * T3F;
					     T3K = W[21];
					     rio[WS(ios, 3)] = FNMS(T3C, T3F, T3B);
					     T3M = T3H * T3L;
					     T3J = T3H * T3I;
					     iio[-WS(ios, 12)] = FMA(T3C, T3A, T3G);
					}
				   }
				   rio[WS(ios, 11)] = FNMS(T3K, T3I, T3M);
				   iio[-WS(ios, 4)] = FMA(T3K, T3L, T3J);
				   {
					E T35, T3g, T3p, T3i, T3h, T3r, T3q;
					T35 = W[28];
					T3v = FNMS(KP923879532, T3f, T38);
					T3g = FMA(KP923879532, T3f, T38);
					T3s = FNMS(KP923879532, T3o, T3l);
					T3p = FMA(KP923879532, T3o, T3l);
					T3i = W[29];
					T3h = T35 * T3g;
					T3r = W[12];
					T3q = T35 * T3p;
					T3u = W[13];
					iio[0] = FMA(T3i, T3p, T3h);
					T3w = T3r * T3v;
					T3t = T3r * T3s;
					rio[WS(ios, 15)] = FNMS(T3i, T3g, T3q);
				   }
			      }
			 }
		    }
	       }
	  }
	  iio[-WS(ios, 8)] = FMA(T3u, T3s, T3w);
	  rio[WS(ios, 7)] = FNMS(T3u, T3v, T3t);
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 16},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 16, "hb_16", twinstr, &GENUS, {104, 30, 70, 0}, 0, 0, 0 };

void X(codelet_hb_16) (planner *p) {
     X(khc2hc_register) (p, hb_16, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc -compact -variables 4 -pipeline-latency 4 -sign 1 -n 16 -dif -name hb_16 -include hb.h */

/*
 * This function contains 174 FP additions, 84 FP multiplications,
 * (or, 136 additions, 46 multiplications, 38 fused multiply/add),
 * 50 stack variables, and 64 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hb.h"

static const R *hb_16(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 30, MAKE_VOLATILE_STRIDE(ios)) {
	  E T7, T2K, T30, Tw, T1a, T2e, T2k, T1B, Te, TD, T1C, T13, T2n, T2Z, T2b;
	  E T2L, Tm, T1v, TN, T10, T1W, T2p, T2P, T2W, Tt, T1w, TW, T11, T23, T2q;
	  E T2S, T2X;
	  {
	       E T3, T2c, T16, T2j, T6, T2i, T19, T2d;
	       {
		    E T1, T2, T14, T15;
		    T1 = rio[0];
		    T2 = iio[-WS(ios, 8)];
		    T3 = T1 + T2;
		    T2c = T1 - T2;
		    T14 = iio[0];
		    T15 = rio[WS(ios, 8)];
		    T16 = T14 - T15;
		    T2j = T14 + T15;
	       }
	       {
		    E T4, T5, T17, T18;
		    T4 = rio[WS(ios, 4)];
		    T5 = iio[-WS(ios, 12)];
		    T6 = T4 + T5;
		    T2i = T4 - T5;
		    T17 = iio[-WS(ios, 4)];
		    T18 = rio[WS(ios, 12)];
		    T19 = T17 - T18;
		    T2d = T17 + T18;
	       }
	       T7 = T3 + T6;
	       T2K = T2j - T2i;
	       T30 = T2c + T2d;
	       Tw = T3 - T6;
	       T1a = T16 - T19;
	       T2e = T2c - T2d;
	       T2k = T2i + T2j;
	       T1B = T16 + T19;
	  }
	  {
	       E Ta, T29, Tz, T28, Td, T25, TC, T26;
	       {
		    E T8, T9, Tx, Ty;
		    T8 = rio[WS(ios, 2)];
		    T9 = iio[-WS(ios, 10)];
		    Ta = T8 + T9;
		    T29 = T8 - T9;
		    Tx = iio[-WS(ios, 2)];
		    Ty = rio[WS(ios, 10)];
		    Tz = Tx - Ty;
		    T28 = Tx + Ty;
	       }
	       {
		    E Tb, Tc, TA, TB;
		    Tb = iio[-WS(ios, 14)];
		    Tc = rio[WS(ios, 6)];
		    Td = Tb + Tc;
		    T25 = Tb - Tc;
		    TA = iio[-WS(ios, 6)];
		    TB = rio[WS(ios, 14)];
		    TC = TA - TB;
		    T26 = TB + TA;
	       }
	       Te = Ta + Td;
	       TD = Tz - TC;
	       T1C = Tz + TC;
	       T13 = Td - Ta;
	       {
		    E T2l, T2m, T27, T2a;
		    T2l = T29 + T28;
		    T2m = T25 + T26;
		    T2n = KP707106781 * (T2l - T2m);
		    T2Z = KP707106781 * (T2l + T2m);
		    T27 = T25 - T26;
		    T2a = T28 - T29;
		    T2b = KP707106781 * (T27 - T2a);
		    T2L = KP707106781 * (T2a + T27);
	       }
	  }
	  {
	       E Ti, T1Q, TI, T1U, Tl, T1T, TL, T1R, TF, TM;
	       {
		    E Tg, Th, TG, TH;
		    Tg = rio[WS(ios, 1)];
		    Th = iio[-WS(ios, 9)];
		    Ti = Tg + Th;
		    T1Q = Tg - Th;
		    TG = iio[-WS(ios, 1)];
		    TH = rio[WS(ios, 9)];
		    TI = TG - TH;
		    T1U = TG + TH;
	       }
	       {
		    E Tj, Tk, TJ, TK;
		    Tj = rio[WS(ios, 5)];
		    Tk = iio[-WS(ios, 13)];
		    Tl = Tj + Tk;
		    T1T = Tj - Tk;
		    TJ = iio[-WS(ios, 5)];
		    TK = rio[WS(ios, 13)];
		    TL = TJ - TK;
		    T1R = TJ + TK;
	       }
	       Tm = Ti + Tl;
	       T1v = TI + TL;
	       TF = Ti - Tl;
	       TM = TI - TL;
	       TN = TF + TM;
	       T10 = TM - TF;
	       {
		    E T1S, T1V, T2N, T2O;
		    T1S = T1Q - T1R;
		    T1V = T1T + T1U;
		    T1W = FNMS(KP382683432, T1V, KP923879532 * T1S);
		    T2p = FMA(KP923879532, T1V, KP382683432 * T1S);
		    T2N = T1U - T1T;
		    T2O = T1Q + T1R;
		    T2P = FNMS(KP382683432, T2O, KP923879532 * T2N);
		    T2W = FMA(KP382683432, T2N, KP923879532 * T2O);
	       }
	  }
	  {
	       E Tp, T1X, TR, T21, Ts, T20, TU, T1Y, TO, TV;
	       {
		    E Tn, To, TP, TQ;
		    Tn = iio[-WS(ios, 15)];
		    To = rio[WS(ios, 7)];
		    Tp = Tn + To;
		    T1X = Tn - To;
		    TP = iio[-WS(ios, 7)];
		    TQ = rio[WS(ios, 15)];
		    TR = TP - TQ;
		    T21 = TQ + TP;
	       }
	       {
		    E Tq, Tr, TS, TT;
		    Tq = rio[WS(ios, 3)];
		    Tr = iio[-WS(ios, 11)];
		    Ts = Tq + Tr;
		    T20 = Tq - Tr;
		    TS = iio[-WS(ios, 3)];
		    TT = rio[WS(ios, 11)];
		    TU = TS - TT;
		    T1Y = TS + TT;
	       }
	       Tt = Tp + Ts;
	       T1w = TU + TR;
	       TO = Tp - Ts;
	       TV = TR - TU;
	       TW = TO - TV;
	       T11 = TO + TV;
	       {
		    E T1Z, T22, T2Q, T2R;
		    T1Z = T1X - T1Y;
		    T22 = T20 - T21;
		    T23 = FMA(KP923879532, T1Z, KP382683432 * T22);
		    T2q = FNMS(KP382683432, T1Z, KP923879532 * T22);
		    T2Q = T1X + T1Y;
		    T2R = T20 + T21;
		    T2S = FNMS(KP923879532, T2R, KP382683432 * T2Q);
		    T2X = FMA(KP923879532, T2Q, KP382683432 * T2R);
	       }
	  }
	  {
	       E Tf, Tu, T1K, T1M, T1N, T1O, T1J, T1L;
	       Tf = T7 + Te;
	       Tu = Tm + Tt;
	       T1K = Tf - Tu;
	       T1M = T1C + T1B;
	       T1N = T1v + T1w;
	       T1O = T1M - T1N;
	       rio[0] = Tf + Tu;
	       iio[-WS(ios, 15)] = T1N + T1M;
	       T1J = W[14];
	       T1L = W[15];
	       rio[WS(ios, 8)] = FNMS(T1L, T1O, T1J * T1K);
	       iio[-WS(ios, 7)] = FMA(T1L, T1K, T1J * T1O);
	  }
	  {
	       E T2U, T36, T32, T34;
	       {
		    E T2M, T2T, T2Y, T31;
		    T2M = T2K + T2L;
		    T2T = T2P + T2S;
		    T2U = T2M + T2T;
		    T36 = T2M - T2T;
		    T2Y = T2W + T2X;
		    T31 = T2Z + T30;
		    T32 = T2Y + T31;
		    T34 = T31 - T2Y;
	       }
	       {
		    E T2J, T2V, T33, T35;
		    T2J = W[28];
		    T2V = W[29];
		    iio[0] = FMA(T2J, T2U, T2V * T32);
		    rio[WS(ios, 15)] = FNMS(T2V, T2U, T2J * T32);
		    T33 = W[12];
		    T35 = W[13];
		    rio[WS(ios, 7)] = FNMS(T35, T36, T33 * T34);
		    iio[-WS(ios, 8)] = FMA(T33, T36, T35 * T34);
	       }
	  }
	  {
	       E TY, T1e, T1c, T1g;
	       {
		    E TE, TX, T12, T1b;
		    TE = Tw + TD;
		    TX = KP707106781 * (TN + TW);
		    TY = TE + TX;
		    T1e = TE - TX;
		    T12 = KP707106781 * (T10 + T11);
		    T1b = T13 + T1a;
		    T1c = T12 + T1b;
		    T1g = T1b - T12;
	       }
	       {
		    E Tv, TZ, T1d, T1f;
		    Tv = W[26];
		    TZ = W[27];
		    rio[WS(ios, 14)] = FNMS(TZ, T1c, Tv * TY);
		    iio[-WS(ios, 1)] = FMA(TZ, TY, Tv * T1c);
		    T1d = W[10];
		    T1f = W[11];
		    rio[WS(ios, 6)] = FNMS(T1f, T1g, T1d * T1e);
		    iio[-WS(ios, 9)] = FMA(T1f, T1e, T1d * T1g);
	       }
	  }
	  {
	       E T2g, T2w, T2s, T2u;
	       {
		    E T24, T2f, T2o, T2r;
		    T24 = T1W + T23;
		    T2f = T2b + T2e;
		    T2g = T24 + T2f;
		    T2w = T2f - T24;
		    T2o = T2k + T2n;
		    T2r = T2p + T2q;
		    T2s = T2o + T2r;
		    T2u = T2o - T2r;
	       }
	       {
		    E T1P, T2h, T2t, T2v;
		    T1P = W[0];
		    T2h = W[1];
		    rio[WS(ios, 1)] = FNMS(T2h, T2s, T1P * T2g);
		    iio[-WS(ios, 14)] = FMA(T1P, T2s, T2h * T2g);
		    T2t = W[16];
		    T2v = W[17];
		    iio[-WS(ios, 6)] = FMA(T2t, T2u, T2v * T2w);
		    rio[WS(ios, 9)] = FNMS(T2v, T2u, T2t * T2w);
	       }
	  }
	  {
	       E T1k, T1q, T1o, T1s;
	       {
		    E T1i, T1j, T1m, T1n;
		    T1i = Tw - TD;
		    T1j = KP707106781 * (T11 - T10);
		    T1k = T1i + T1j;
		    T1q = T1i - T1j;
		    T1m = KP707106781 * (TN - TW);
		    T1n = T1a - T13;
		    T1o = T1m + T1n;
		    T1s = T1n - T1m;
	       }
	       {
		    E T1h, T1l, T1p, T1r;
		    T1h = W[2];
		    T1l = W[3];
		    rio[WS(ios, 2)] = FNMS(T1l, T1o, T1h * T1k);
		    iio[-WS(ios, 13)] = FMA(T1l, T1k, T1h * T1o);
		    T1p = W[18];
		    T1r = W[19];
		    rio[WS(ios, 10)] = FNMS(T1r, T1s, T1p * T1q);
		    iio[-WS(ios, 5)] = FMA(T1r, T1q, T1p * T1s);
	       }
	  }
	  {
	       E T2A, T2I, T2E, T2G;
	       {
		    E T2y, T2z, T2C, T2D;
		    T2y = T2k - T2n;
		    T2z = T23 - T1W;
		    T2A = T2y + T2z;
		    T2I = T2y - T2z;
		    T2C = T2p - T2q;
		    T2D = T2e - T2b;
		    T2E = T2C + T2D;
		    T2G = T2D - T2C;
	       }
	       {
		    E T2x, T2B, T2F, T2H;
		    T2x = W[24];
		    T2B = W[25];
		    iio[-WS(ios, 2)] = FMA(T2x, T2A, T2B * T2E);
		    rio[WS(ios, 13)] = FNMS(T2B, T2A, T2x * T2E);
		    T2F = W[8];
		    T2H = W[9];
		    rio[WS(ios, 5)] = FNMS(T2H, T2I, T2F * T2G);
		    iio[-WS(ios, 10)] = FMA(T2F, T2I, T2H * T2G);
	       }
	  }
	  {
	       E T1y, T1G, T1E, T1I;
	       {
		    E T1u, T1x, T1A, T1D;
		    T1u = T7 - Te;
		    T1x = T1v - T1w;
		    T1y = T1u + T1x;
		    T1G = T1u - T1x;
		    T1A = Tt - Tm;
		    T1D = T1B - T1C;
		    T1E = T1A + T1D;
		    T1I = T1D - T1A;
	       }
	       {
		    E T1t, T1z, T1F, T1H;
		    T1t = W[22];
		    T1z = W[23];
		    rio[WS(ios, 12)] = FNMS(T1z, T1E, T1t * T1y);
		    iio[-WS(ios, 3)] = FMA(T1z, T1y, T1t * T1E);
		    T1F = W[6];
		    T1H = W[7];
		    rio[WS(ios, 4)] = FNMS(T1H, T1I, T1F * T1G);
		    iio[-WS(ios, 11)] = FMA(T1H, T1G, T1F * T1I);
	       }
	  }
	  {
	       E T3a, T3i, T3e, T3g;
	       {
		    E T38, T39, T3c, T3d;
		    T38 = T2S - T2P;
		    T39 = T30 - T2Z;
		    T3a = T38 + T39;
		    T3i = T39 - T38;
		    T3c = T2K - T2L;
		    T3d = T2W - T2X;
		    T3e = T3c + T3d;
		    T3g = T3c - T3d;
	       }
	       {
		    E T37, T3b, T3f, T3h;
		    T37 = W[4];
		    T3b = W[5];
		    rio[WS(ios, 3)] = FNMS(T3b, T3e, T37 * T3a);
		    iio[-WS(ios, 12)] = FMA(T37, T3e, T3b * T3a);
		    T3f = W[20];
		    T3h = W[21];
		    iio[-WS(ios, 4)] = FMA(T3f, T3g, T3h * T3i);
		    rio[WS(ios, 11)] = FNMS(T3h, T3g, T3f * T3i);
	       }
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 16},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 16, "hb_16", twinstr, &GENUS, {136, 46, 38, 0}, 0, 0, 0 };

void X(codelet_hb_16) (planner *p) {
     X(khc2hc_register) (p, hb_16, &desc);
}
#endif				/* HAVE_FMA */