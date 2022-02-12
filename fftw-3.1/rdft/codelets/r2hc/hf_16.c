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
/* Generated on Fri Jan 27 20:19:22 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 16 -dit -name hf_16 -include hf.h */

/*
 * This function contains 174 FP additions, 100 FP multiplications,
 * (or, 104 additions, 30 multiplications, 70 fused multiply/add),
 * 97 stack variables, and 64 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hf.h"

static const R *hf_16(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 30, MAKE_VOLATILE_STRIDE(ios)) {
	  E T3G, T3F;
	  {
	       E T3z, T3o, T8, T1I, T2p, T35, T2r, T1s, T3k, T1N, T3A, Tl, T2w, T36, T2k;
	       E T1F, T1T, T2V, T1U, Tz, T2h, T31, T2a, T1e, TQ, TT, T21, T2W, T1W, TM;
	       E TR, T25, TW, TZ, TV, TS, TY;
	       {
		    E T1u, T1x, T1v, T2s, T1A, T1D, T1z, T1w, T1C;
		    {
			 E T1h, T1k, T1n, T2l, T1i, T1q, T1m, T1j, T1p;
			 {
			      E T1, T3n, T3, T6, T2, T5;
			      T1 = rio[0];
			      T3n = iio[-WS(ios, 15)];
			      T3 = rio[WS(ios, 8)];
			      T6 = iio[-WS(ios, 7)];
			      T2 = W[14];
			      T5 = W[15];
			      {
				   E T3l, T4, T1g, T3m, T7;
				   T1h = rio[WS(ios, 15)];
				   T1k = iio[0];
				   T3l = T2 * T6;
				   T4 = T2 * T3;
				   T1g = W[28];
				   T1n = rio[WS(ios, 7)];
				   T3m = FNMS(T5, T3, T3l);
				   T7 = FMA(T5, T6, T4);
				   T2l = T1g * T1k;
				   T1i = T1g * T1h;
				   T3z = T3n - T3m;
				   T3o = T3m + T3n;
				   T8 = T1 + T7;
				   T1I = T1 - T7;
				   T1q = iio[-WS(ios, 8)];
				   T1m = W[12];
			      }
			      T1j = W[29];
			      T1p = W[13];
			 }
			 {
			      E Ta, Td, Tb, T1J, Tg, Tj, Tf, Tc, Ti;
			      {
				   E T2m, T1l, T2o, T1r, T2n, T1o, T9;
				   Ta = rio[WS(ios, 4)];
				   T2n = T1m * T1q;
				   T1o = T1m * T1n;
				   T2m = FNMS(T1j, T1h, T2l);
				   T1l = FMA(T1j, T1k, T1i);
				   T2o = FNMS(T1p, T1n, T2n);
				   T1r = FMA(T1p, T1q, T1o);
				   Td = iio[-WS(ios, 11)];
				   T9 = W[6];
				   T2p = T2m - T2o;
				   T35 = T2m + T2o;
				   T2r = T1l - T1r;
				   T1s = T1l + T1r;
				   Tb = T9 * Ta;
				   T1J = T9 * Td;
			      }
			      Tg = rio[WS(ios, 12)];
			      Tj = iio[-WS(ios, 3)];
			      Tf = W[22];
			      Tc = W[7];
			      Ti = W[23];
			      {
				   E T1K, Te, T1M, Tk, T1L, Th, T1t;
				   T1u = rio[WS(ios, 3)];
				   T1L = Tf * Tj;
				   Th = Tf * Tg;
				   T1K = FNMS(Tc, Ta, T1J);
				   Te = FMA(Tc, Td, Tb);
				   T1M = FNMS(Ti, Tg, T1L);
				   Tk = FMA(Ti, Tj, Th);
				   T1x = iio[-WS(ios, 12)];
				   T1t = W[4];
				   T3k = T1K + T1M;
				   T1N = T1K - T1M;
				   T3A = Te - Tk;
				   Tl = Te + Tk;
				   T1v = T1t * T1u;
				   T2s = T1t * T1x;
			      }
			      T1A = rio[WS(ios, 11)];
			      T1D = iio[-WS(ios, 4)];
			      T1z = W[20];
			      T1w = W[5];
			      T1C = W[21];
			 }
		    }
		    {
			 E T13, T16, T14, T2d, T19, T1c, T18, T15, T1b;
			 {
			      E To, Tr, Tp, T1P, Tu, Tx, Tt, Tq, Tw;
			      {
				   E T2t, T1y, T2v, T1E, T2u, T1B, Tn;
				   To = rio[WS(ios, 2)];
				   T2u = T1z * T1D;
				   T1B = T1z * T1A;
				   T2t = FNMS(T1w, T1u, T2s);
				   T1y = FMA(T1w, T1x, T1v);
				   T2v = FNMS(T1C, T1A, T2u);
				   T1E = FMA(T1C, T1D, T1B);
				   Tr = iio[-WS(ios, 13)];
				   Tn = W[2];
				   T2w = T2t - T2v;
				   T36 = T2t + T2v;
				   T2k = T1E - T1y;
				   T1F = T1y + T1E;
				   Tp = Tn * To;
				   T1P = Tn * Tr;
			      }
			      Tu = rio[WS(ios, 10)];
			      Tx = iio[-WS(ios, 5)];
			      Tt = W[18];
			      Tq = W[3];
			      Tw = W[19];
			      {
				   E T1Q, Ts, T1S, Ty, T1R, Tv, T12;
				   T13 = rio[WS(ios, 5)];
				   T1R = Tt * Tx;
				   Tv = Tt * Tu;
				   T1Q = FNMS(Tq, To, T1P);
				   Ts = FMA(Tq, Tr, Tp);
				   T1S = FNMS(Tw, Tu, T1R);
				   Ty = FMA(Tw, Tx, Tv);
				   T16 = iio[-WS(ios, 10)];
				   T12 = W[8];
				   T1T = T1Q - T1S;
				   T2V = T1Q + T1S;
				   T1U = Ts - Ty;
				   Tz = Ts + Ty;
				   T14 = T12 * T13;
				   T2d = T12 * T16;
			      }
			      T19 = rio[WS(ios, 13)];
			      T1c = iio[-WS(ios, 2)];
			      T18 = W[24];
			      T15 = W[9];
			      T1b = W[25];
			 }
			 {
			      E TB, TE, TC, T1X, TH, TK, TG, TD, TJ;
			      {
				   E T2e, T17, T2g, T1d, T2f, T1a, TA;
				   TB = rio[WS(ios, 14)];
				   T2f = T18 * T1c;
				   T1a = T18 * T19;
				   T2e = FNMS(T15, T13, T2d);
				   T17 = FMA(T15, T16, T14);
				   T2g = FNMS(T1b, T19, T2f);
				   T1d = FMA(T1b, T1c, T1a);
				   TE = iio[-WS(ios, 1)];
				   TA = W[26];
				   T2h = T2e - T2g;
				   T31 = T2e + T2g;
				   T2a = T17 - T1d;
				   T1e = T17 + T1d;
				   TC = TA * TB;
				   T1X = TA * TE;
			      }
			      TH = rio[WS(ios, 6)];
			      TK = iio[-WS(ios, 9)];
			      TG = W[10];
			      TD = W[27];
			      TJ = W[11];
			      {
				   E T1Y, TF, T20, TL, T1Z, TI, TP;
				   TQ = rio[WS(ios, 1)];
				   T1Z = TG * TK;
				   TI = TG * TH;
				   T1Y = FNMS(TD, TB, T1X);
				   TF = FMA(TD, TE, TC);
				   T20 = FNMS(TJ, TH, T1Z);
				   TL = FMA(TJ, TK, TI);
				   TT = iio[-WS(ios, 14)];
				   TP = W[0];
				   T21 = T1Y - T20;
				   T2W = T1Y + T20;
				   T1W = TF - TL;
				   TM = TF + TL;
				   TR = TP * TQ;
				   T25 = TP * TT;
			      }
			      TW = rio[WS(ios, 9)];
			      TZ = iio[-WS(ios, 6)];
			      TV = W[16];
			      TS = W[1];
			      TY = W[17];
			 }
		    }
	       }
	       {
		    E T2U, T3t, T2X, T29, T2c, T3e, TO, T3u, T2Z, T34, T32, T3f, T3s, T3q, T3r;
		    E T1H, T3g, T37;
		    {
			 E T3j, T30, T11, T3p, T1f, T1G;
			 {
			      E Tm, T26, TU, T28, T10, TN, T27, TX;
			      T2U = T8 - Tl;
			      Tm = T8 + Tl;
			      T27 = TV * TZ;
			      TX = TV * TW;
			      T26 = FNMS(TS, TQ, T25);
			      TU = FMA(TS, TT, TR);
			      T28 = FNMS(TY, TW, T27);
			      T10 = FMA(TY, TZ, TX);
			      TN = Tz + TM;
			      T3t = TM - Tz;
			      T2X = T2V - T2W;
			      T3j = T2V + T2W;
			      T29 = T26 - T28;
			      T30 = T26 + T28;
			      T2c = TU - T10;
			      T11 = TU + T10;
			      T3e = Tm - TN;
			      TO = Tm + TN;
			      T3p = T3k + T3o;
			      T3u = T3o - T3k;
			 }
			 T2Z = T11 - T1e;
			 T1f = T11 + T1e;
			 T1G = T1s + T1F;
			 T34 = T1s - T1F;
			 T32 = T30 - T31;
			 T3f = T30 + T31;
			 T3s = T3p - T3j;
			 T3q = T3j + T3p;
			 T3r = T1G - T1f;
			 T1H = T1f + T1G;
			 T3g = T35 + T36;
			 T37 = T35 - T36;
		    }
		    {
			 E T3a, T2Y, T3x, T3v, T3b, T33, T3i, T3h;
			 iio[-WS(ios, 4)] = T3r + T3s;
			 rio[WS(ios, 12)] = T3r - T3s;
			 rio[0] = TO + T1H;
			 iio[-WS(ios, 8)] = TO - T1H;
			 T3i = T3f + T3g;
			 T3h = T3f - T3g;
			 iio[0] = T3i + T3q;
			 rio[WS(ios, 8)] = T3i - T3q;
			 rio[WS(ios, 4)] = T3e + T3h;
			 iio[-WS(ios, 12)] = T3e - T3h;
			 T3a = T2U - T2X;
			 T2Y = T2U + T2X;
			 T3x = T3u - T3t;
			 T3v = T3t + T3u;
			 T3b = T32 - T2Z;
			 T33 = T2Z + T32;
			 {
			      E T2E, T1O, T3B, T3H, T2x, T2q, T3C, T23, T2S, T2O, T2K, T2J, T3I, T2H, T2B;
			      E T2j;
			      {
				   E T2F, T1V, T22, T2G, T3c, T38;
				   T2E = T1I + T1N;
				   T1O = T1I - T1N;
				   T3B = T3z - T3A;
				   T3H = T3A + T3z;
				   T3c = T34 + T37;
				   T38 = T34 - T37;
				   T2F = T1U + T1T;
				   T1V = T1T - T1U;
				   {
					E T3d, T3w, T3y, T39;
					T3d = T3b - T3c;
					T3w = T3b + T3c;
					T3y = T38 - T33;
					T39 = T33 + T38;
					rio[WS(ios, 6)] = FMA(KP707106781, T3d, T3a);
					iio[-WS(ios, 14)] = FNMS(KP707106781, T3d, T3a);
					iio[-WS(ios, 2)] = FMA(KP707106781, T3w, T3v);
					rio[WS(ios, 10)] = FMS(KP707106781, T3w, T3v);
					iio[-WS(ios, 6)] = FMA(KP707106781, T3y, T3x);
					rio[WS(ios, 14)] = FMS(KP707106781, T3y, T3x);
					rio[WS(ios, 2)] = FMA(KP707106781, T39, T2Y);
					iio[-WS(ios, 10)] = FNMS(KP707106781, T39, T2Y);
					T22 = T1W + T21;
					T2G = T1W - T21;
				   }
				   {
					E T2M, T2N, T2b, T2i;
					T2x = T2r - T2w;
					T2M = T2r + T2w;
					T2N = T2p + T2k;
					T2q = T2k - T2p;
					T3C = T1V + T22;
					T23 = T1V - T22;
					T2S = FMA(KP414213562, T2M, T2N);
					T2O = FNMS(KP414213562, T2N, T2M);
					T2K = T29 - T2a;
					T2b = T29 + T2a;
					T2i = T2c - T2h;
					T2J = T2c + T2h;
					T3I = T2G - T2F;
					T2H = T2F + T2G;
					T2B = FNMS(KP414213562, T2b, T2i);
					T2j = FMA(KP414213562, T2i, T2b);
				   }
			      }
			      {
				   E T2R, T2L, T3L, T3M;
				   {
					E T2A, T24, T2C, T2y, T3J, T3K, T2D, T2z;
					T2A = FNMS(KP707106781, T23, T1O);
					T24 = FMA(KP707106781, T23, T1O);
					T2R = FNMS(KP414213562, T2J, T2K);
					T2L = FMA(KP414213562, T2K, T2J);
					T2C = FNMS(KP414213562, T2q, T2x);
					T2y = FMA(KP414213562, T2x, T2q);
					T3J = FMA(KP707106781, T3I, T3H);
					T3L = FNMS(KP707106781, T3I, T3H);
					T3K = T2C - T2B;
					T2D = T2B + T2C;
					T3M = T2y - T2j;
					T2z = T2j + T2y;
					iio[-WS(ios, 3)] = FMA(KP923879532, T3K, T3J);
					rio[WS(ios, 11)] = FMS(KP923879532, T3K, T3J);
					rio[WS(ios, 3)] = FMA(KP923879532, T2z, T24);
					iio[-WS(ios, 11)] = FNMS(KP923879532, T2z, T24);
					iio[-WS(ios, 15)] = FMA(KP923879532, T2D, T2A);
					rio[WS(ios, 7)] = FNMS(KP923879532, T2D, T2A);
				   }
				   {
					E T2Q, T3D, T3E, T2T, T2I, T2P;
					T2Q = FNMS(KP707106781, T2H, T2E);
					T2I = FMA(KP707106781, T2H, T2E);
					T2P = T2L + T2O;
					T3G = T2O - T2L;
					T3F = FNMS(KP707106781, T3C, T3B);
					T3D = FMA(KP707106781, T3C, T3B);
					iio[-WS(ios, 7)] = FMA(KP923879532, T3M, T3L);
					rio[WS(ios, 15)] = FMS(KP923879532, T3M, T3L);
					rio[WS(ios, 1)] = FMA(KP923879532, T2P, T2I);
					iio[-WS(ios, 9)] = FNMS(KP923879532, T2P, T2I);
					T3E = T2R + T2S;
					T2T = T2R - T2S;
					iio[-WS(ios, 1)] = FMA(KP923879532, T3E, T3D);
					rio[WS(ios, 9)] = FMS(KP923879532, T3E, T3D);
					rio[WS(ios, 5)] = FMA(KP923879532, T2T, T2Q);
					iio[-WS(ios, 13)] = FNMS(KP923879532, T2T, T2Q);
				   }
			      }
			 }
		    }
	       }
	  }
	  iio[-WS(ios, 5)] = FMA(KP923879532, T3G, T3F);
	  rio[WS(ios, 13)] = FMS(KP923879532, T3G, T3F);
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 16},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 16, "hf_16", twinstr, &GENUS, {104, 30, 70, 0}, 0, 0, 0 };

void X(codelet_hf_16) (planner *p) {
     X(khc2hc_register) (p, hf_16, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc -compact -variables 4 -pipeline-latency 4 -n 16 -dit -name hf_16 -include hf.h */

/*
 * This function contains 174 FP additions, 84 FP multiplications,
 * (or, 136 additions, 46 multiplications, 38 fused multiply/add),
 * 52 stack variables, and 64 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hf.h"

static const R *hf_16(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 30, MAKE_VOLATILE_STRIDE(ios)) {
	  E T7, T37, T1t, T2U, Ti, T38, T1w, T2R, Tu, T2s, T1C, T2c, TF, T2t, T1H;
	  E T2d, TS, T13, T2w, T2x, T2y, T2z, T1O, T2g, T1T, T2h, T1f, T1q, T2B, T2C;
	  E T2D, T2E, T1Z, T2j, T24, T2k;
	  {
	       E T1, T2T, T6, T2S;
	       T1 = rio[0];
	       T2T = iio[-WS(ios, 15)];
	       {
		    E T3, T5, T2, T4;
		    T3 = rio[WS(ios, 8)];
		    T5 = iio[-WS(ios, 7)];
		    T2 = W[14];
		    T4 = W[15];
		    T6 = FMA(T2, T3, T4 * T5);
		    T2S = FNMS(T4, T3, T2 * T5);
	       }
	       T7 = T1 + T6;
	       T37 = T2T - T2S;
	       T1t = T1 - T6;
	       T2U = T2S + T2T;
	  }
	  {
	       E Tc, T1u, Th, T1v;
	       {
		    E T9, Tb, T8, Ta;
		    T9 = rio[WS(ios, 4)];
		    Tb = iio[-WS(ios, 11)];
		    T8 = W[6];
		    Ta = W[7];
		    Tc = FMA(T8, T9, Ta * Tb);
		    T1u = FNMS(Ta, T9, T8 * Tb);
	       }
	       {
		    E Te, Tg, Td, Tf;
		    Te = rio[WS(ios, 12)];
		    Tg = iio[-WS(ios, 3)];
		    Td = W[22];
		    Tf = W[23];
		    Th = FMA(Td, Te, Tf * Tg);
		    T1v = FNMS(Tf, Te, Td * Tg);
	       }
	       Ti = Tc + Th;
	       T38 = Tc - Th;
	       T1w = T1u - T1v;
	       T2R = T1u + T1v;
	  }
	  {
	       E To, T1y, Tt, T1z, T1A, T1B;
	       {
		    E Tl, Tn, Tk, Tm;
		    Tl = rio[WS(ios, 2)];
		    Tn = iio[-WS(ios, 13)];
		    Tk = W[2];
		    Tm = W[3];
		    To = FMA(Tk, Tl, Tm * Tn);
		    T1y = FNMS(Tm, Tl, Tk * Tn);
	       }
	       {
		    E Tq, Ts, Tp, Tr;
		    Tq = rio[WS(ios, 10)];
		    Ts = iio[-WS(ios, 5)];
		    Tp = W[18];
		    Tr = W[19];
		    Tt = FMA(Tp, Tq, Tr * Ts);
		    T1z = FNMS(Tr, Tq, Tp * Ts);
	       }
	       Tu = To + Tt;
	       T2s = T1y + T1z;
	       T1A = T1y - T1z;
	       T1B = To - Tt;
	       T1C = T1A - T1B;
	       T2c = T1B + T1A;
	  }
	  {
	       E Tz, T1E, TE, T1F, T1D, T1G;
	       {
		    E Tw, Ty, Tv, Tx;
		    Tw = rio[WS(ios, 14)];
		    Ty = iio[-WS(ios, 1)];
		    Tv = W[26];
		    Tx = W[27];
		    Tz = FMA(Tv, Tw, Tx * Ty);
		    T1E = FNMS(Tx, Tw, Tv * Ty);
	       }
	       {
		    E TB, TD, TA, TC;
		    TB = rio[WS(ios, 6)];
		    TD = iio[-WS(ios, 9)];
		    TA = W[10];
		    TC = W[11];
		    TE = FMA(TA, TB, TC * TD);
		    T1F = FNMS(TC, TB, TA * TD);
	       }
	       TF = Tz + TE;
	       T2t = T1E + T1F;
	       T1D = Tz - TE;
	       T1G = T1E - T1F;
	       T1H = T1D + T1G;
	       T2d = T1D - T1G;
	  }
	  {
	       E TM, T1K, T12, T1R, TR, T1L, TX, T1Q;
	       {
		    E TJ, TL, TI, TK;
		    TJ = rio[WS(ios, 1)];
		    TL = iio[-WS(ios, 14)];
		    TI = W[0];
		    TK = W[1];
		    TM = FMA(TI, TJ, TK * TL);
		    T1K = FNMS(TK, TJ, TI * TL);
	       }
	       {
		    E TZ, T11, TY, T10;
		    TZ = rio[WS(ios, 13)];
		    T11 = iio[-WS(ios, 2)];
		    TY = W[24];
		    T10 = W[25];
		    T12 = FMA(TY, TZ, T10 * T11);
		    T1R = FNMS(T10, TZ, TY * T11);
	       }
	       {
		    E TO, TQ, TN, TP;
		    TO = rio[WS(ios, 9)];
		    TQ = iio[-WS(ios, 6)];
		    TN = W[16];
		    TP = W[17];
		    TR = FMA(TN, TO, TP * TQ);
		    T1L = FNMS(TP, TO, TN * TQ);
	       }
	       {
		    E TU, TW, TT, TV;
		    TU = rio[WS(ios, 5)];
		    TW = iio[-WS(ios, 10)];
		    TT = W[8];
		    TV = W[9];
		    TX = FMA(TT, TU, TV * TW);
		    T1Q = FNMS(TV, TU, TT * TW);
	       }
	       TS = TM + TR;
	       T13 = TX + T12;
	       T2w = TS - T13;
	       T2x = T1K + T1L;
	       T2y = T1Q + T1R;
	       T2z = T2x - T2y;
	       {
		    E T1M, T1N, T1P, T1S;
		    T1M = T1K - T1L;
		    T1N = TX - T12;
		    T1O = T1M + T1N;
		    T2g = T1M - T1N;
		    T1P = TM - TR;
		    T1S = T1Q - T1R;
		    T1T = T1P - T1S;
		    T2h = T1P + T1S;
	       }
	  }
	  {
	       E T19, T20, T1p, T1X, T1e, T21, T1k, T1W;
	       {
		    E T16, T18, T15, T17;
		    T16 = rio[WS(ios, 15)];
		    T18 = iio[0];
		    T15 = W[28];
		    T17 = W[29];
		    T19 = FMA(T15, T16, T17 * T18);
		    T20 = FNMS(T17, T16, T15 * T18);
	       }
	       {
		    E T1m, T1o, T1l, T1n;
		    T1m = rio[WS(ios, 11)];
		    T1o = iio[-WS(ios, 4)];
		    T1l = W[20];
		    T1n = W[21];
		    T1p = FMA(T1l, T1m, T1n * T1o);
		    T1X = FNMS(T1n, T1m, T1l * T1o);
	       }
	       {
		    E T1b, T1d, T1a, T1c;
		    T1b = rio[WS(ios, 7)];
		    T1d = iio[-WS(ios, 8)];
		    T1a = W[12];
		    T1c = W[13];
		    T1e = FMA(T1a, T1b, T1c * T1d);
		    T21 = FNMS(T1c, T1b, T1a * T1d);
	       }
	       {
		    E T1h, T1j, T1g, T1i;
		    T1h = rio[WS(ios, 3)];
		    T1j = iio[-WS(ios, 12)];
		    T1g = W[4];
		    T1i = W[5];
		    T1k = FMA(T1g, T1h, T1i * T1j);
		    T1W = FNMS(T1i, T1h, T1g * T1j);
	       }
	       T1f = T19 + T1e;
	       T1q = T1k + T1p;
	       T2B = T1f - T1q;
	       T2C = T20 + T21;
	       T2D = T1W + T1X;
	       T2E = T2C - T2D;
	       {
		    E T1V, T1Y, T22, T23;
		    T1V = T19 - T1e;
		    T1Y = T1W - T1X;
		    T1Z = T1V - T1Y;
		    T2j = T1V + T1Y;
		    T22 = T20 - T21;
		    T23 = T1k - T1p;
		    T24 = T22 + T23;
		    T2k = T22 - T23;
	       }
	  }
	  {
	       E T1J, T27, T3g, T3i, T26, T3h, T2a, T3d;
	       {
		    E T1x, T1I, T3e, T3f;
		    T1x = T1t - T1w;
		    T1I = KP707106781 * (T1C - T1H);
		    T1J = T1x + T1I;
		    T27 = T1x - T1I;
		    T3e = KP707106781 * (T2d - T2c);
		    T3f = T38 + T37;
		    T3g = T3e + T3f;
		    T3i = T3f - T3e;
	       }
	       {
		    E T1U, T25, T28, T29;
		    T1U = FMA(KP923879532, T1O, KP382683432 * T1T);
		    T25 = FNMS(KP923879532, T24, KP382683432 * T1Z);
		    T26 = T1U + T25;
		    T3h = T25 - T1U;
		    T28 = FNMS(KP923879532, T1T, KP382683432 * T1O);
		    T29 = FMA(KP382683432, T24, KP923879532 * T1Z);
		    T2a = T28 - T29;
		    T3d = T28 + T29;
	       }
	       iio[-WS(ios, 11)] = T1J - T26;
	       rio[WS(ios, 11)] = T3d - T3g;
	       rio[WS(ios, 3)] = T1J + T26;
	       iio[-WS(ios, 3)] = T3d + T3g;
	       iio[-WS(ios, 15)] = T27 - T2a;
	       rio[WS(ios, 15)] = T3h - T3i;
	       rio[WS(ios, 7)] = T27 + T2a;
	       iio[-WS(ios, 7)] = T3h + T3i;
	  }
	  {
	       E T2v, T2H, T32, T34, T2G, T33, T2K, T2Z;
	       {
		    E T2r, T2u, T30, T31;
		    T2r = T7 - Ti;
		    T2u = T2s - T2t;
		    T2v = T2r + T2u;
		    T2H = T2r - T2u;
		    T30 = TF - Tu;
		    T31 = T2U - T2R;
		    T32 = T30 + T31;
		    T34 = T31 - T30;
	       }
	       {
		    E T2A, T2F, T2I, T2J;
		    T2A = T2w + T2z;
		    T2F = T2B - T2E;
		    T2G = KP707106781 * (T2A + T2F);
		    T33 = KP707106781 * (T2F - T2A);
		    T2I = T2z - T2w;
		    T2J = T2B + T2E;
		    T2K = KP707106781 * (T2I - T2J);
		    T2Z = KP707106781 * (T2I + T2J);
	       }
	       iio[-WS(ios, 10)] = T2v - T2G;
	       rio[WS(ios, 10)] = T2Z - T32;
	       rio[WS(ios, 2)] = T2v + T2G;
	       iio[-WS(ios, 2)] = T2Z + T32;
	       iio[-WS(ios, 14)] = T2H - T2K;
	       rio[WS(ios, 14)] = T33 - T34;
	       rio[WS(ios, 6)] = T2H + T2K;
	       iio[-WS(ios, 6)] = T33 + T34;
	  }
	  {
	       E T2f, T2n, T3a, T3c, T2m, T3b, T2q, T35;
	       {
		    E T2b, T2e, T36, T39;
		    T2b = T1t + T1w;
		    T2e = KP707106781 * (T2c + T2d);
		    T2f = T2b + T2e;
		    T2n = T2b - T2e;
		    T36 = KP707106781 * (T1C + T1H);
		    T39 = T37 - T38;
		    T3a = T36 + T39;
		    T3c = T39 - T36;
	       }
	       {
		    E T2i, T2l, T2o, T2p;
		    T2i = FMA(KP382683432, T2g, KP923879532 * T2h);
		    T2l = FNMS(KP382683432, T2k, KP923879532 * T2j);
		    T2m = T2i + T2l;
		    T3b = T2l - T2i;
		    T2o = FNMS(KP382683432, T2h, KP923879532 * T2g);
		    T2p = FMA(KP923879532, T2k, KP382683432 * T2j);
		    T2q = T2o - T2p;
		    T35 = T2o + T2p;
	       }
	       iio[-WS(ios, 9)] = T2f - T2m;
	       rio[WS(ios, 9)] = T35 - T3a;
	       rio[WS(ios, 1)] = T2f + T2m;
	       iio[-WS(ios, 1)] = T35 + T3a;
	       iio[-WS(ios, 13)] = T2n - T2q;
	       rio[WS(ios, 13)] = T3b - T3c;
	       rio[WS(ios, 5)] = T2n + T2q;
	       iio[-WS(ios, 5)] = T3b + T3c;
	  }
	  {
	       E TH, T2L, T2W, T2Y, T1s, T2X, T2O, T2P;
	       {
		    E Tj, TG, T2Q, T2V;
		    Tj = T7 + Ti;
		    TG = Tu + TF;
		    TH = Tj + TG;
		    T2L = Tj - TG;
		    T2Q = T2s + T2t;
		    T2V = T2R + T2U;
		    T2W = T2Q + T2V;
		    T2Y = T2V - T2Q;
	       }
	       {
		    E T14, T1r, T2M, T2N;
		    T14 = TS + T13;
		    T1r = T1f + T1q;
		    T1s = T14 + T1r;
		    T2X = T1r - T14;
		    T2M = T2x + T2y;
		    T2N = T2C + T2D;
		    T2O = T2M - T2N;
		    T2P = T2M + T2N;
	       }
	       iio[-WS(ios, 8)] = TH - T1s;
	       rio[WS(ios, 8)] = T2P - T2W;
	       rio[0] = TH + T1s;
	       iio[0] = T2P + T2W;
	       iio[-WS(ios, 12)] = T2L - T2O;
	       rio[WS(ios, 12)] = T2X - T2Y;
	       rio[WS(ios, 4)] = T2L + T2O;
	       iio[-WS(ios, 4)] = T2X + T2Y;
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 16},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 16, "hf_16", twinstr, &GENUS, {136, 46, 38, 0}, 0, 0, 0 };

void X(codelet_hf_16) (planner *p) {
     X(khc2hc_register) (p, hf_16, &desc);
}
#endif				/* HAVE_FMA */
