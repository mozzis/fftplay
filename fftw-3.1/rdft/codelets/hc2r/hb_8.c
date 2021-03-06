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
/* Generated on Fri Jan 27 20:42:22 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -dif -name hb_8 -include hb.h */

/*
 * This function contains 66 FP additions, 36 FP multiplications,
 * (or, 44 additions, 14 multiplications, 22 fused multiply/add),
 * 48 stack variables, and 32 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hb.h"

static const R *hb_8(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 14, MAKE_VOLATILE_STRIDE(ios)) {
	  E TF, TC, TB, TE, TD, TG;
	  {
	       E T15, T1i, Tg, T7, TS, T1n, TL, Ty, Te, Tr, Tn, TM, T1j, TZ, T1o;
	       E T18;
	       {
		    E T4, TQ, T3, T14, Tu, T5, Tv, Tw;
		    {
			 E T1, T2, Ts, Tt;
			 T1 = rio[0];
			 T2 = iio[-WS(ios, 4)];
			 Ts = iio[0];
			 Tt = rio[WS(ios, 4)];
			 T4 = rio[WS(ios, 2)];
			 TQ = T1 - T2;
			 T3 = T1 + T2;
			 T14 = Ts + Tt;
			 Tu = Ts - Tt;
			 T5 = iio[-WS(ios, 6)];
			 Tv = iio[-WS(ios, 2)];
			 Tw = rio[WS(ios, 6)];
		    }
		    {
			 E Tb, TX, Ta, TW, Tj, Tc, Tk, Tl;
			 {
			      E T8, T9, Th, Ti;
			      T8 = rio[WS(ios, 1)];
			      {
				   E T13, T6, TR, Tx;
				   T13 = T4 - T5;
				   T6 = T4 + T5;
				   TR = Tv + Tw;
				   Tx = Tv - Tw;
				   T15 = T13 + T14;
				   T1i = T14 - T13;
				   Tg = T3 - T6;
				   T7 = T3 + T6;
				   TS = TQ - TR;
				   T1n = TQ + TR;
				   TL = Tu + Tx;
				   Ty = Tu - Tx;
				   T9 = iio[-WS(ios, 5)];
			      }
			      Th = iio[-WS(ios, 1)];
			      Ti = rio[WS(ios, 5)];
			      Tb = iio[-WS(ios, 7)];
			      TX = T8 - T9;
			      Ta = T8 + T9;
			      TW = Th + Ti;
			      Tj = Th - Ti;
			      Tc = rio[WS(ios, 3)];
			      Tk = iio[-WS(ios, 3)];
			      Tl = rio[WS(ios, 7)];
			 }
			 {
			      E TY, T16, TT, Td, TU, Tm, T17, TV;
			      TY = TW - TX;
			      T16 = TX + TW;
			      TT = Tb - Tc;
			      Td = Tb + Tc;
			      TU = Tl + Tk;
			      Tm = Tk - Tl;
			      Te = Ta + Td;
			      Tr = Td - Ta;
			      T17 = TT + TU;
			      TV = TT - TU;
			      Tn = Tj - Tm;
			      TM = Tj + Tm;
			      T1j = TY + TV;
			      TZ = TV - TY;
			      T1o = T16 + T17;
			      T18 = T16 - T17;
			 }
		    }
	       }
	       {
		    E T1s, T1v, T1u, T1w, T1t;
		    rio[0] = T7 + Te;
		    iio[-WS(ios, 7)] = TM + TL;
		    {
			 E T1p, T1k, T1h, T1m, T1q, T1l, T1r;
			 T1s = FNMS(KP707106781, T1o, T1n);
			 T1p = FMA(KP707106781, T1o, T1n);
			 T1k = FMA(KP707106781, T1j, T1i);
			 T1v = FNMS(KP707106781, T1j, T1i);
			 T1h = W[12];
			 T1m = W[13];
			 T1q = T1h * T1p;
			 T1l = T1h * T1k;
			 T1r = W[4];
			 T1u = W[5];
			 rio[WS(ios, 7)] = FNMS(T1m, T1k, T1q);
			 iio[0] = FMA(T1m, T1p, T1l);
			 T1w = T1r * T1v;
			 T1t = T1r * T1s;
		    }
		    {
			 E Tz, Tf, Tq, Tp, TA;
			 {
			      E TN, TI, TH, TK, To, TJ, TO;
			      TN = TL - TM;
			      TI = T7 - Te;
			      iio[-WS(ios, 4)] = FMA(T1u, T1s, T1w);
			      rio[WS(ios, 3)] = FNMS(T1u, T1v, T1t);
			      TH = W[6];
			      TK = W[7];
			      TF = Ty - Tr;
			      Tz = Tr + Ty;
			      To = Tg + Tn;
			      TC = Tg - Tn;
			      TJ = TH * TI;
			      TO = TK * TI;
			      Tf = W[10];
			      Tq = W[11];
			      rio[WS(ios, 4)] = FNMS(TK, TN, TJ);
			      iio[-WS(ios, 3)] = FMA(TH, TN, TO);
			      Tp = Tf * To;
			      TA = Tq * To;
			 }
			 {
			      E T1f, T10, T19, T1c, TP, T12;
			      T1f = FNMS(KP707106781, TZ, TS);
			      T10 = FMA(KP707106781, TZ, TS);
			      T19 = FMA(KP707106781, T18, T15);
			      T1c = FNMS(KP707106781, T18, T15);
			      rio[WS(ios, 6)] = FNMS(Tq, Tz, Tp);
			      iio[-WS(ios, 1)] = FMA(Tf, Tz, TA);
			      TP = W[0];
			      T12 = W[1];
			      {
				   E T1e, T1g, T1d, T1a, T11, T1b;
				   T1a = TP * T19;
				   T11 = TP * T10;
				   T1b = W[8];
				   T1e = W[9];
				   iio[-WS(ios, 6)] = FMA(T12, T10, T1a);
				   rio[WS(ios, 1)] = FNMS(T12, T19, T11);
				   T1g = T1b * T1f;
				   T1d = T1b * T1c;
				   rio[WS(ios, 5)] = FNMS(T1e, T1c, T1g);
				   iio[-WS(ios, 2)] = FMA(T1e, T1f, T1d);
				   TB = W[2];
				   TE = W[3];
			      }
			 }
		    }
	       }
	  }
	  TD = TB * TC;
	  TG = TE * TC;
	  rio[WS(ios, 2)] = FNMS(TE, TF, TD);
	  iio[-WS(ios, 5)] = FMA(TB, TF, TG);
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 8},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 8, "hb_8", twinstr, &GENUS, {44, 14, 22, 0}, 0, 0, 0 };

void X(codelet_hb_8) (planner *p) {
     X(khc2hc_register) (p, hb_8, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2hc -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -dif -name hb_8 -include hb.h */

/*
 * This function contains 66 FP additions, 32 FP multiplications,
 * (or, 52 additions, 18 multiplications, 14 fused multiply/add),
 * 30 stack variables, and 32 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2hc.ml,v 1.15 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hb.h"

static const R *hb_8(R *rio, R *iio, const R *W, stride ios, INT m, INT dist)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     INT i;
     for (i = m - 2; i > 0; i = i - 2, rio = rio + dist, iio = iio - dist, W = W + 14, MAKE_VOLATILE_STRIDE(ios)) {
	  E T7, T18, T1d, Tg, Tx, TT, TY, TG, Te, TZ, T10, Tn, Tq, TM, TP;
	  E TH;
	  {
	       E T3, TR, Tt, TX, T6, TW, Tw, TS;
	       {
		    E T1, T2, Tr, Ts;
		    T1 = rio[0];
		    T2 = iio[-WS(ios, 4)];
		    T3 = T1 + T2;
		    TR = T1 - T2;
		    Tr = iio[0];
		    Ts = rio[WS(ios, 4)];
		    Tt = Tr - Ts;
		    TX = Tr + Ts;
	       }
	       {
		    E T4, T5, Tu, Tv;
		    T4 = rio[WS(ios, 2)];
		    T5 = iio[-WS(ios, 6)];
		    T6 = T4 + T5;
		    TW = T4 - T5;
		    Tu = iio[-WS(ios, 2)];
		    Tv = rio[WS(ios, 6)];
		    Tw = Tu - Tv;
		    TS = Tu + Tv;
	       }
	       T7 = T3 + T6;
	       T18 = TX - TW;
	       T1d = TR + TS;
	       Tg = T3 - T6;
	       Tx = Tt - Tw;
	       TT = TR - TS;
	       TY = TW + TX;
	       TG = Tt + Tw;
	  }
	  {
	       E Ta, TO, Tj, TN, Td, TK, Tm, TL;
	       {
		    E T8, T9, Th, Ti;
		    T8 = rio[WS(ios, 1)];
		    T9 = iio[-WS(ios, 5)];
		    Ta = T8 + T9;
		    TO = T8 - T9;
		    Th = iio[-WS(ios, 1)];
		    Ti = rio[WS(ios, 5)];
		    Tj = Th - Ti;
		    TN = Th + Ti;
	       }
	       {
		    E Tb, Tc, Tk, Tl;
		    Tb = iio[-WS(ios, 7)];
		    Tc = rio[WS(ios, 3)];
		    Td = Tb + Tc;
		    TK = Tb - Tc;
		    Tk = iio[-WS(ios, 3)];
		    Tl = rio[WS(ios, 7)];
		    Tm = Tk - Tl;
		    TL = Tl + Tk;
	       }
	       Te = Ta + Td;
	       TZ = TO + TN;
	       T10 = TK + TL;
	       Tn = Tj - Tm;
	       Tq = Td - Ta;
	       TM = TK - TL;
	       TP = TN - TO;
	       TH = Tj + Tm;
	  }
	  rio[0] = T7 + Te;
	  iio[-WS(ios, 7)] = TH + TG;
	  {
	       E To, Ty, Tf, Tp;
	       To = Tg + Tn;
	       Ty = Tq + Tx;
	       Tf = W[10];
	       Tp = W[11];
	       rio[WS(ios, 6)] = FNMS(Tp, Ty, Tf * To);
	       iio[-WS(ios, 1)] = FMA(Tp, To, Tf * Ty);
	  }
	  {
	       E TE, TI, TD, TF;
	       TE = T7 - Te;
	       TI = TG - TH;
	       TD = W[6];
	       TF = W[7];
	       rio[WS(ios, 4)] = FNMS(TF, TI, TD * TE);
	       iio[-WS(ios, 3)] = FMA(TF, TE, TD * TI);
	  }
	  {
	       E TA, TC, Tz, TB;
	       TA = Tg - Tn;
	       TC = Tx - Tq;
	       Tz = W[2];
	       TB = W[3];
	       rio[WS(ios, 2)] = FNMS(TB, TC, Tz * TA);
	       iio[-WS(ios, 5)] = FMA(TB, TA, Tz * TC);
	  }
	  {
	       E TU, T16, T12, T14, TQ, T11;
	       TQ = KP707106781 * (TM - TP);
	       TU = TQ + TT;
	       T16 = TT - TQ;
	       T11 = KP707106781 * (TZ - T10);
	       T12 = TY + T11;
	       T14 = TY - T11;
	       {
		    E TJ, TV, T13, T15;
		    TJ = W[0];
		    TV = W[1];
		    rio[WS(ios, 1)] = FNMS(TV, T12, TJ * TU);
		    iio[-WS(ios, 6)] = FMA(TJ, T12, TV * TU);
		    T13 = W[8];
		    T15 = W[9];
		    iio[-WS(ios, 2)] = FMA(T13, T14, T15 * T16);
		    rio[WS(ios, 5)] = FNMS(T15, T14, T13 * T16);
	       }
	  }
	  {
	       E T1a, T1i, T1e, T1g, T19, T1c;
	       T19 = KP707106781 * (TP + TM);
	       T1a = T18 + T19;
	       T1i = T18 - T19;
	       T1c = KP707106781 * (TZ + T10);
	       T1e = T1c + T1d;
	       T1g = T1d - T1c;
	       {
		    E T17, T1b, T1f, T1h;
		    T17 = W[12];
		    T1b = W[13];
		    iio[0] = FMA(T17, T1a, T1b * T1e);
		    rio[WS(ios, 7)] = FNMS(T1b, T1a, T17 * T1e);
		    T1f = W[4];
		    T1h = W[5];
		    rio[WS(ios, 3)] = FNMS(T1h, T1i, T1f * T1g);
		    iio[-WS(ios, 4)] = FMA(T1f, T1i, T1h * T1g);
	       }
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     {TW_FULL, 0, 8},
     {TW_NEXT, 1, 0}
};

static const hc2hc_desc desc = { 8, "hb_8", twinstr, &GENUS, {52, 18, 14, 0}, 0, 0, 0 };

void X(codelet_hb_8) (planner *p) {
     X(khc2hc_register) (p, hb_8, &desc);
}
#endif				/* HAVE_FMA */
