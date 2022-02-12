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
/* Generated on Fri Jan 27 20:52:40 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2r -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 16 -name hc2rIII_16 -dft-III -include hc2rIII.h */

/*
 * This function contains 66 FP additions, 36 FP multiplications,
 * (or, 46 additions, 16 multiplications, 20 fused multiply/add),
 * 55 stack variables, and 32 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2rIII.h"

static void hc2rIII_16(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP668178637, +0.668178637919298919997757686523080761552472251);
     DK(KP1_662939224, +1.662939224605090474157576755235811513477121624);
     DK(KP198912367, +0.198912367379658006911597622644676228597850501);
     DK(KP1_961570560, +1.961570560806460898252364472268478073947867462);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E TA, TD, Tv, TG, TE, TF;
	  {
	       E TK, TP, T7, T13, TW, TH, Tj, TC, To, Te, TX, TS, T12, Tt, TB;
	       {
		    E T4, Tf, T3, TU, Tz, T5, Tg, Th;
		    {
			 E T1, T2, Tx, Ty;
			 T1 = ri[0];
			 T2 = ri[WS(ris, 7)];
			 Tx = ii[0];
			 Ty = ii[WS(iis, 7)];
			 T4 = ri[WS(ris, 4)];
			 Tf = T1 - T2;
			 T3 = T1 + T2;
			 TU = Ty - Tx;
			 Tz = Tx + Ty;
			 T5 = ri[WS(ris, 3)];
			 Tg = ii[WS(iis, 4)];
			 Th = ii[WS(iis, 3)];
		    }
		    {
			 E Tb, Tk, Ta, TR, Tn, Tc, Tq, Tr;
			 {
			      E T8, T9, Tl, Tm;
			      T8 = ri[WS(ris, 2)];
			      {
				   E Tw, T6, TV, Ti;
				   Tw = T4 - T5;
				   T6 = T4 + T5;
				   TV = Th - Tg;
				   Ti = Tg + Th;
				   TK = Tw - Tz;
				   TA = Tw + Tz;
				   TP = T3 - T6;
				   T7 = T3 + T6;
				   T13 = TV + TU;
				   TW = TU - TV;
				   TH = Tf + Ti;
				   Tj = Tf - Ti;
				   T9 = ri[WS(ris, 5)];
			      }
			      Tl = ii[WS(iis, 2)];
			      Tm = ii[WS(iis, 5)];
			      Tb = ri[WS(ris, 1)];
			      Tk = T8 - T9;
			      Ta = T8 + T9;
			      TR = Tl - Tm;
			      Tn = Tl + Tm;
			      Tc = ri[WS(ris, 6)];
			      Tq = ii[WS(iis, 1)];
			      Tr = ii[WS(iis, 6)];
			 }
			 TC = Tk + Tn;
			 To = Tk - Tn;
			 {
			      E Tp, Td, TQ, Ts;
			      Tp = Tb - Tc;
			      Td = Tb + Tc;
			      TQ = Tr - Tq;
			      Ts = Tq + Tr;
			      Te = Ta + Td;
			      TX = Ta - Td;
			      TS = TQ - TR;
			      T12 = TR + TQ;
			      Tt = Tp - Ts;
			      TB = Tp + Ts;
			 }
		    }
	       }
	       {
		    E T10, TT, TY, TZ;
		    O[0] = KP2_000000000 * (T7 + Te);
		    O[WS(os, 8)] = KP2_000000000 * (T13 - T12);
		    T10 = TP - TS;
		    TT = TP + TS;
		    TY = TW - TX;
		    TZ = TX + TW;
		    {
			 E T11, T14, TI, TL, Tu;
			 T11 = T7 - Te;
			 T14 = T12 + T13;
			 O[WS(os, 10)] = KP1_847759065 * (FNMS(KP414213562, TT, TY));
			 O[WS(os, 2)] = KP1_847759065 * (FMA(KP414213562, TY, TT));
			 O[WS(os, 12)] = KP1_414213562 * (T14 - T11);
			 O[WS(os, 4)] = KP1_414213562 * (T11 + T14);
			 TD = TB - TC;
			 TI = TC + TB;
			 TL = To - Tt;
			 Tu = To + Tt;
			 {
			      E TO, TJ, TN, TM;
			      O[WS(os, 14)] = -(KP1_847759065 * (FNMS(KP414213562, TZ, T10)));
			      O[WS(os, 6)] = KP1_847759065 * (FMA(KP414213562, T10, TZ));
			      TO = FMA(KP707106781, TI, TH);
			      TJ = FNMS(KP707106781, TI, TH);
			      TN = FMA(KP707106781, TL, TK);
			      TM = FNMS(KP707106781, TL, TK);
			      Tv = FMA(KP707106781, Tu, Tj);
			      TG = FNMS(KP707106781, Tu, Tj);
			      O[WS(os, 7)] = KP1_961570560 * (FMA(KP198912367, TO, TN));
			      O[WS(os, 15)] = -(KP1_961570560 * (FNMS(KP198912367, TN, TO)));
			      O[WS(os, 11)] = KP1_662939224 * (FNMS(KP668178637, TJ, TM));
			      O[WS(os, 3)] = KP1_662939224 * (FMA(KP668178637, TM, TJ));
			 }
		    }
	       }
	  }
	  TE = FNMS(KP707106781, TD, TA);
	  TF = FMA(KP707106781, TD, TA);
	  O[WS(os, 5)] = -(KP1_662939224 * (FNMS(KP668178637, TG, TF)));
	  O[WS(os, 13)] = -(KP1_662939224 * (FMA(KP668178637, TF, TG)));
	  O[WS(os, 9)] = -(KP1_961570560 * (FMA(KP198912367, Tv, TE)));
	  O[WS(os, 1)] = KP1_961570560 * (FNMS(KP198912367, TE, Tv));
     }
}

static const khc2r_desc desc = { 16, "hc2rIII_16", {46, 16, 20, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2rIII_16) (planner *p) {
     X(khc2rIII_register) (p, hc2rIII_16, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2r -compact -variables 4 -pipeline-latency 4 -sign 1 -n 16 -name hc2rIII_16 -dft-III -include hc2rIII.h */

/*
 * This function contains 66 FP additions, 32 FP multiplications,
 * (or, 54 additions, 20 multiplications, 12 fused multiply/add),
 * 40 stack variables, and 32 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2rIII.h"

static void hc2rIII_16(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP1_961570560, +1.961570560806460898252364472268478073947867462);
     DK(KP390180644, +0.390180644032256535696569736954044481855383236);
     DK(KP1_111140466, +1.111140466039204449485661627897065748749874382);
     DK(KP1_662939224, +1.662939224605090474157576755235811513477121624);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP765366864, +0.765366864730179543456919968060797733522689125);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E T7, TW, T13, Tj, TD, TK, TP, TH, Te, TX, T12, To, Tt, Tx, TS;
	  E Tw, TT, TY;
	  {
	       E T3, Tf, TC, TV, T6, Tz, Ti, TU;
	       {
		    E T1, T2, TA, TB;
		    T1 = ri[0];
		    T2 = ri[WS(ris, 7)];
		    T3 = T1 + T2;
		    Tf = T1 - T2;
		    TA = ii[0];
		    TB = ii[WS(iis, 7)];
		    TC = TA + TB;
		    TV = TB - TA;
	       }
	       {
		    E T4, T5, Tg, Th;
		    T4 = ri[WS(ris, 4)];
		    T5 = ri[WS(ris, 3)];
		    T6 = T4 + T5;
		    Tz = T4 - T5;
		    Tg = ii[WS(iis, 4)];
		    Th = ii[WS(iis, 3)];
		    Ti = Tg + Th;
		    TU = Tg - Th;
	       }
	       T7 = T3 + T6;
	       TW = TU + TV;
	       T13 = TV - TU;
	       Tj = Tf - Ti;
	       TD = Tz + TC;
	       TK = Tz - TC;
	       TP = T3 - T6;
	       TH = Tf + Ti;
	  }
	  {
	       E Ta, Tk, Tn, TR, Td, Tp, Ts, TQ;
	       {
		    E T8, T9, Tl, Tm;
		    T8 = ri[WS(ris, 2)];
		    T9 = ri[WS(ris, 5)];
		    Ta = T8 + T9;
		    Tk = T8 - T9;
		    Tl = ii[WS(iis, 2)];
		    Tm = ii[WS(iis, 5)];
		    Tn = Tl + Tm;
		    TR = Tl - Tm;
	       }
	       {
		    E Tb, Tc, Tq, Tr;
		    Tb = ri[WS(ris, 1)];
		    Tc = ri[WS(ris, 6)];
		    Td = Tb + Tc;
		    Tp = Tb - Tc;
		    Tq = ii[WS(iis, 1)];
		    Tr = ii[WS(iis, 6)];
		    Ts = Tq + Tr;
		    TQ = Tr - Tq;
	       }
	       Te = Ta + Td;
	       TX = Ta - Td;
	       T12 = TR + TQ;
	       To = Tk - Tn;
	       Tt = Tp - Ts;
	       Tx = Tp + Ts;
	       TS = TQ - TR;
	       Tw = Tk + Tn;
	  }
	  O[0] = KP2_000000000 * (T7 + Te);
	  O[WS(os, 8)] = KP2_000000000 * (T13 - T12);
	  TT = TP + TS;
	  TY = TW - TX;
	  O[WS(os, 2)] = FMA(KP1_847759065, TT, KP765366864 * TY);
	  O[WS(os, 10)] = FNMS(KP765366864, TT, KP1_847759065 * TY);
	  {
	       E T11, T14, TZ, T10;
	       T11 = T7 - Te;
	       T14 = T12 + T13;
	       O[WS(os, 4)] = KP1_414213562 * (T11 + T14);
	       O[WS(os, 12)] = KP1_414213562 * (T14 - T11);
	       TZ = TP - TS;
	       T10 = TX + TW;
	       O[WS(os, 6)] = FMA(KP765366864, TZ, KP1_847759065 * T10);
	       O[WS(os, 14)] = FNMS(KP1_847759065, TZ, KP765366864 * T10);
	  }
	  {
	       E TJ, TN, TM, TO, TI, TL;
	       TI = KP707106781 * (Tw + Tx);
	       TJ = TH - TI;
	       TN = TH + TI;
	       TL = KP707106781 * (To - Tt);
	       TM = TK - TL;
	       TO = TL + TK;
	       O[WS(os, 3)] = FMA(KP1_662939224, TJ, KP1_111140466 * TM);
	       O[WS(os, 15)] = FNMS(KP1_961570560, TN, KP390180644 * TO);
	       O[WS(os, 11)] = FNMS(KP1_111140466, TJ, KP1_662939224 * TM);
	       O[WS(os, 7)] = FMA(KP390180644, TN, KP1_961570560 * TO);
	  }
	  {
	       E Tv, TF, TE, TG, Tu, Ty;
	       Tu = KP707106781 * (To + Tt);
	       Tv = Tj + Tu;
	       TF = Tj - Tu;
	       Ty = KP707106781 * (Tw - Tx);
	       TE = Ty + TD;
	       TG = Ty - TD;
	       O[WS(os, 1)] = FNMS(KP390180644, TE, KP1_961570560 * Tv);
	       O[WS(os, 13)] = FNMS(KP1_662939224, TF, KP1_111140466 * TG);
	       O[WS(os, 9)] = -(FMA(KP390180644, Tv, KP1_961570560 * TE));
	       O[WS(os, 5)] = FMA(KP1_111140466, TF, KP1_662939224 * TG);
	  }
     }
}

static const khc2r_desc desc = { 16, "hc2rIII_16", {54, 20, 12, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2rIII_16) (planner *p) {
     X(khc2rIII_register) (p, hc2rIII_16, &desc);
}

#endif				/* HAVE_FMA */
