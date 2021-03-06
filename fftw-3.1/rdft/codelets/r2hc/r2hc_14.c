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
/* Generated on Fri Jan 27 20:16:32 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_r2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 14 -name r2hc_14 -include r2hc.h */

/*
 * This function contains 62 FP additions, 36 FP multiplications,
 * (or, 32 additions, 6 multiplications, 30 fused multiply/add),
 * 45 stack variables, and 28 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_r2hc.ml,v 1.17 2006-01-05 03:04:27 stevenj Exp $
 */

#include "r2hc.h"

static void r2hc_14(const R *I, R *ro, R *io, stride is, stride ros, stride ios, INT v, INT ivs, INT ovs)
{
     DK(KP900968867, +0.900968867902419126236102319507445051165919162);
     DK(KP692021471, +0.692021471630095869627814897002069140197260599);
     DK(KP801937735, +0.801937735804838252472204639014890102331838324);
     DK(KP974927912, +0.974927912181823607018131682993931217232785801);
     DK(KP356895867, +0.356895867892209443894399510021300583399127187);
     DK(KP554958132, +0.554958132087371191422194871006410481067288862);
     INT i;
     for (i = v; i > 0; i = i - 1, I = I + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(ros), MAKE_VOLATILE_STRIDE(ios)) {
	  E TN, T3, TG, TQ, Tx, To, TH, Td, TD, TO, Tw, Ta, TL, Ty, TT;
	  E TI, Tg, Tr, Te, Tf, TP, TJ;
	  {
	       E Tl, TE, Tk, Tm;
	       {
		    E T1, T2, Ti, Tj;
		    T1 = I[0];
		    T2 = I[WS(is, 7)];
		    Ti = I[WS(is, 6)];
		    Tj = I[WS(is, 13)];
		    Tl = I[WS(is, 8)];
		    TN = T1 + T2;
		    T3 = T1 - T2;
		    TE = Ti + Tj;
		    Tk = Ti - Tj;
		    Tm = I[WS(is, 1)];
	       }
	       {
		    E T7, TC, T6, T8;
		    {
			 E T4, T5, TF, Tn;
			 T4 = I[WS(is, 2)];
			 T5 = I[WS(is, 9)];
			 T7 = I[WS(is, 12)];
			 TF = Tl + Tm;
			 Tn = Tl - Tm;
			 TC = T4 + T5;
			 T6 = T4 - T5;
			 TG = TE - TF;
			 TQ = TE + TF;
			 Tx = Tn - Tk;
			 To = Tk + Tn;
			 T8 = I[WS(is, 5)];
		    }
		    {
			 E Tb, Tc, TB, T9;
			 Tb = I[WS(is, 4)];
			 Tc = I[WS(is, 11)];
			 Te = I[WS(is, 10)];
			 TB = T7 + T8;
			 T9 = T7 - T8;
			 TH = Tb + Tc;
			 Td = Tb - Tc;
			 TD = TB - TC;
			 TO = TC + TB;
			 Tw = T6 - T9;
			 Ta = T6 + T9;
			 Tf = I[WS(is, 3)];
		    }
	       }
	  }
	  TL = FNMS(KP554958132, TG, TD);
	  Ty = FNMS(KP554958132, Tx, Tw);
	  TT = FNMS(KP356895867, TO, TQ);
	  TI = Te + Tf;
	  Tg = Te - Tf;
	  Tr = FNMS(KP356895867, Ta, To);
	  TP = TH + TI;
	  TJ = TH - TI;
	  {
	       E Th, Tv, TK, TM;
	       Th = Td + Tg;
	       Tv = Tg - Td;
	       TK = FMA(KP554958132, TJ, TG);
	       TM = FMA(KP554958132, TD, TJ);
	       io[WS(ios, 6)] = KP974927912 * (FNMS(KP801937735, TL, TJ));
	       {
		    E TR, TV, TU, Tz;
		    TR = FNMS(KP356895867, TQ, TP);
		    TV = FNMS(KP356895867, TP, TO);
		    TU = FNMS(KP692021471, TT, TP);
		    ro[0] = TN + TO + TP + TQ;
		    Tz = FMA(KP554958132, Tv, Tx);
		    io[WS(ios, 1)] = KP974927912 * (FNMS(KP801937735, Ty, Tv));
		    {
			 E TA, Ts, Tt, Tp;
			 TA = FMA(KP554958132, Tw, Tv);
			 Ts = FNMS(KP692021471, Tr, Th);
			 Tt = FNMS(KP356895867, Th, Ta);
			 Tp = FNMS(KP356895867, To, Th);
			 ro[WS(ros, 7)] = T3 + Ta + Th + To;
			 io[WS(ios, 2)] = KP974927912 * (FMA(KP801937735, TK, TD));
			 io[WS(ios, 4)] = KP974927912 * (FNMS(KP801937735, TM, TG));
			 {
			      E TS, TW, Tu, Tq;
			      TS = FNMS(KP692021471, TR, TO);
			      TW = FNMS(KP692021471, TV, TQ);
			      ro[WS(ros, 2)] = FNMS(KP900968867, TU, TN);
			      io[WS(ios, 5)] = KP974927912 * (FMA(KP801937735, Tz, Tw));
			      io[WS(ios, 3)] = KP974927912 * (FNMS(KP801937735, TA, Tx));
			      ro[WS(ros, 5)] = FNMS(KP900968867, Ts, T3);
			      Tu = FNMS(KP692021471, Tt, To);
			      Tq = FNMS(KP692021471, Tp, Ta);
			      ro[WS(ros, 4)] = FNMS(KP900968867, TS, TN);
			      ro[WS(ros, 6)] = FNMS(KP900968867, TW, TN);
			      ro[WS(ros, 1)] = FNMS(KP900968867, Tu, T3);
			      ro[WS(ros, 3)] = FNMS(KP900968867, Tq, T3);
			 }
		    }
	       }
	  }
     }
}

static const kr2hc_desc desc = { 14, "r2hc_14", {32, 6, 30, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_r2hc_14) (planner *p) {
     X(kr2hc_register) (p, r2hc_14, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_r2hc -compact -variables 4 -pipeline-latency 4 -n 14 -name r2hc_14 -include r2hc.h */

/*
 * This function contains 62 FP additions, 36 FP multiplications,
 * (or, 38 additions, 12 multiplications, 24 fused multiply/add),
 * 29 stack variables, and 28 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_r2hc.ml,v 1.17 2006-01-05 03:04:27 stevenj Exp $
 */

#include "r2hc.h"

static void r2hc_14(const R *I, R *ro, R *io, stride is, stride ros, stride ios, INT v, INT ivs, INT ovs)
{
     DK(KP900968867, +0.900968867902419126236102319507445051165919162);
     DK(KP222520933, +0.222520933956314404288902564496794759466355569);
     DK(KP623489801, +0.623489801858733530525004884004239810632274731);
     DK(KP433883739, +0.433883739117558120475768332848358754609990728);
     DK(KP974927912, +0.974927912181823607018131682993931217232785801);
     DK(KP781831482, +0.781831482468029808708444526674057750232334519);
     INT i;
     for (i = v; i > 0; i = i - 1, I = I + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(ros), MAKE_VOLATILE_STRIDE(ios)) {
	  E T3, TB, T6, Tv, Tn, Ts, Tk, Tt, Td, Ty, T9, Tw, Tg, Tz, T1;
	  E T2;
	  T1 = I[0];
	  T2 = I[WS(is, 7)];
	  T3 = T1 - T2;
	  TB = T1 + T2;
	  {
	       E T4, T5, Tl, Tm;
	       T4 = I[WS(is, 4)];
	       T5 = I[WS(is, 11)];
	       T6 = T4 - T5;
	       Tv = T4 + T5;
	       Tl = I[WS(is, 12)];
	       Tm = I[WS(is, 5)];
	       Tn = Tl - Tm;
	       Ts = Tl + Tm;
	  }
	  {
	       E Ti, Tj, Tb, Tc;
	       Ti = I[WS(is, 2)];
	       Tj = I[WS(is, 9)];
	       Tk = Ti - Tj;
	       Tt = Ti + Tj;
	       Tb = I[WS(is, 6)];
	       Tc = I[WS(is, 13)];
	       Td = Tb - Tc;
	       Ty = Tb + Tc;
	  }
	  {
	       E T7, T8, Te, Tf;
	       T7 = I[WS(is, 10)];
	       T8 = I[WS(is, 3)];
	       T9 = T7 - T8;
	       Tw = T7 + T8;
	       Te = I[WS(is, 8)];
	       Tf = I[WS(is, 1)];
	       Tg = Te - Tf;
	       Tz = Te + Tf;
	  }
	  {
	       E Tp, Tr, Tq, Ta, To, Th;
	       Tp = Tn - Tk;
	       Tr = Tg - Td;
	       Tq = T9 - T6;
	       io[WS(ios, 1)] = FMA(KP781831482, Tp, KP974927912 * Tq) + (KP433883739 * Tr);
	       io[WS(ios, 5)] = FMA(KP433883739, Tq, KP781831482 * Tr) - (KP974927912 * Tp);
	       io[WS(ios, 3)] = FMA(KP433883739, Tp, KP974927912 * Tr) - (KP781831482 * Tq);
	       Ta = T6 + T9;
	       To = Tk + Tn;
	       Th = Td + Tg;
	       ro[WS(ros, 3)] = FMA(KP623489801, Ta, T3) + FNMA(KP222520933, Th, KP900968867 * To);
	       ro[WS(ros, 7)] = T3 + To + Ta + Th;
	       ro[WS(ros, 1)] = FMA(KP623489801, To, T3) + FNMA(KP900968867, Th, KP222520933 * Ta);
	       ro[WS(ros, 5)] = FMA(KP623489801, Th, T3) + FNMA(KP900968867, Ta, KP222520933 * To);
	  }
	  {
	       E Tu, TA, Tx, TC, TE, TD;
	       Tu = Ts - Tt;
	       TA = Ty - Tz;
	       Tx = Tv - Tw;
	       io[WS(ios, 2)] = FMA(KP974927912, Tu, KP433883739 * Tx) + (KP781831482 * TA);
	       io[WS(ios, 6)] = FMA(KP974927912, Tx, KP433883739 * TA) - (KP781831482 * Tu);
	       io[WS(ios, 4)] = FNMS(KP781831482, Tx, KP974927912 * TA) - (KP433883739 * Tu);
	       TC = Tt + Ts;
	       TE = Tv + Tw;
	       TD = Ty + Tz;
	       ro[WS(ros, 6)] = FMA(KP623489801, TC, TB) + FNMA(KP900968867, TD, KP222520933 * TE);
	       ro[WS(ros, 2)] = FMA(KP623489801, TD, TB) + FNMA(KP900968867, TE, KP222520933 * TC);
	       ro[WS(ros, 4)] = FMA(KP623489801, TE, TB) + FNMA(KP222520933, TD, KP900968867 * TC);
	       ro[0] = TB + TC + TE + TD;
	  }
     }
}

static const kr2hc_desc desc = { 14, "r2hc_14", {38, 12, 24, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_r2hc_14) (planner *p) {
     X(kr2hc_register) (p, r2hc_14, &desc);
}

#endif				/* HAVE_FMA */
