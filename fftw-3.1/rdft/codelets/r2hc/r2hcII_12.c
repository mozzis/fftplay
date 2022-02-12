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
/* Generated on Fri Jan 27 20:31:48 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_r2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 12 -name r2hcII_12 -dft-II -include r2hcII.h */

/*
 * This function contains 45 FP additions, 24 FP multiplications,
 * (or, 21 additions, 0 multiplications, 24 fused multiply/add),
 * 37 stack variables, and 24 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_r2hc.ml,v 1.17 2006-01-05 03:04:27 stevenj Exp $
 */

#include "r2hcII.h"

static void r2hcII_12(const R *I, R *ro, R *io, stride is, stride ros, stride ios, INT v, INT ivs, INT ovs)
{
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, I = I + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(ros), MAKE_VOLATILE_STRIDE(ios)) {
	  E TD, TB, Tp, T9, Tq, Tr, TE, To, Ts, TC;
	  {
	       E T8, T1, Tv, Tm, TF, Tz, Tl, Ta, Tb, Tt, TA, T4, Tc;
	       {
		    E Tx, Th, Ti, Tj, Ty, T6, T7, T2, T3, Tk;
		    Tx = I[WS(is, 6)];
		    T6 = I[WS(is, 10)];
		    T7 = I[WS(is, 2)];
		    Th = I[WS(is, 9)];
		    Ti = I[WS(is, 5)];
		    Tj = I[WS(is, 1)];
		    Ty = T6 + T7;
		    T8 = T6 - T7;
		    T1 = I[0];
		    Tv = Ti - Tj - Th;
		    Tk = Ti - Tj;
		    Tm = Ti + Tj;
		    TF = Tx - Ty;
		    Tz = FMA(KP500000000, Ty, Tx);
		    T2 = I[WS(is, 4)];
		    T3 = I[WS(is, 8)];
		    Tl = FMA(KP500000000, Tk, Th);
		    Ta = I[WS(is, 3)];
		    Tb = I[WS(is, 7)];
		    Tt = T1 + T3 - T2;
		    TA = T3 + T2;
		    T4 = T2 - T3;
		    Tc = I[WS(is, 11)];
	       }
	       {
		    E Tn, Tg, T5, Tu;
		    TD = FNMS(KP866025403, TA, Tz);
		    TB = FMA(KP866025403, TA, Tz);
		    T5 = FMA(KP500000000, T4, T1);
		    Tu = Ta + Tc - Tb;
		    {
			 E Td, Tf, TG, Tw, Te;
			 Td = Tb - Tc;
			 Tf = Tc + Tb;
			 Tp = FMA(KP866025403, T8, T5);
			 T9 = FNMS(KP866025403, T8, T5);
			 TG = Tv - Tu;
			 Tw = Tu + Tv;
			 Te = FMA(KP500000000, Td, Ta);
			 Tq = FMA(KP866025403, Tm, Tl);
			 Tn = FNMS(KP866025403, Tm, Tl);
			 io[WS(ios, 1)] = FMA(KP707106781, TG, TF);
			 io[WS(ios, 4)] = FMS(KP707106781, TG, TF);
			 ro[WS(ros, 4)] = FMA(KP707106781, Tw, Tt);
			 ro[WS(ros, 1)] = FNMS(KP707106781, Tw, Tt);
			 Tg = FNMS(KP866025403, Tf, Te);
			 Tr = FMA(KP866025403, Tf, Te);
		    }
		    TE = Tg + Tn;
		    To = Tg - Tn;
	       }
	  }
	  io[WS(ios, 2)] = FMS(KP707106781, TE, TD);
	  io[WS(ios, 3)] = FMA(KP707106781, TE, TD);
	  ro[0] = FMA(KP707106781, To, T9);
	  ro[WS(ros, 5)] = FNMS(KP707106781, To, T9);
	  Ts = Tq - Tr;
	  TC = Tr + Tq;
	  io[0] = -(FMA(KP707106781, TC, TB));
	  io[WS(ios, 5)] = FNMS(KP707106781, TC, TB);
	  ro[WS(ros, 2)] = FMA(KP707106781, Ts, Tp);
	  ro[WS(ros, 3)] = FNMS(KP707106781, Ts, Tp);
     }
}

static const kr2hc_desc desc = { 12, "r2hcII_12", {21, 0, 24, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_r2hcII_12) (planner *p) {
     X(kr2hcII_register) (p, r2hcII_12, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_r2hc -compact -variables 4 -pipeline-latency 4 -n 12 -name r2hcII_12 -dft-II -include r2hcII.h */

/*
 * This function contains 43 FP additions, 12 FP multiplications,
 * (or, 39 additions, 8 multiplications, 4 fused multiply/add),
 * 28 stack variables, and 24 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_r2hc.ml,v 1.17 2006-01-05 03:04:27 stevenj Exp $
 */

#include "r2hcII.h"

static void r2hcII_12(const R *I, R *ro, R *io, stride is, stride ros, stride ios, INT v, INT ivs, INT ovs)
{
     DK(KP353553390, +0.353553390593273762200422181052424519642417969);
     DK(KP707106781, +0.707106781186547524400844362104849039284835938);
     DK(KP612372435, +0.612372435695794524549321018676472847991486870);
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DK(KP866025403, +0.866025403784438646763723170752936183471402627);
     INT i;
     for (i = v; i > 0; i = i - 1, I = I + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(ros), MAKE_VOLATILE_STRIDE(ios)) {
	  E Tx, Tg, T4, Tz, Ty, Tj, TA, T9, Tm, Tl, Te, Tp, To, Tf, TE;
	  E TF;
	  {
	       E T1, T3, T2, Th, Ti;
	       T1 = I[0];
	       T3 = I[WS(is, 4)];
	       T2 = I[WS(is, 8)];
	       Tx = KP866025403 * (T2 + T3);
	       Tg = FMA(KP500000000, T3 - T2, T1);
	       T4 = T1 + T2 - T3;
	       Tz = I[WS(is, 6)];
	       Th = I[WS(is, 10)];
	       Ti = I[WS(is, 2)];
	       Ty = Th + Ti;
	       Tj = KP866025403 * (Th - Ti);
	       TA = FMA(KP500000000, Ty, Tz);
	  }
	  {
	       E T5, T6, T7, T8;
	       T5 = I[WS(is, 3)];
	       T6 = I[WS(is, 11)];
	       T7 = I[WS(is, 7)];
	       T8 = T6 - T7;
	       T9 = T5 + T8;
	       Tm = KP612372435 * (T6 + T7);
	       Tl = FNMS(KP353553390, T8, KP707106781 * T5);
	  }
	  {
	       E Td, Ta, Tb, Tc;
	       Td = I[WS(is, 9)];
	       Ta = I[WS(is, 5)];
	       Tb = I[WS(is, 1)];
	       Tc = Ta - Tb;
	       Te = Tc - Td;
	       Tp = FMA(KP353553390, Tc, KP707106781 * Td);
	       To = KP612372435 * (Ta + Tb);
	  }
	  Tf = KP707106781 * (T9 + Te);
	  ro[WS(ros, 1)] = T4 - Tf;
	  ro[WS(ros, 4)] = T4 + Tf;
	  TE = KP707106781 * (Te - T9);
	  TF = Tz - Ty;
	  io[WS(ios, 4)] = TE - TF;
	  io[WS(ios, 1)] = TE + TF;
	  {
	       E Tk, TB, Tr, Tw, Tn, Tq;
	       Tk = Tg - Tj;
	       TB = Tx - TA;
	       Tn = Tl - Tm;
	       Tq = To - Tp;
	       Tr = Tn + Tq;
	       Tw = Tn - Tq;
	       ro[WS(ros, 5)] = Tk - Tr;
	       io[WS(ios, 2)] = Tw + TB;
	       ro[0] = Tk + Tr;
	       io[WS(ios, 3)] = Tw - TB;
	  }
	  {
	       E Ts, TD, Tv, TC, Tt, Tu;
	       Ts = Tg + Tj;
	       TD = Tx + TA;
	       Tt = To + Tp;
	       Tu = Tm + Tl;
	       Tv = Tt - Tu;
	       TC = Tu + Tt;
	       ro[WS(ros, 3)] = Ts - Tv;
	       io[WS(ios, 5)] = TD - TC;
	       ro[WS(ros, 2)] = Ts + Tv;
	       io[0] = -(TC + TD);
	  }
     }
}

static const kr2hc_desc desc = { 12, "r2hcII_12", {39, 8, 4, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_r2hcII_12) (planner *p) {
     X(kr2hcII_register) (p, r2hcII_12, &desc);
}

#endif				/* HAVE_FMA */