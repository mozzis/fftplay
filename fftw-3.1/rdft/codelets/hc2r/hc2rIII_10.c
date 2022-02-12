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
/* Generated on Fri Jan 27 20:52:26 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2r -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 10 -name hc2rIII_10 -dft-III -include hc2rIII.h */

/*
 * This function contains 32 FP additions, 28 FP multiplications,
 * (or, 14 additions, 10 multiplications, 18 fused multiply/add),
 * 38 stack variables, and 20 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2rIII.h"

static void hc2rIII_10(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DK(KP618033988, +0.618033988749894848204586834365638117720309180);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E Tq, Ti, Tk, Tu, Tw, Tp, Tb, Tj, Tr, Tv;
	  {
	       E T1, To, Ts, Tt, T8, Ta, Te, Tl, Tm, Th, Tn, T9;
	       T1 = ri[WS(ris, 2)];
	       To = ii[WS(iis, 2)];
	       {
		    E T2, T3, T5, T6;
		    T2 = ri[WS(ris, 4)];
		    T3 = ri[0];
		    T5 = ri[WS(ris, 3)];
		    T6 = ri[WS(ris, 1)];
		    {
			 E Tc, T4, T7, Td, Tf, Tg;
			 Tc = ii[WS(iis, 3)];
			 Ts = T2 - T3;
			 T4 = T2 + T3;
			 Tt = T5 - T6;
			 T7 = T5 + T6;
			 Td = ii[WS(iis, 1)];
			 Tf = ii[WS(iis, 4)];
			 Tg = ii[0];
			 T8 = T4 + T7;
			 Ta = T7 - T4;
			 Te = Tc - Td;
			 Tl = Tc + Td;
			 Tm = Tf + Tg;
			 Th = Tf - Tg;
		    }
	       }
	       O[0] = KP2_000000000 * (T1 + T8);
	       Tn = Tl - Tm;
	       Tq = Tl + Tm;
	       Ti = FMA(KP618033988, Th, Te);
	       Tk = FNMS(KP618033988, Te, Th);
	       O[WS(os, 5)] = KP2_000000000 * (Tn - To);
	       T9 = FMS(KP250000000, T8, T1);
	       Tu = FMA(KP618033988, Tt, Ts);
	       Tw = FNMS(KP618033988, Ts, Tt);
	       Tp = FMA(KP250000000, Tn, To);
	       Tb = FNMS(KP559016994, Ta, T9);
	       Tj = FMA(KP559016994, Ta, T9);
	  }
	  Tr = FMA(KP559016994, Tq, Tp);
	  Tv = FNMS(KP559016994, Tq, Tp);
	  O[WS(os, 4)] = -(KP2_000000000 * (FNMS(KP951056516, Tk, Tj)));
	  O[WS(os, 6)] = KP2_000000000 * (FMA(KP951056516, Tk, Tj));
	  O[WS(os, 8)] = -(KP2_000000000 * (FNMS(KP951056516, Ti, Tb)));
	  O[WS(os, 2)] = KP2_000000000 * (FMA(KP951056516, Ti, Tb));
	  O[WS(os, 3)] = KP2_000000000 * (FMA(KP951056516, Tw, Tv));
	  O[WS(os, 7)] = KP2_000000000 * (FNMS(KP951056516, Tw, Tv));
	  O[WS(os, 9)] = -(KP2_000000000 * (FNMS(KP951056516, Tu, Tr)));
	  O[WS(os, 1)] = -(KP2_000000000 * (FMA(KP951056516, Tu, Tr)));
     }
}

static const khc2r_desc desc = { 10, "hc2rIII_10", {14, 10, 18, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2rIII_10) (planner *p) {
     X(khc2rIII_register) (p, hc2rIII_10, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2r -compact -variables 4 -pipeline-latency 4 -sign 1 -n 10 -name hc2rIII_10 -dft-III -include hc2rIII.h */

/*
 * This function contains 32 FP additions, 16 FP multiplications,
 * (or, 26 additions, 10 multiplications, 6 fused multiply/add),
 * 22 stack variables, and 20 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2rIII.h"

static void hc2rIII_10(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP500000000, +0.500000000000000000000000000000000000000000000);
     DK(KP1_902113032, +1.902113032590307144232878666758764286811397268);
     DK(KP1_175570504, +1.175570504584946258337411909278145537195304875);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     DK(KP1_118033988, +1.118033988749894848204586834365638117720309180);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E T1, To, T8, Tq, Ta, Tp, Te, Ts, Th, Tn;
	  T1 = ri[WS(ris, 2)];
	  To = ii[WS(iis, 2)];
	  {
	       E T2, T3, T4, T5, T6, T7;
	       T2 = ri[WS(ris, 4)];
	       T3 = ri[0];
	       T4 = T2 + T3;
	       T5 = ri[WS(ris, 3)];
	       T6 = ri[WS(ris, 1)];
	       T7 = T5 + T6;
	       T8 = T4 + T7;
	       Tq = T5 - T6;
	       Ta = KP1_118033988 * (T7 - T4);
	       Tp = T2 - T3;
	  }
	  {
	       E Tc, Td, Tm, Tf, Tg, Tl;
	       Tc = ii[WS(iis, 4)];
	       Td = ii[0];
	       Tm = Tc + Td;
	       Tf = ii[WS(iis, 1)];
	       Tg = ii[WS(iis, 3)];
	       Tl = Tg + Tf;
	       Te = Tc - Td;
	       Ts = KP1_118033988 * (Tl + Tm);
	       Th = Tf - Tg;
	       Tn = Tl - Tm;
	  }
	  O[0] = KP2_000000000 * (T1 + T8);
	  O[WS(os, 5)] = KP2_000000000 * (Tn - To);
	  {
	       E Ti, Tj, Tb, Tk, T9;
	       Ti = FNMS(KP1_902113032, Th, KP1_175570504 * Te);
	       Tj = FMA(KP1_175570504, Th, KP1_902113032 * Te);
	       T9 = FNMS(KP2_000000000, T1, KP500000000 * T8);
	       Tb = T9 - Ta;
	       Tk = T9 + Ta;
	       O[WS(os, 2)] = Tb + Ti;
	       O[WS(os, 6)] = Tk + Tj;
	       O[WS(os, 8)] = Ti - Tb;
	       O[WS(os, 4)] = Tj - Tk;
	  }
	  {
	       E Tr, Tv, Tu, Tw, Tt;
	       Tr = FMA(KP1_902113032, Tp, KP1_175570504 * Tq);
	       Tv = FNMS(KP1_175570504, Tp, KP1_902113032 * Tq);
	       Tt = FMA(KP500000000, Tn, KP2_000000000 * To);
	       Tu = Ts + Tt;
	       Tw = Tt - Ts;
	       O[WS(os, 1)] = -(Tr + Tu);
	       O[WS(os, 7)] = Tw - Tv;
	       O[WS(os, 9)] = Tr - Tu;
	       O[WS(os, 3)] = Tv + Tw;
	  }
     }
}

static const khc2r_desc desc = { 10, "hc2rIII_10", {26, 10, 6, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2rIII_10) (planner *p) {
     X(khc2rIII_register) (p, hc2rIII_10, &desc);
}

#endif				/* HAVE_FMA */
