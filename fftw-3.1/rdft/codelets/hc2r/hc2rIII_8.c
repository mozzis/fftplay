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
/* Generated on Fri Jan 27 20:52:17 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2r -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -name hc2rIII_8 -dft-III -include hc2rIII.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 18 additions, 8 multiplications, 4 fused multiply/add),
 * 23 stack variables, and 16 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2rIII.h"

static void hc2rIII_8(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E T4, T7, T3, Tl, Tf, T5, T8, T9, T6, Tc;
	  {
	       E T1, T2, Td, Te;
	       T1 = ri[0];
	       T2 = ri[WS(ris, 3)];
	       Td = ii[0];
	       Te = ii[WS(iis, 3)];
	       T4 = ri[WS(ris, 2)];
	       T7 = T1 - T2;
	       T3 = T1 + T2;
	       Tl = Te - Td;
	       Tf = Td + Te;
	       T5 = ri[WS(ris, 1)];
	       T8 = ii[WS(iis, 2)];
	       T9 = ii[WS(iis, 1)];
	  }
	  T6 = T4 + T5;
	  Tc = T4 - T5;
	  {
	       E Ta, Tk, Tg, Th;
	       Ta = T8 + T9;
	       Tk = T8 - T9;
	       Tg = Tc + Tf;
	       Th = Tc - Tf;
	       {
		    E Tj, Tm, Tb, Ti;
		    Tj = T3 - T6;
		    O[0] = KP2_000000000 * (T3 + T6);
		    Tm = Tk + Tl;
		    O[WS(os, 4)] = KP2_000000000 * (Tl - Tk);
		    Tb = T7 - Ta;
		    Ti = T7 + Ta;
		    O[WS(os, 6)] = KP1_414213562 * (Tm - Tj);
		    O[WS(os, 2)] = KP1_414213562 * (Tj + Tm);
		    O[WS(os, 7)] = -(KP1_847759065 * (FNMS(KP414213562, Th, Ti)));
		    O[WS(os, 3)] = KP1_847759065 * (FMA(KP414213562, Ti, Th));
		    O[WS(os, 5)] = -(KP1_847759065 * (FMA(KP414213562, Tb, Tg)));
		    O[WS(os, 1)] = KP1_847759065 * (FNMS(KP414213562, Tg, Tb));
	       }
	  }
     }
}

static const khc2r_desc desc = { 8, "hc2rIII_8", {18, 8, 4, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2rIII_8) (planner *p) {
     X(khc2rIII_register) (p, hc2rIII_8, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2r -compact -variables 4 -pipeline-latency 4 -sign 1 -n 8 -name hc2rIII_8 -dft-III -include hc2rIII.h */

/*
 * This function contains 22 FP additions, 12 FP multiplications,
 * (or, 18 additions, 8 multiplications, 4 fused multiply/add),
 * 19 stack variables, and 16 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2rIII.h"

static void hc2rIII_8(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP1_414213562, +1.414213562373095048801688724209698078569671875);
     DK(KP765366864, +0.765366864730179543456919968060797733522689125);
     DK(KP1_847759065, +1.847759065022573512256366378793576573644833252);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E T3, T7, Tf, Tl, T6, Tc, Ta, Tk, Tb, Tg;
	  {
	       E T1, T2, Td, Te;
	       T1 = ri[0];
	       T2 = ri[WS(ris, 3)];
	       T3 = T1 + T2;
	       T7 = T1 - T2;
	       Td = ii[0];
	       Te = ii[WS(iis, 3)];
	       Tf = Td + Te;
	       Tl = Te - Td;
	  }
	  {
	       E T4, T5, T8, T9;
	       T4 = ri[WS(ris, 2)];
	       T5 = ri[WS(ris, 1)];
	       T6 = T4 + T5;
	       Tc = T4 - T5;
	       T8 = ii[WS(iis, 2)];
	       T9 = ii[WS(iis, 1)];
	       Ta = T8 + T9;
	       Tk = T8 - T9;
	  }
	  O[0] = KP2_000000000 * (T3 + T6);
	  O[WS(os, 4)] = KP2_000000000 * (Tl - Tk);
	  Tb = T7 - Ta;
	  Tg = Tc + Tf;
	  O[WS(os, 1)] = FNMS(KP765366864, Tg, KP1_847759065 * Tb);
	  O[WS(os, 5)] = -(FMA(KP765366864, Tb, KP1_847759065 * Tg));
	  {
	       E Th, Ti, Tj, Tm;
	       Th = T7 + Ta;
	       Ti = Tc - Tf;
	       O[WS(os, 3)] = FMA(KP765366864, Th, KP1_847759065 * Ti);
	       O[WS(os, 7)] = FNMS(KP1_847759065, Th, KP765366864 * Ti);
	       Tj = T3 - T6;
	       Tm = Tk + Tl;
	       O[WS(os, 2)] = KP1_414213562 * (Tj + Tm);
	       O[WS(os, 6)] = KP1_414213562 * (Tm - Tj);
	  }
     }
}

static const khc2r_desc desc = { 8, "hc2rIII_8", {18, 8, 4, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2rIII_8) (planner *p) {
     X(khc2rIII_register) (p, hc2rIII_8, &desc);
}

#endif				/* HAVE_FMA */
