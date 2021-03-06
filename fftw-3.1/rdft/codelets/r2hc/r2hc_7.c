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
/* Generated on Fri Jan 27 20:16:13 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_r2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 7 -name r2hc_7 -include r2hc.h */

/*
 * This function contains 24 FP additions, 18 FP multiplications,
 * (or, 9 additions, 3 multiplications, 15 fused multiply/add),
 * 25 stack variables, and 14 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_r2hc.ml,v 1.17 2006-01-05 03:04:27 stevenj Exp $
 */

#include "r2hc.h"

static void r2hc_7(const R *I, R *ro, R *io, stride is, stride ros, stride ios, INT v, INT ivs, INT ovs)
{
     DK(KP900968867, +0.900968867902419126236102319507445051165919162);
     DK(KP801937735, +0.801937735804838252472204639014890102331838324);
     DK(KP974927912, +0.974927912181823607018131682993931217232785801);
     DK(KP692021471, +0.692021471630095869627814897002069140197260599);
     DK(KP554958132, +0.554958132087371191422194871006410481067288862);
     DK(KP356895867, +0.356895867892209443894399510021300583399127187);
     INT i;
     for (i = v; i > 0; i = i - 1, I = I + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(ros), MAKE_VOLATILE_STRIDE(ios)) {
	  E T1, Tg, Tc;
	  {
	       E Th, T4, Ti, Ta, Tj, T7, Td, T5, T6, Tl, Tk;
	       T1 = I[0];
	       {
		    E T2, T3, T8, T9;
		    T2 = I[WS(is, 1)];
		    T3 = I[WS(is, 6)];
		    T8 = I[WS(is, 3)];
		    T9 = I[WS(is, 4)];
		    T5 = I[WS(is, 2)];
		    Th = T3 - T2;
		    T4 = T2 + T3;
		    T6 = I[WS(is, 5)];
		    Ti = T9 - T8;
		    Ta = T8 + T9;
	       }
	       Tj = T6 - T5;
	       T7 = T5 + T6;
	       Td = FNMS(KP356895867, T4, Ta);
	       Tl = FMA(KP554958132, Ti, Th);
	       Tk = FMA(KP554958132, Tj, Ti);
	       {
		    E Tm, Tf, Tb, Te;
		    Tm = FNMS(KP554958132, Th, Tj);
		    ro[0] = T1 + T4 + T7 + Ta;
		    Tf = FNMS(KP356895867, T7, T4);
		    Tb = FNMS(KP356895867, Ta, T7);
		    Te = FNMS(KP692021471, Td, T7);
		    io[WS(ios, 2)] = KP974927912 * (FNMS(KP801937735, Tk, Th));
		    io[WS(ios, 3)] = KP974927912 * (FNMS(KP801937735, Tm, Ti));
		    Tg = FNMS(KP692021471, Tf, Ta);
		    Tc = FNMS(KP692021471, Tb, T4);
		    ro[WS(ros, 2)] = FNMS(KP900968867, Te, T1);
		    io[WS(ios, 1)] = KP974927912 * (FMA(KP801937735, Tl, Tj));
	       }
	  }
	  ro[WS(ros, 1)] = FNMS(KP900968867, Tg, T1);
	  ro[WS(ros, 3)] = FNMS(KP900968867, Tc, T1);
     }
}

static const kr2hc_desc desc = { 7, "r2hc_7", {9, 3, 15, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_r2hc_7) (planner *p) {
     X(kr2hc_register) (p, r2hc_7, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_r2hc -compact -variables 4 -pipeline-latency 4 -n 7 -name r2hc_7 -include r2hc.h */

/*
 * This function contains 24 FP additions, 18 FP multiplications,
 * (or, 12 additions, 6 multiplications, 12 fused multiply/add),
 * 20 stack variables, and 14 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_r2hc.ml,v 1.17 2006-01-05 03:04:27 stevenj Exp $
 */

#include "r2hc.h"

static void r2hc_7(const R *I, R *ro, R *io, stride is, stride ros, stride ios, INT v, INT ivs, INT ovs)
{
     DK(KP222520933, +0.222520933956314404288902564496794759466355569);
     DK(KP900968867, +0.900968867902419126236102319507445051165919162);
     DK(KP623489801, +0.623489801858733530525004884004239810632274731);
     DK(KP433883739, +0.433883739117558120475768332848358754609990728);
     DK(KP781831482, +0.781831482468029808708444526674057750232334519);
     DK(KP974927912, +0.974927912181823607018131682993931217232785801);
     INT i;
     for (i = v; i > 0; i = i - 1, I = I + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(ros), MAKE_VOLATILE_STRIDE(ios)) {
	  E T1, Ta, Tb, T4, Td, T7, Tc, T8, T9;
	  T1 = I[0];
	  T8 = I[WS(is, 1)];
	  T9 = I[WS(is, 6)];
	  Ta = T8 + T9;
	  Tb = T9 - T8;
	  {
	       E T2, T3, T5, T6;
	       T2 = I[WS(is, 2)];
	       T3 = I[WS(is, 5)];
	       T4 = T2 + T3;
	       Td = T3 - T2;
	       T5 = I[WS(is, 3)];
	       T6 = I[WS(is, 4)];
	       T7 = T5 + T6;
	       Tc = T6 - T5;
	  }
	  io[WS(ios, 2)] = FNMS(KP781831482, Tc, KP974927912 * Tb) - (KP433883739 * Td);
	  io[WS(ios, 1)] = FMA(KP781831482, Tb, KP974927912 * Td) + (KP433883739 * Tc);
	  ro[WS(ros, 2)] = FMA(KP623489801, T7, T1) + FNMA(KP900968867, T4, KP222520933 * Ta);
	  io[WS(ios, 3)] = FMA(KP433883739, Tb, KP974927912 * Tc) - (KP781831482 * Td);
	  ro[WS(ros, 3)] = FMA(KP623489801, T4, T1) + FNMA(KP222520933, T7, KP900968867 * Ta);
	  ro[WS(ros, 1)] = FMA(KP623489801, Ta, T1) + FNMA(KP900968867, T7, KP222520933 * T4);
	  ro[0] = T1 + Ta + T4 + T7;
     }
}

static const kr2hc_desc desc = { 7, "r2hc_7", {12, 6, 12, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_r2hc_7) (planner *p) {
     X(kr2hc_register) (p, r2hc_7, &desc);
}

#endif				/* HAVE_FMA */
