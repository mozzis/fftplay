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
/* Generated on Fri Jan 27 20:39:05 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_hc2r -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -sign 1 -n 3 -name hc2r_3 -include hc2r.h */

/*
 * This function contains 4 FP additions, 3 FP multiplications,
 * (or, 1 additions, 0 multiplications, 3 fused multiply/add),
 * 7 stack variables, and 6 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2r.h"

static void hc2r_3(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP1_732050807, +1.732050807568877293527446341505872366942805254);
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E T4, T1, T2, T3;
	  T4 = ii[WS(iis, 1)];
	  T1 = ri[0];
	  T2 = ri[WS(ris, 1)];
	  O[0] = FMA(KP2_000000000, T2, T1);
	  T3 = T1 - T2;
	  O[WS(os, 1)] = FNMS(KP1_732050807, T4, T3);
	  O[WS(os, 2)] = FMA(KP1_732050807, T4, T3);
     }
}

static const khc2r_desc desc = { 3, "hc2r_3", {1, 0, 3, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2r_3) (planner *p) {
     X(khc2r_register) (p, hc2r_3, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_hc2r -compact -variables 4 -pipeline-latency 4 -sign 1 -n 3 -name hc2r_3 -include hc2r.h */

/*
 * This function contains 4 FP additions, 2 FP multiplications,
 * (or, 3 additions, 1 multiplications, 1 fused multiply/add),
 * 8 stack variables, and 6 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_hc2r.ml,v 1.18 2006-01-05 03:04:27 stevenj Exp $
 */

#include "hc2r.h"

static void hc2r_3(const R *ri, const R *ii, R *O, stride ris, stride iis, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP2_000000000, +2.000000000000000000000000000000000000000000000);
     DK(KP1_732050807, +1.732050807568877293527446341505872366942805254);
     INT i;
     for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, O = O + ovs, MAKE_VOLATILE_STRIDE(ris), MAKE_VOLATILE_STRIDE(iis), MAKE_VOLATILE_STRIDE(os)) {
	  E T5, T1, T2, T3, T4;
	  T4 = ii[WS(iis, 1)];
	  T5 = KP1_732050807 * T4;
	  T1 = ri[0];
	  T2 = ri[WS(ris, 1)];
	  T3 = T1 - T2;
	  O[0] = FMA(KP2_000000000, T2, T1);
	  O[WS(os, 2)] = T3 + T5;
	  O[WS(os, 1)] = T3 - T5;
     }
}

static const khc2r_desc desc = { 3, "hc2r_3", {3, 1, 1, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_hc2r_3) (planner *p) {
     X(khc2r_register) (p, hc2r_3, &desc);
}

#endif				/* HAVE_FMA */
