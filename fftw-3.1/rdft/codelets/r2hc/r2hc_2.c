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
/* Generated on Fri Jan 27 20:16:12 EST 2006 */

#include "codelet-rdft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_r2hc -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 2 -name r2hc_2 -include r2hc.h */

/*
 * This function contains 2 FP additions, 0 FP multiplications,
 * (or, 2 additions, 0 multiplications, 0 fused multiply/add),
 * 3 stack variables, and 4 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_r2hc.ml,v 1.17 2006-01-05 03:04:27 stevenj Exp $
 */

#include "r2hc.h"

static void r2hc_2(const R *I, R *ro, R *io, stride is, stride ros, stride ios, INT v, INT ivs, INT ovs)
{
     INT i;
     for (i = v; i > 0; i = i - 1, I = I + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(ros), MAKE_VOLATILE_STRIDE(ios)) {
	  E T1, T2;
	  T1 = I[0];
	  T2 = I[WS(is, 1)];
	  ro[0] = T1 + T2;
	  ro[WS(ros, 1)] = T1 - T2;
     }
}

static const kr2hc_desc desc = { 2, "r2hc_2", {2, 0, 0, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_r2hc_2) (planner *p) {
     X(kr2hc_register) (p, r2hc_2, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_r2hc -compact -variables 4 -pipeline-latency 4 -n 2 -name r2hc_2 -include r2hc.h */

/*
 * This function contains 2 FP additions, 0 FP multiplications,
 * (or, 2 additions, 0 multiplications, 0 fused multiply/add),
 * 3 stack variables, and 4 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_r2hc.ml,v 1.17 2006-01-05 03:04:27 stevenj Exp $
 */

#include "r2hc.h"

static void r2hc_2(const R *I, R *ro, R *io, stride is, stride ros, stride ios, INT v, INT ivs, INT ovs)
{
     INT i;
     for (i = v; i > 0; i = i - 1, I = I + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(ros), MAKE_VOLATILE_STRIDE(ios)) {
	  E T1, T2;
	  T1 = I[0];
	  T2 = I[WS(is, 1)];
	  ro[WS(ros, 1)] = T1 - T2;
	  ro[0] = T1 + T2;
     }
}

static const kr2hc_desc desc = { 2, "r2hc_2", {2, 0, 0, 0}, &GENUS, 0, 0, 0, 0, 0 };

void X(codelet_r2hc_2) (planner *p) {
     X(kr2hc_register) (p, r2hc_2, &desc);
}

#endif				/* HAVE_FMA */
