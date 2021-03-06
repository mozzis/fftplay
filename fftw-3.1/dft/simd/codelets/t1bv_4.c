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
/* Generated on Fri Jan 27 19:59:27 EST 2006 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_twiddle_c -fma -reorder-insns -schedule-for-pipeline -simd -compact -variables 4 -pipeline-latency 8 -n 4 -name t1bv_4 -include t1b.h -sign 1 */

/*
 * This function contains 11 FP additions, 8 FP multiplications,
 * (or, 9 additions, 6 multiplications, 2 fused multiply/add),
 * 13 stack variables, and 8 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle_c.ml,v 1.13 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t1b.h"

static const R *t1bv_4(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     INT i;
     R *x;
     x = ii;
     for (i = m; i > 0; i = i - VL, x = x + (VL * dist), W = W + (TWVL * 6), MAKE_VOLATILE_STRIDE(ios)) {
	  V T1, T7, T2, T5, T8, T3, T6;
	  T1 = LD(&(x[0]), dist, &(x[0]));
	  T7 = LD(&(x[WS(ios, 3)]), dist, &(x[WS(ios, 1)]));
	  T2 = LD(&(x[WS(ios, 2)]), dist, &(x[0]));
	  T5 = LD(&(x[WS(ios, 1)]), dist, &(x[WS(ios, 1)]));
	  T8 = BYTW(&(W[TWVL * 4]), T7);
	  T3 = BYTW(&(W[TWVL * 2]), T2);
	  T6 = BYTW(&(W[0]), T5);
	  {
	       V Ta, T4, Tb, T9;
	       Ta = VADD(T1, T3);
	       T4 = VSUB(T1, T3);
	       Tb = VADD(T6, T8);
	       T9 = VSUB(T6, T8);
	       ST(&(x[0]), VADD(Ta, Tb), dist, &(x[0]));
	       ST(&(x[WS(ios, 2)]), VSUB(Ta, Tb), dist, &(x[0]));
	       ST(&(x[WS(ios, 1)]), VFMAI(T9, T4), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 3)]), VFNMSI(T9, T4), dist, &(x[WS(ios, 1)]));
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     VTW(1),
     VTW(2),
     VTW(3),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 4, "t1bv_4", twinstr, &GENUS, {9, 6, 2, 0}, 0, 0, 0 };

void X(codelet_t1bv_4) (planner *p) {
     X(kdft_dit_register) (p, t1bv_4, &desc);
}
#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_twiddle_c -simd -compact -variables 4 -pipeline-latency 8 -n 4 -name t1bv_4 -include t1b.h -sign 1 */

/*
 * This function contains 11 FP additions, 6 FP multiplications,
 * (or, 11 additions, 6 multiplications, 0 fused multiply/add),
 * 13 stack variables, and 8 memory accesses
 */
/*
 * Generator Id's : 
 * $Id: algsimp.ml,v 1.8 2006-01-05 03:04:27 stevenj Exp $
 * $Id: fft.ml,v 1.4 2006-01-05 03:04:27 stevenj Exp $
 * $Id: gen_twiddle_c.ml,v 1.13 2006-01-05 03:04:27 stevenj Exp $
 */

#include "t1b.h"

static const R *t1bv_4(R *ri, R *ii, const R *W, stride ios, INT m, INT dist)
{
     INT i;
     R *x;
     x = ii;
     for (i = m; i > 0; i = i - VL, x = x + (VL * dist), W = W + (TWVL * 6), MAKE_VOLATILE_STRIDE(ios)) {
	  V T1, T8, T3, T6, T7, T2, T5;
	  T1 = LD(&(x[0]), dist, &(x[0]));
	  T7 = LD(&(x[WS(ios, 3)]), dist, &(x[WS(ios, 1)]));
	  T8 = BYTW(&(W[TWVL * 4]), T7);
	  T2 = LD(&(x[WS(ios, 2)]), dist, &(x[0]));
	  T3 = BYTW(&(W[TWVL * 2]), T2);
	  T5 = LD(&(x[WS(ios, 1)]), dist, &(x[WS(ios, 1)]));
	  T6 = BYTW(&(W[0]), T5);
	  {
	       V T4, T9, Ta, Tb;
	       T4 = VSUB(T1, T3);
	       T9 = VBYI(VSUB(T6, T8));
	       ST(&(x[WS(ios, 3)]), VSUB(T4, T9), dist, &(x[WS(ios, 1)]));
	       ST(&(x[WS(ios, 1)]), VADD(T4, T9), dist, &(x[WS(ios, 1)]));
	       Ta = VADD(T1, T3);
	       Tb = VADD(T6, T8);
	       ST(&(x[WS(ios, 2)]), VSUB(Ta, Tb), dist, &(x[0]));
	       ST(&(x[0]), VADD(Ta, Tb), dist, &(x[0]));
	  }
     }
     return W;
}

static const tw_instr twinstr[] = {
     VTW(1),
     VTW(2),
     VTW(3),
     {TW_NEXT, VL, 0}
};

static const ct_desc desc = { 4, "t1bv_4", twinstr, &GENUS, {11, 6, 0, 0}, 0, 0, 0 };

void X(codelet_t1bv_4) (planner *p) {
     X(kdft_dit_register) (p, t1bv_4, &desc);
}
#endif				/* HAVE_FMA */
